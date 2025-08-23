# streamlit_app.py
# -------------------------------------------------------------
# Gene2Therapy
# Gene list → KEGG enrichment (counts-only) → Disease links (OpenTargets)
# → Drug repurposing + PK/PD + approval status → Visualizations
# -------------------------------------------------------------

import time
import requests
import pandas as pd
import streamlit as st
from collections import defaultdict, Counter
from xml.etree import ElementTree as ET
from Bio import Entrez
import plotly.express as px
import plotly.graph_objects as go
import networkx as nx
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# ----------------------------
# App Config / Theming (must be FIRST Streamlit call)
# ----------------------------
st.set_page_config(page_title="Gene2Therapy – BioContext", page_icon="logoo.png", layout="wide")
st.markdown(
    """
    <style>
    .stApp { background-color: #0b1220; color: #e6edf3; }
    .stTabs [data-baseweb="tab"] { font-weight: 600; }
    .stButton>button { border-radius: 12px; font-weight: 600; }
    .metric-card { background: #121a2b; border-radius: 16px; padding: 18px; border: 1px solid #1f2a44; }
    .soft-card { background: #0f172a; border-radius: 16px; padding: 16px; border: 1px solid #1e293b; }
    .accent { color: #7dd3fc; }
    </style>
    """,
    unsafe_allow_html=True,
)

# Top header
c1, c2 = st.columns([1, 9])
with c1:
    st.image("logoo.png", width=120)
with c2:
    st.markdown("# Gene2Therapy")
    st.caption("Gene → Enrichment → Disease → Drug repurposing")

# ----------------------------
# HTTP session with retries/backoff
# ----------------------------
RETRY = Retry(
    total=5, connect=3, read=3,
    backoff_factor=0.5,
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=["GET", "POST"],
)
SESSION = requests.Session()
SESSION.mount("https://", HTTPAdapter(max_retries=RETRY))
SESSION.mount("http://",  HTTPAdapter(max_retries=RETRY))

# ----------------------------
# Helper: load genes from any supported input
# ----------------------------
def load_genes_from_any(uploaded_file) -> list[str]:
    """
    Read genes from CSV/TSV/XLSX/TXT.
    Prefer columns 'Gene.symbol' or 'Symbol' (case-insensitive).
    Returns ≤200 unique, uppercased symbols.
    """
    name = (uploaded_file.name or "").lower()

    def _clean(series: pd.Series) -> list[str]:
        vals = series.dropna().astype(str).str.strip().str.upper()
        seen, out = set(), []
        for v in vals:
            if v and v not in seen:
                seen.add(v)
                out.append(v)
            if len(out) >= 200:
                break
        return out

    try:
        if name.endswith((".csv", ".csv.gz")):
            df = pd.read_csv(uploaded_file, compression="infer")
        elif name.endswith((".tsv", ".tsv.gz")):
            df = pd.read_csv(uploaded_file, sep="\t", compression="infer")
        elif name.endswith((".xlsx", ".xls")):
            df = pd.read_excel(uploaded_file)
        else:
            df = None

        if isinstance(df, pd.DataFrame) and not df.empty:
            lower_map = {str(c).lower(): c for c in df.columns}
            target_col = None
            for key in ("gene.symbol", "symbol"):
                if key in lower_map:
                    target_col = lower_map[key]
                    break
            if target_col is None:
                target_col = df.columns[0]
            return _clean(df[target_col])
    except Exception:
        pass

    try:
        raw = uploaded_file.read()
        text = raw.decode("utf-8", errors="ignore") if isinstance(raw, (bytes, bytearray)) else str(raw)
        lines = [ln for ln in (t.strip() for t in text.splitlines()) if ln]
        return _clean(pd.Series(lines))
    except Exception:
        return []

# ----------------------------
# Sidebar
# ----------------------------
with st.sidebar:
    st.markdown("### 🧪 BioContext")
    st.caption("Gene metadata, pathway enrichment, disease links & drug repurposing")
    st.markdown("---")
    st.markdown("**Tips**")
    st.markdown("- Keep gene lists modest (≤300) to avoid API throttling.\n- Re-run if APIs rate-limit (we cache results for 1h).")

# ----------------------------
# Caching helpers – KEGG / NCBI
# ----------------------------
@st.cache_data(ttl=3600)
def kegg_get(path: str) -> str:
    r = SESSION.get(f"https://rest.kegg.jp{path}", timeout=30)
    r.raise_for_status()
    return r.text

@st.cache_data(ttl=3600)
def ncbi_esearch_gene_ids(gene_symbol: str, organism_entrez: str) -> list[str]:
    handle = Entrez.esearch(
        db="gene",
        term=f"{gene_symbol}[Gene] AND {organism_entrez}[Organism]",
        retmode="xml",
    )
    record = Entrez.read(handle)
    handle.close()
    return record.get("IdList", [])

@st.cache_data(ttl=3600)
def ncbi_esummary_description(gene_id: str) -> str:
    handle = Entrez.esummary(db="gene", id=gene_id, retmode="xml")
    raw_xml = handle.read()
    handle.close()
    root = ET.fromstring(raw_xml)
    docsum = root.find(".//DocumentSummary")
    return (docsum.findtext("Description", default="") or "").strip()

@st.cache_data(ttl=3600)
def kegg_ncbi_to_kegg_gene_id(ncbi_gene_id: str, kegg_org_prefix: str) -> str | None:
    txt = kegg_get(f"/conv/genes/ncbi-geneid:{ncbi_gene_id}")
    if not txt.strip():
        return None
    for line in txt.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) == 2 and parts[0].endswith(f"{ncbi_gene_id}") and parts[1].startswith(f"{kegg_org_prefix}:"):
            return parts[1].strip()
    return None

@st.cache_data(ttl=3600)
def kegg_gene_pathways(kegg_gene_id: str) -> list[str]:
    txt = kegg_get(f"/link/pathway/{kegg_gene_id}")
    if not txt.strip():
        return []
    pids = []
    for line in txt.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) == 2 and parts[1].startswith("path:"):
            pids.append(parts[1])
    return pids

@st.cache_data(ttl=3600)
def kegg_pathway_name(pathway_id: str) -> str | None:
    pid = pathway_id.replace("path:", "")
    txt = kegg_get(f"/get/{pid}")
    for line in txt.split("\n"):
        if line.startswith("NAME"):
            return line.replace("NAME", "").strip()
    return None

# ----------------------------
# OpenTargets (no API key)
# ----------------------------
OT_GQL = "https://api.platform.opentargets.org/api/v4/graphql"

@st.cache_data(ttl=3600)
def ot_query(query: str, variables: dict | None = None) -> dict:
    try:
        r = SESSION.post(OT_GQL, json={"query": query, "variables": variables or {}}, timeout=40)
        data = r.json()
        if r.status_code >= 400 or (isinstance(data, dict) and data.get("errors") and not data.get("data")):
            return {}
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}

@st.cache_data(ttl=3600)
def ot_target_from_symbol(symbol: str, species: str = "Homo sapiens") -> dict | None:
    q = """
    query FindTarget($q: String!) {
      search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 5}) {
        hits { id name entity }
      }
    }
    """
    data = ot_query(q, {"q": symbol})
    hits = (((data or {}).get("data", {})).get("search", {}) or {}).get("hits", [])
    for h in hits:
        if h.get("entity") == "target":
            return {"id": h.get("id"), "approvedSymbol": h.get("name")}
    return None

@st.cache_data(ttl=3600)
def ot_diseases_for_target(ensembl_id: str, size: int = 25) -> pd.DataFrame:
    q = """
    query Associations($id: String!, $size: Int!) {
      target(ensemblId: $id) {
        id
        associatedDiseases(page: {size: $size, index: 0}) {
          rows { disease { id name } score }
        }
      }
    }
    """
    data = ot_query(q, {"id": ensembl_id, "size": size})
    rows = (((data or {}).get("data", {})).get("target", {}) or {}).get("associatedDiseases", {}).get("rows", [])
    out = []
    for r in rows:
        d = r.get("disease", {})
        out.append({
            "target": ensembl_id,
            "drug_target_ensembl": ensembl_id,
            "disease_id": d.get("id"),
            "disease_name": d.get("name"),
            "association_score": r.get("score"),
        })
    return pd.DataFrame(out)

@st.cache_data(ttl=3600)
def ot_drugs_for_target(ensembl_id: str, size: int = 50) -> pd.DataFrame:
    q = """
    query KnownDrugs($id: String!, $size: Int!) {
      target(ensemblId: $id) {
        id
        knownDrugs(size: $size) {
          rows {
            phase
            mechanismOfAction
            drug { id name }
            disease { id name }
          }
        }
      }
    }
    """
    data = ot_query(q, {"id": ensembl_id, "size": size})
    rows = (((data or {}).get("data", {})).get("target", {}) or {}).get("knownDrugs", {}).get("rows", [])
    out = []
    for r in rows:
        drug_obj = r.get("drug") or {}
        disease_obj = r.get("disease") or {}
        out.append({
            "target": ensembl_id,
            "drug_id": drug_obj.get("id"),
            "drug_name": drug_obj.get("name"),
            "phase": r.get("phase"),
            "moa": r.get("mechanismOfAction"),  # PD proxy
            "diseases": "; ".join(filter(None, [disease_obj.get("name")])),
        })
    return pd.DataFrame(out)

# ----------------------------
# PK / Approval enrichment
# ----------------------------
@st.cache_data(ttl=3600)
def chembl_find_molecule(drug_name: str) -> dict | None:
    """Find ChEMBL molecule by name/synonym; return key fields."""
    try:
        r = SESSION.get(
            "https://www.ebi.ac.uk/chembl/api/data/molecule.json",
            params={"pref_name__iexact": drug_name, "limit": 1},
            timeout=30,
        )
        js = r.json()
        if js.get("page_meta", {}).get("total_count", 0) == 0:
            r = SESSION.get(
                "https://www.ebi.ac.uk/chembl/api/data/molecule/search.json",
                params={"q": drug_name, "limit": 1},
                timeout=30,
            )
            js = r.json()
            hits = js.get("molecules", [])
            if not hits:
                return None
            mid = hits[0]["molecule_chembl_id"]
        else:
            mid = js["molecules"][0]["molecule_chembl_id"]

        r = SESSION.get(f"https://www.ebi.ac.uk/chembl/api/data/molecule/{mid}.json", timeout=30)
        m = r.json()
        return {
            "chembl_id": mid,
            "chembl_max_phase": m.get("max_phase"),
            "chembl_first_approval": m.get("first_approval"),
            "chembl_atc": "; ".join(m.get("atc_classifications") or []),
            "chembl_roa": "; ".join(m.get("dosage_form") or []) if isinstance(m.get("dosage_form"), list) else m.get("dosage_form"),
        }
    except Exception:
        return None

@st.cache_data(ttl=3600)
def pkdb_half_life(drug_name: str) -> float | None:
    """
    Try to fetch an elimination half-life (hours) from PK-DB.
    Returns a single representative value (median) if any, else None.
    """
    try:
        r = SESSION.get("https://pk-db.com/api/compounds/", params={"search": drug_name}, timeout=30)
        data = r.json()
        results = data.get("results") or []
        if not results:
            return None
        comp_id = results[0].get("id")
        if not comp_id:
            return None

        r = SESSION.get(
            "https://pk-db.com/api/parameters/",
            params={"compound": comp_id, "name": "half-life"},
            timeout=30,
        )
        params = r.json().get("results") or []
        vals = []
        for p in params:
            val = p.get("value")
            unit = (p.get("unit") or "").lower()
            if val is None:
                continue
            if "hour" in unit:
                vals.append(float(val))
            elif "min" in unit:
                vals.append(float(val) / 60.0)
        if not vals:
            return None
        return float(pd.Series(vals).median())
    except Exception:
        return None

@st.cache_data(ttl=3600)
def enrich_drugs_with_pk_and_approval(df_drugs: pd.DataFrame) -> pd.DataFrame:
    """
    Merge per-drug ChEMBL approval info and PK-DB half-life into the per-target drug rows.
    """
    if df_drugs.empty:
        return df_drugs

    unique_drugs = sorted(set(df_drugs["drug_name"].dropna().astype(str)))
    records = []
    for name in unique_drugs:
        chem = chembl_find_molecule(name) or {}
        t12 = pkdb_half_life(name)
        records.append({
            "drug_name": name,
            "chembl_id": chem.get("chembl_id"),
            "chembl_max_phase": chem.get("chembl_max_phase"),
            "chembl_first_approval": chem.get("chembl_first_approval"),
            "chembl_atc": chem.get("chembl_atc"),
            "chembl_roa": chem.get("chembl_roa"),
            "pk_half_life_h": t12,
        })
        time.sleep(0.05)

    info = pd.DataFrame.from_records(records)
    out = df_drugs.merge(info, on="drug_name", how="left")

    def _num(x):
        return pd.to_numeric(x, errors="coerce")

    def _approved(row):
        p_ot = _num(row.get("phase"))
        p_ch = _num(row.get("chembl_max_phase"))
        has_phase4 = (pd.notna(p_ot) and p_ot >= 4) or (pd.notna(p_ch) and p_ch >= 4)
        has_approval = bool(row.get("chembl_first_approval"))
        return bool(has_phase4 or has_approval)

    def _clin(row):
        p_ot = _num(row.get("phase"))
        p_ch = _num(row.get("chembl_max_phase"))
        return bool((pd.notna(p_ot) and p_ot >= 1) or (pd.notna(p_ch) and p_ch >= 1))

    out["approved"] = out.apply(_approved, axis=1)
    out["clinically_tested"] = out.apply(_clin, axis=1)

    return out

# ----------------------------
# Core functions
# ----------------------------
def fetch_gene_metadata_and_kegg(gene_list: list[str], organism_entrez: str, kegg_org_prefix: str, progress=None):
    results = []
    pathway_to_genes = defaultdict(set)
    for i, gene in enumerate(gene_list, start=1):
        if progress:
            progress.progress(min(i / max(len(gene_list), 1), 1.0))
        try:
            ids = ncbi_esearch_gene_ids(gene, organism_entrez)
            if not ids:
                results.append({"Gene": gene, "NCBI_ID": None, "Description": "No match found", "KEGG_Pathways": None})
                continue
            gene_id = ids[0]
            description = ncbi_esummary_description(gene_id)
            kegg_id = kegg_ncbi_to_kegg_gene_id(gene_id, kegg_org_prefix)
            if not kegg_id:
                results.append({"Gene": gene, "NCBI_ID": gene_id, "Description": description, "KEGG_Pathways": None})
                continue
            pids = kegg_gene_pathways(kegg_id)
            pairs = []
            for pid in pids:
                name = kegg_pathway_name(pid) or ""
                pairs.append(f"{pid} - {name}")
                pathway_to_genes[pid].add(gene)
            pathways_str = "; ".join(pairs) if pairs else None
            results.append({"Gene": gene, "NCBI_ID": gene_id, "Description": description, "KEGG_Pathways": pathways_str})
            time.sleep(0.20)
        except Exception as e:
            results.append({"Gene": gene, "NCBI_ID": None, "Description": f"Error: {e}", "KEGG_Pathways": None})
    return pd.DataFrame(results), pathway_to_genes

def compute_enrichment_counts_only(pathway_to_genes: dict):
    """
    Fast enrichment summary WITHOUT p-values/FDR.
    Just counts how many of your genes hit each KEGG pathway.
    """
    rows = []
    for pid, genes in sorted(pathway_to_genes.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        pname = kegg_pathway_name(pid) or ""
        rows.append({
            "Pathway_ID": pid.replace("path:", ""),
            "Pathway_Name": pname,
            "Count": len(genes),
            "Genes": ";".join(sorted(genes)),
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["Count", "Pathway_Name"], ascending=[False, True]).reset_index(drop=True)
    return df

# ----------------------------
# Cohort-level OpenTargets: mapping, diseases, drugs
# ----------------------------
@st.cache_data(ttl=3600)
def build_gene_to_ot_target_map(genes: list[str], species: str = "Homo sapiens") -> dict:
    g2t = {}
    for g in genes:
        hit = ot_target_from_symbol(g, species)
        if hit:
            g2t[g] = hit
        time.sleep(0.05)
    return g2t

@st.cache_data(ttl=3600)
def collect_disease_links(gene_to_target: dict) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for g, tgt in gene_to_target.items():
        tid = tgt.get("id")
        if not tid:
            continue
        df = ot_diseases_for_target(tid)
        if isinstance(df, pd.DataFrame) and not df.empty:
            df.insert(0, "gene", g)
            frames.append(df)
        time.sleep(0.05)
    if frames:
        return pd.concat(frames, ignore_index=True)
    return pd.DataFrame(columns=["gene", "drug_target_ensembl", "disease_id", "disease_name", "association_score"])

@st.cache_data(ttl=3600)
def collect_drug_suggestions(gene_to_target: dict) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for g, tgt in gene_to_target.items():
        tid = tgt.get("id")
        if not tid:
            continue
        df = ot_drugs_for_target(tid)
        if isinstance(df, pd.DataFrame) and not df.empty:
            df.insert(0, "gene", g)
            frames.append(df)
        time.sleep(0.05)
    if frames:
        return pd.concat(frames, ignore_index=True)
    return pd.DataFrame(columns=["gene", "target", "drug_id", "drug_name", "phase", "moa", "diseases"])

# ----------------------------
# UI – Inputs
# ----------------------------
st.markdown("Upload a gene list (CSV/TSV/XLSX/TXT) or paste genes, then download annotations and enrichment. Explore disease links and repurposable drugs. Now with PK/PD and approval flags.")

email = st.text_input("NCBI Entrez email (required)", value="", help="NCBI asks for a contact email for E-Utilities.")
if email:
    Entrez.email = email

organisms = {
    "Homo sapiens (human)": {"entrez": "Homo sapiens", "kegg": "hsa"},
    "Mus musculus (mouse)": {"entrez": "Mus musculus", "kegg": "mmu"},
    "Rattus norvegicus (rat)": {"entrez": "Rattus norvegicus", "kegg": "rno"},
}
org_label = st.selectbox("Organism", list(organisms.keys()), index=0)
organism_entrez = organisms[org_label]["entrez"]
kegg_org_prefix = organisms[org_label]["kegg"]

st.markdown("### Input Options")
uploaded = st.file_uploader(
    "Upload gene list (.csv, .tsv, .txt, .xlsx). If a table, I'll use the 'Gene.symbol' or 'Symbol' column.",
    type=["csv", "tsv", "txt", "xlsx"],
)
manual_input = st.text_area(
    "Or paste gene symbols here (comma, space, or newline separated):",
    placeholder="e.g. TP53, BRCA1, EGFR, MYC",
)

genes_from_input: list[str] = []
if manual_input.strip():
    raw = manual_input.replace(",", "\n").replace(" ", "\n")
    seen: set[str] = set()
    cleaned: list[str] = []
    for g in raw.splitlines():
        gg = g.strip().upper()
        if gg and gg not in seen:
            seen.add(gg)
            cleaned.append(gg)
            if len(cleaned) >= 200:
                break
    genes_from_input = cleaned
elif uploaded is not None:
    genes_from_input = load_genes_from_any(uploaded)

run_btn = st.button("Analyze", type="primary", disabled=(not genes_from_input or not email))

# ----------------------------
# Tabs
# ----------------------------
meta_tab, enrich_tab, disease_tab, drug_tab, viz_tab = st.tabs([
    "1) Metadata", "2) Enrichment", "3) Disease Links", "4) Drug Suggestions (PK/PD + Approval)", "5) Visualize"
])

if run_btn:
    # Step 1: Metadata
    with meta_tab:
        st.subheader("Step 1 — NCBI + KEGG annotations")
        progress = st.progress(0.0)
        with st.spinner("Querying NCBI and KEGG..."):
            df_meta, pathway_to_genes = fetch_gene_metadata_and_kegg(
                genes_from_input, organism_entrez, kegg_org_prefix, progress=progress
            )
        st.success("Metadata retrieval complete.")
        st.dataframe(df_meta, use_container_width=True, hide_index=True)
        st.download_button(
            "Download metadata CSV",
            data=df_meta.to_csv(index=False).encode("utf-8"),
            file_name="gene_metadata_with_kegg.csv",
            mime="text/csv",
        )

    # Step 2: Enrichment (counts-only; no p-values)
    with enrich_tab:
        st.subheader("Step 2 — Pathway Enrichment (counts-only)")
        with st.spinner("Summarizing pathway hits..."):
            df_enrich = compute_enrichment_counts_only(pathway_to_genes)

        if df_enrich.empty:
            st.info("No pathways found for enrichment with the current gene list.")
        else:
            show = df_enrich.copy()
            show.insert(0, "#", range(1, len(show) + 1))
            st.dataframe(show.reset_index(drop=True), use_container_width=True, hide_index=True)
            st.download_button(
                "Download enrichment CSV (counts-only)",
                data=df_enrich.to_csv(index=False).encode("utf-8"),
                file_name="pathway_enrichment_counts_only.csv",
                mime="text/csv",
            )
            try:
                topN = df_enrich.head(15).copy()
                fig = px.bar(topN, x="Count", y="Pathway_Name", orientation="h", title="Top pathways by hit count")
                fig.update_layout(height=600)
                st.plotly_chart(fig, use_container_width=True)
            except Exception:
                pass

    # Step 3: Disease Links
    with disease_tab:
        st.subheader("Step 3 — Disease Associations (OpenTargets)")
        with st.spinner("Mapping symbols to targets and fetching disease links..."):
            g2t = build_gene_to_ot_target_map(genes_from_input, species="Homo sapiens")
            df_dis = collect_disease_links(g2t)
        if df_dis.empty:
            st.info("No disease associations retrieved (try human genes or a smaller list).")
        else:
            agg = (
                df_dis.groupby(["disease_id", "disease_name"])
                .agg(n_genes=("gene", lambda s: len(set(s))),
                     max_score=("association_score", "max"))
                .reset_index()
                .sort_values(["n_genes", "max_score"], ascending=[False, False])
            )
            show = agg.copy()
            show.insert(0, "#", range(1, len(show) + 1))
            st.dataframe(show.reset_index(drop=True), use_container_width=True, hide_index=True)
            st.download_button(
                "Download disease associations (per gene)",
                data=df_dis.to_csv(index=False).encode("utf-8"),
                file_name="gene_disease_links_opentargets.csv",
                mime="text/csv",
            )
            st.download_button(
                "Download disease summary (aggregated)",
                data=agg.to_csv(index=False).encode("utf-8"),
                file_name="disease_summary_aggregated.csv",
                mime="text/csv",
            )

    # Step 4: Drug Suggestions + PK/PD + Approval
    with drug_tab:
        st.subheader("Step 4 — Repurposable Drugs (with PK/PD and Approval Status)")
        with st.spinner("Fetching known drugs targeting your genes..."):
            df_drugs_raw = collect_drug_suggestions(g2t)
            df_drugs = enrich_drugs_with_pk_and_approval(df_drugs_raw)

        if df_drugs.empty:
            st.info("No drugs found for the mapped targets.")
        else:
            df_drugs["phase_rank"] = pd.to_numeric(df_drugs["phase"], errors="coerce").fillna(0).astype(int)
            df_drugs["chembl_phase_rank"] = pd.to_numeric(df_drugs["chembl_max_phase"], errors="coerce").fillna(0).astype(int)
            df_drugs = df_drugs.sort_values(
                ["approved", "phase_rank", "chembl_phase_rank", "drug_name"],
                ascending=[False, False, False, True]
            )

            drug_sum = (
                df_drugs.groupby(["drug_id", "drug_name"]).agg(
                    targets=("target", lambda s: ";".join(sorted(set(s)))),
                    genes=("gene", lambda s: ";".join(sorted(set(s)))),
                    indications=("diseases", lambda s: "; ".join(sorted({x for x in "; ".join(s).split("; ") if x}))),
                    moa=("moa", lambda s: "; ".join(sorted({x for x in s if x}))),  # PD
                    max_phase=("phase", "max"),
                    approved=("approved", "max"),
                    clinically_tested=("clinically_tested", "max"),
                    chembl_id=("chembl_id", "first"),
                    chembl_first_approval=("chembl_first_approval", "first"),
                    chembl_max_phase=("chembl_max_phase", "max"),
                    chembl_atc=("chembl_atc", "first"),
                    chembl_roa=("chembl_roa", "first"),
                    pk_half_life_h=("pk_half_life_h", "median"),
                ).reset_index()
            )

            drug_sum["phase_rank"] = pd.to_numeric(drug_sum["max_phase"], errors="coerce").fillna(0).astype(int)
            drug_sum["chembl_phase_rank"] = pd.to_numeric(drug_sum["chembl_max_phase"], errors="coerce").fillna(0).astype(int)
            drug_sum = drug_sum.sort_values(
                ["approved", "clinically_tested", "phase_rank", "chembl_phase_rank", "drug_name"],
                ascending=[False, False, False, False, True]
            )

            colf1, colf2, colf3 = st.columns(3)
            with colf1:
                only_approved = st.checkbox("Only approved drugs", value=True)
            with colf2:
                include_clinical = st.checkbox("Include clinically tested (Phase ≥1)", value=True)
            with colf3:
                show_top_n = st.number_input("Show top N", min_value=5, max_value=200, value=50, step=5)

            filt = drug_sum.copy()
            if only_approved:
                filt = filt[filt["approved"] == True]
            elif include_clinical:
                filt = filt[(filt["approved"] == True) | (filt["clinically_tested"] == True)]

            st.markdown("**Recommended drugs (ranked)**")
            st.dataframe(filt.head(int(show_top_n)).reset_index(drop=True), use_container_width=True, hide_index=True)

            st.download_button(
                "Download drug suggestions (per target, with PK/PD + approval)",
                data=df_drugs.to_csv(index=False).encode("utf-8"),
                file_name="drug_suggestions_per_target_enriched.csv",
                mime="text/csv",
            )
            st.download_button(
                "Download drug suggestions (aggregated, with PK/PD + approval)",
                data=drug_sum.to_csv(index=False).encode("utf-8"),
                file_name="drug_suggestions_aggregated_enriched.csv",
                mime="text/csv",
            )

    # Step 5: Visualizations
    with viz_tab:
        st.subheader("Step 5 — Visualize the landscape")
        colA, colB = st.columns(2)

        # Sankey: Genes → Diseases (top 10) → Drugs (top ~15)
        with colA:
            try:
                if 'df_dis' in locals() and not df_dis.empty:
                    aggD = (
                        df_dis.groupby("disease_name").agg(n_genes=("gene", lambda s: len(set(s)))).reset_index()
                        .sort_values("n_genes", ascending=False).head(10)
                    )
                    top_dis = set(aggD["disease_name"].tolist())
                    genes_set = sorted(set(df_dis[df_dis["disease_name"].isin(top_dis)]["gene"]))
                    dis_list = sorted(top_dis)

                    drugs_set = []
                    if 'df_drugs' in locals() and not df_drugs.empty:
                        tmp = df_drugs[df_drugs['gene'].isin(genes_set)].copy()
                        tmp['phase_rank'] = pd.to_numeric(tmp['phase'], errors='coerce').fillna(0).astype(int)
                        tmp = tmp.sort_values('phase_rank', ascending=False)
                        drugs_set = sorted(set(tmp.head(100)['drug_name']))[:15]

                    nodes = [f"G: {g}" for g in genes_set] + [f"D: {d}" for d in dis_list] + [f"Rx: {d}" for d in drugs_set]
                    node_index = {n:i for i,n in enumerate(nodes)}

                    links = []
                    for d in dis_list:
                        sub = df_dis[df_dis["disease_name"] == d]
                        for g, cnt in Counter(sub["gene"]).items():
                            links.append((node_index[f"G: {g}"], node_index[f"D: {d}"], max(cnt, 1)))
                    if drugs_set and 'df_drugs' in locals() and not df_drugs.empty:
                        tmp = df_drugs[df_drugs['drug_name'].isin(drugs_set) & df_drugs['gene'].isin(genes_set)]
                        for _, row in tmp.iterrows():
                            s = node_index[f"G: {row['gene']}"]
                            t = node_index.get(f"Rx: {row['drug_name']}")
                            if t is not None:
                                val = int(pd.to_numeric(row.get('phase'), errors='coerce') or 0)
                                links.append((s, t, max(val, 1)))

                    sources = [s for s,_,_ in links]
                    targets = [t for _,t,_ in links]
                    values  = [v for *_,v in links]
                    fig_sankey = go.Figure(data=[go.Sankey(
                        node=dict(pad=12, thickness=14, label=nodes),
                        link=dict(source=sources, target=targets, value=values),
                    )])
                    fig_sankey.update_layout(title_text="Gene → Disease (→ Drug) connections", height=700)
                    st.plotly_chart(fig_sankey, use_container_width=True)
            except Exception as e:
                st.warning(f"Sankey could not be drawn: {e}")

        # Network: Pathway ↔ Genes
        with colB:
            try:
                if 'df_enrich' in locals() and not df_enrich.empty:
                    top_paths = df_enrich.head(8).copy()
                    edges = []
                    for _, r in top_paths.iterrows():
                        p = r["Pathway_Name"] or r["Pathway_ID"]
                        for g in (r["Genes"] or "").split(";"):
                            if g:
                                edges.append((p, g))
                    if edges:
                        G = nx.Graph()
                        G.add_edges_from(edges)
                        pos = nx.spring_layout(G, seed=42, k=0.7)
                        xe, ye = [], []
                        for u, v in G.edges():
                            xe += [pos[u][0], pos[v][0], None]
                            ye += [pos[u][1], pos[v][1], None]
                        fig_net = go.Figure()
                        fig_net.add_trace(go.Scatter(x=xe, y=ye, mode='lines', opacity=0.5))
                        fig_net.add_trace(go.Scatter(
                            x=[pos[n][0] for n in G.nodes()],
                            y=[pos[n][1] for n in G.nodes()],
                            mode='markers+text',
                            text=list(G.nodes()),
                            textposition='top center',
                        ))
                        fig_net.update_layout(title="Gene–Pathway network (top pathways by hit count)", height=700, showlegend=False)
                        st.plotly_chart(fig_net, use_container_width=True)
            except Exception as e:
                st.warning(f"Network could not be drawn: {e}")

st.markdown("---")
st.caption("APIs: NCBI E-utilities, KEGG REST, OpenTargets GraphQL, ChEMBL, PK-DB. Data is fetched live and cached for 1h. Validate findings with primary sources.")
