# streamlit_app.py
# -------------------------------------------------------------
# BioContext – Gene2Therapy
# Gene list → KEGG enrichment (counts-only) → Disease links (OpenTargets)
# → Drug repurposing (OpenTargets + ChEMBL, optional DrugCentral) → Visualizations
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

# ----------------------------
# App Config / Theming (MUST be first Streamlit call)
# ----------------------------
st.set_page_config(
    page_title="Gene2Therapy – BioContext",
    page_icon="logoo.png",
    layout="wide",
)

st.markdown(
    """
    <style>
    .stApp { background-color: #0b1220; color: #e6edf3; }
    .stTabs [data-baseweb="tab"] { font-weight: 600; }
    .stButton>button { border-radius: 12px; font-weight: 600; }
    .metric-card { background: #121a2b; border-radius: 16px; padding: 18px; border: 1px solid #1f2a44; }
    .soft-card { background: #0f172a; border-radius: 16px; padding: 16px; border: 1px solid #1e293b; }
    .accent { color: #7dd3fc; }
    .section-title { font-size: 1.1rem; margin: 0.25rem 0 0.5rem 0; }
    </style>
    """,
    unsafe_allow_html=True,
)

# Top header (branding)
left, right = st.columns([1, 9])
with left:
    st.image("logoo.png", width=110)
with right:
    st.markdown("### Gene2Therapy")
    st.caption("Gene → Enrichment → Disease → Drug repurposing")

st.divider()

# ----------------------------
# Helper: load genes from any supported input
# ----------------------------
def load_genes_from_any(uploaded_file) -> list[str]:
    """
    Read genes from CSV/TSV/XLSX/TXT.
    Prefer columns named 'Gene.symbol' or 'Symbol' (case-insensitive).
    Returns ≤200 unique, uppercased symbols.
    """
    name = (uploaded_file.name or "").lower()

    def _clean_series_to_genes(series: pd.Series) -> list[str]:
        # accept comma/space/newline/semicolon separated values
        vals = (
            series.dropna().astype(str)
            .str.replace(r"[,;|\t ]+", "\n", regex=True)
            .str.split("\n").explode()
        )
        vals = vals.astype(str).str.strip().str.upper()
        vals = vals[vals != ""]
        seen, out = set(), []
        for v in vals:
            if v and v not in seen:
                seen.add(v)
                out.append(v)
            if len(out) >= 200:
                break
        return out

    # Try dataframe-like formats first
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
                target_col = df.columns[0]  # fallback
            return _clean_series_to_genes(df[target_col])
    except Exception:
        pass  # fall back to plain text parsing

    # Plain text list
    try:
        try:
            uploaded_file.seek(0)
        except Exception:
            pass
        raw = uploaded_file.read()
        text = raw.decode("utf-8", errors="ignore") if isinstance(raw, (bytes, bytearray)) else str(raw)
        return _clean_series_to_genes(pd.Series([text]))
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
    st.markdown("- Keep gene lists modest (≤100) to avoid API throttling.\n- Re-run if APIs rate-limit (we cache results for 1h).")
    st.markdown("---")
    approved_only_toggle = st.checkbox("✅ Show only approved (ChEMBL Phase 4)", value=True,
                                       help="Filters final drug list to ChEMBL Phase 4 or with a first approval date.")
    dc_toggle = st.checkbox("🔍 DrugCentral cross-check (experimental)", value=False,
                            help="Attempts to detect your drugs in DrugCentral (best effort).")

# ----------------------------
# Caching helpers – KEGG / NCBI
# ----------------------------
@st.cache_data(ttl=3600)
def kegg_get(path: str) -> str:
    r = requests.get(f"https://rest.kegg.jp{path}", timeout=30)
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
        r = requests.post(OT_GQL, json={"query": query, "variables": variables or {}}, timeout=40)
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
            "ot_phase": r.get("phase"),  # keep OT phase separate
            "moa": r.get("mechanismOfAction"),
            "diseases": "; ".join(filter(None, [disease_obj.get("name")])),
        })
    df = pd.DataFrame(out)
    # Drop blank/None drug names early
    if not df.empty:
        df = df[df["drug_name"].notna() & df["drug_name"].astype(str).str.strip().ne("")]
        df = df.reset_index(drop=True)
    return df

# ----------------------------
# ChEMBL enrichment (approval filtering)
# ----------------------------
@st.cache_data(ttl=3600)
def chembl_find_molecule(drug_name: str) -> dict | None:
    """Find ChEMBL molecule by preferred name or full-text; return key fields."""
    try:
        # Try exact preferred name first
        r = requests.get(
            "https://www.ebi.ac.uk/chembl/api/data/molecule.json",
            params={"pref_name__iexact": drug_name, "limit": 1},
            timeout=30,
        )
        js = r.json()
        if js.get("page_meta", {}).get("total_count", 0) == 0:
            # Fallback: full-text search
            r = requests.get(
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

        # Get details
        r = requests.get(f"https://www.ebi.ac.uk/chembl/api/data/molecule/{mid}.json", timeout=30)
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
def enrich_with_chembl(df_drugs: pd.DataFrame) -> pd.DataFrame:
    """Attach ChEMBL approval info to OpenTargets drugs."""
    if df_drugs.empty:
        return df_drugs

    unique_drugs = sorted(set(df_drugs["drug_name"].dropna().astype(str)))
    records = []
    for name in unique_drugs:
        chem = chembl_find_molecule(name) or {}
        records.append({
            "drug_name": name,
            "chembl_id": chem.get("chembl_id"),
            "chembl_max_phase": chem.get("chembl_max_phase"),
            "chembl_first_approval": chem.get("chembl_first_approval"),
            "chembl_atc": chem.get("chembl_atc"),
            "chembl_roa": chem.get("chembl_roa"),
        })
        time.sleep(0.05)  # be gentle
    info = pd.DataFrame.from_records(records)
    out = df_drugs.merge(info, on="drug_name", how="left")

    def _num(x):
        return pd.to_numeric(x, errors="coerce")
    out["approved"] = (
        (_num(out["chembl_max_phase"]).fillna(0) >= 4) |
        (out["chembl_first_approval"].notna() & (out["chembl_first_approval"] != ""))
    )
    return out

# ----------------------------
# Optional DrugCentral presence check (best-effort)
# ----------------------------
@st.cache_data(ttl=3600)
def drugcentral_hit(drug_name: str) -> bool | None:
    """Return True if DrugCentral returns a hit for the drug_name; None on failure."""
    try:
        for url in [
            "https://drugcentral.org/api/activeingredient",
            "https://drugcentral.org/api/drug"
        ]:
            r = requests.get(url, params={"name": drug_name}, timeout=30)
            if r.status_code == 200:
                js = r.json()
                if (isinstance(js, list) and len(js) > 0) or (isinstance(js, dict) and js):
                    return True
        return False
    except Exception:
        return None

@st.cache_data(ttl=3600)
def enrich_with_drugcentral_presence(df_drugs: pd.DataFrame) -> pd.DataFrame:
    if df_drugs.empty:
        return df_drugs
    unique_drugs = sorted(set(df_drugs["drug_name"].dropna().astype(str)))
    rows = []
    for name in unique_drugs:
        hit = drugcentral_hit(name)
        rows.append({"drug_name": name, "drugcentral_hit": hit})
        time.sleep(0.05)
    dc = pd.DataFrame(rows)
    return df_drugs.merge(dc, on="drug_name", how="left")

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

def compute_enrichment_counts_only(pathway_to_genes: dict) -> pd.DataFrame:
    """Counts-only enrichment summary (no p-values)."""
    rows = []
    for pid, genes in sorted(pathway_to_genes.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        rows.append({
            "Pathway_ID": pid.replace("path:", ""),
            "Pathway_Name": kegg_pathway_name(pid) or "",
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
    frames = []
    for g, tgt in gene_to_target.items():
        tid = tgt.get("id")
        if not tid:
            continue
        df = ot_diseases_for_target(tid)
        if not df.empty:
            df.insert(0, "gene", g)
            frames.append(df)
        time.sleep(0.05)
    if frames:
        return pd.concat(frames, ignore_index=True)
    return pd.DataFrame(columns=["gene", "target", "disease_id", "disease_name", "association_score"])

@st.cache_data(ttl=3600)
def collect_drug_suggestions(gene_to_target: dict) -> pd.DataFrame:
    frames = []
    for g, tgt in gene_to_target.items():
        tid = tgt.get("id")
        if not tid:
            continue
        df = ot_drugs_for_target(tid)
        if not df.empty:
            df.insert(0, "gene", g)
            frames.append(df)
        time.sleep(0.05)
    if frames:
        df = pd.concat(frames, ignore_index=True)
        # Drop blank/None drug names here too (extra safety)
        df = df[df["drug_name"].notna() & df["drug_name"].astype(str).str.strip().ne("")]
        df = df.reset_index(drop=True)
        return df
    return pd.DataFrame(columns=["gene", "target", "drug_id", "drug_name", "ot_phase", "moa", "diseases"])

# ----------------------------
# UI – Inputs
# ----------------------------
st.markdown("### 🔧 Input")
st.markdown("Upload a gene list (CSV/TSV/XLSX/TXT) or paste genes, then explore annotations, enrichment, diseases, and drugs.")

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

st.markdown("#### Input Options")
uploaded = st.file_uploader(
    "Upload gene list - csv, tsv, txt, xlsx",
    type=["csv", "tsv", "txt", "xlsx"]
)
manual_input = st.text_area(
    "Or paste gene symbols here (comma, space, or newline separated):",
    placeholder="e.g. TP53, BRCA1, EGFR, MYC"
)

genes_from_input: list[str] = []
if manual_input.strip():
    genes_from_input = (
        pd.Series([manual_input])
        .str.replace(r"[,;|\t ]+", "\n", regex=True)
        .str.split("\n").explode().str.strip().str.upper()
    )
    genes_from_input = [g for g in genes_from_input.tolist() if g][:200]
elif uploaded is not None:
    genes_from_input = load_genes_from_any(uploaded)

run_btn = st.button("▶️ Analyze", type="primary", disabled=(not genes_from_input or not email))
st.divider()

# ----------------------------
# Tabs
# ----------------------------
meta_tab, enrich_tab, disease_tab, drug_tab, viz_tab = st.tabs([
    "1) Metadata", "2) Enrichment", "3) Disease Links", "4) Drug Suggestions", "5) Visualize"
])

if run_btn:
    # -------- Load genes --------
    try:
        genes = genes_from_input
        if not genes:
            st.error("Could not parse any genes. Provide a 'Gene.symbol' / 'Symbol' column or a plain list.")
            st.stop()
        st.success(f"Loaded {len(genes)} genes.")
    except Exception as e:
        st.error(f"Could not read input: {e}")
        st.stop()

    # -------- Step 1: Metadata --------
    with meta_tab:
        st.markdown('<div class="section-title">Step 1 — NCBI + KEGG annotations</div>', unsafe_allow_html=True)
        progress = st.progress(0.0)
        with st.spinner("Querying NCBI and KEGG..."):
            df_meta, pathway_to_genes = fetch_gene_metadata_and_kegg(
                genes, organism_entrez, kegg_org_prefix, progress=progress
            )
        st.success("Metadata retrieval complete.")

        if not df_meta.empty:
            show_meta = df_meta.copy()
            show_meta.insert(0, "#", range(1, len(show_meta) + 1))  # add serial numbers
            st.dataframe(show_meta, use_container_width=True, hide_index=True)
            st.download_button(
                "⬇️ Download metadata CSV",
                data=df_meta.to_csv(index=False).encode("utf-8"),
                file_name="gene_metadata_with_kegg.csv",
                mime="text/csv"
            )
        else:
            st.info("No metadata found for the provided gene list.")

    # -------- Step 2: Enrichment (counts-only) --------
    with enrich_tab:
        st.markdown('<div class="section-title">Step 2 — Pathway Enrichment (counts-only)</div>', unsafe_allow_html=True)
        with st.spinner("Summarizing pathway hits..."):
            df_enrich = compute_enrichment_counts_only(pathway_to_genes)

        if df_enrich.empty:
            st.info("No pathways found for enrichment with the current gene list.")
        else:
            show = df_enrich.copy()
            show.insert(0, "#", range(1, len(show) + 1))
            st.dataframe(show, use_container_width=True, hide_index=True)
            st.download_button(
                "⬇️ Download enrichment CSV (counts-only)",
                data=df_enrich.to_csv(index=False).encode("utf-8"),
                file_name="pathway_enrichment_counts_only.csv",
                mime="text/csv"
            )
            try:
                topN = df_enrich.head(15).copy()
                fig = px.bar(topN, x="Count", y="Pathway_Name", orientation="h", title="Top pathways by gene hits")
                fig.update_layout(height=600)
                st.plotly_chart(fig, use_container_width=True)
            except Exception:
                pass

    # -------- Step 3: Disease Links --------
    with disease_tab:
        st.markdown('<div class="section-title">Step 3 — Disease Associations (OpenTargets)</div>', unsafe_allow_html=True)
        with st.spinner("Mapping symbols to Ensembl IDs and fetching disease links..."):
            g2t = build_gene_to_ot_target_map(genes, species="Homo sapiens")
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
            showD = agg.copy()
            showD.insert(0, "#", range(1, len(showD) + 1))  # proper serial number
            st.dataframe(showD, use_container_width=True, hide_index=True)

            st.download_button(
                "⬇️ Download disease associations (per gene)",
                data=df_dis.to_csv(index=False).encode("utf-8"),
                file_name="gene_disease_links_opentargets.csv",
                mime="text/csv"
            )
            st.download_button(
                "⬇️ Download disease summary (aggregated)",
                data=agg.to_csv(index=False).encode("utf-8"),
                file_name="disease_summary_aggregated.csv",
                mime="text/csv"
            )

            try:
                topD = agg.head(20)
                figd = px.bar(topD, x="n_genes", y="disease_name", orientation="h",
                              title="Top disease associations (by # genes)")
                figd.update_layout(height=650)
                st.plotly_chart(figd, use_container_width=True)
            except Exception:
                pass

    # -------- Step 4: Drug Suggestions --------
    with drug_tab:
        st.markdown('<div class="section-title">Step 4 — Repurposable Drug Suggestions</div>', unsafe_allow_html=True)
        with st.spinner("Fetching known drugs targeting your genes..."):
            df_drugs_raw = collect_drug_suggestions(g2t)

        if df_drugs_raw.empty:
            st.info("No drugs found for the mapped targets.")
        else:
            # Enrich with ChEMBL approval info
            df_drugs = enrich_with_chembl(df_drugs_raw)

            # Optional DrugCentral presence check
            if dc_toggle:
                df_drugs = enrich_with_drugcentral_presence(df_drugs)

            # Aggregate per drug
            drug_sum = (
                df_drugs.groupby(["drug_id", "drug_name"]).agg(
                    targets=("target", lambda s: ";".join(sorted(set(s)))),
                    genes=("gene", lambda s: ";".join(sorted(set(s)))),
                    indications=("diseases", lambda s: "; ".join(sorted({x for x in "; ".join(s).split("; ") if x}))),
                    moa=("moa", lambda s: "; ".join(sorted({x for x in s if x}))),
                    max_ot_phase=("ot_phase", "max"),
                    chembl_id=("chembl_id", "first"),
                    chembl_max_phase=("chembl_max_phase", "max"),
                    chembl_first_approval=("chembl_first_approval", "first"),
                    chembl_atc=("chembl_atc", "first"),
                    chembl_roa=("chembl_roa", "first"),
                    approved=("approved", "max"),
                    drugcentral_hit=("drugcentral_hit", "max") if "drugcentral_hit" in df_drugs.columns else ("approved", "max"),
                ).reset_index()
            )

            # Clean numeric/sort
            drug_sum["chembl_max_phase"] = pd.to_numeric(drug_sum["chembl_max_phase"], errors="coerce").fillna(0).astype(int)
            drug_sum["max_ot_phase"] = pd.to_numeric(drug_sum["max_ot_phase"], errors="coerce").fillna(0).astype(int)

            # Apply sidebar choice
            view = drug_sum.copy()
            if approved_only_toggle:
                view = view[view["approved"] == True]

            # Sort: approved first, higher ChEMBL phase, then OT phase, then name
            view = view.sort_values(
                ["approved", "chembl_max_phase", "max_ot_phase", "drug_name"],
                ascending=[False, False, False, True]
            )

            # Display with serial '#'
            showRx = view.copy()
            showRx.insert(0, "#", range(1, len(showRx) + 1))

            # Preferred column order
            col_order = [
                "#", "drug_id", "drug_name", "approved",
                "chembl_max_phase", "chembl_first_approval", "chembl_id",
                "targets", "genes", "indications", "moa",
            ]
            if "drugcentral_hit" in showRx.columns:
                col_order.insert(4, "drugcentral_hit")  # show after approved

            cols_order = [c for c in col_order if c in showRx.columns]
            other_cols = [c for c in showRx.columns if c not in cols_order]
            st.markdown("**Drugs that hit at least one of your targets**")
            st.dataframe(showRx[cols_order + other_cols], use_container_width=True, hide_index=True)

            st.download_button(
                "⬇️ Download drug suggestions (per target, enriched with ChEMBL)",
                data=df_drugs.to_csv(index=False).encode("utf-8"),
                file_name="drug_suggestions_per_target_chembl.csv",
                mime="text/csv"
            )
            st.download_button(
                "⬇️ Download drug suggestions (aggregated, approved flag)",
                data=drug_sum.to_csv(index=False).encode("utf-8"),
                file_name="drug_suggestions_aggregated_chembl.csv",
                mime="text/csv"
            )

    # -------- Step 5: Visualizations --------
    with viz_tab:
        st.markdown('<div class="section-title">Step 5 — Visualize the landscape</div>', unsafe_allow_html=True)
        colA, colB = st.columns(2)

        # Sankey: Genes → Diseases (top 10) → Drugs (optional)
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
                    if 'df_drugs_raw' in locals() and not df_drugs_raw.empty:
                        tmp = df_drugs_raw[df_drugs_raw['gene'].isin(genes_set)].copy()
                        tmp['ot_phase_rank'] = pd.to_numeric(tmp['ot_phase'], errors='coerce').fillna(0).astype(int)
                        tmp = tmp.sort_values('ot_phase_rank', ascending=False)
                        drugs_set = sorted(set(tmp.head(100)['drug_name']))[:15]

                    nodes = [f"G: {g}" for g in genes_set] + [f"D: {d}" for d in dis_list] + [f"Rx: {d}" for d in drugs_set]
                    node_index = {n: i for i, n in enumerate(nodes)}

                    links = []
                    for d in dis_list:
                        sub = df_dis[df_dis["disease_name"] == d]
                        for g, cnt in Counter(sub["gene"]).items():
                            links.append((node_index[f"G: {g}"], node_index[f"D: {d}"], max(cnt, 1)))
                    if drugs_set and 'df_drugs_raw' in locals() and not df_drugs_raw.empty:
                        tmp = df_drugs_raw[df_drugs_raw['drug_name'].isin(drugs_set) & df_drugs_raw['gene'].isin(genes_set)]
                        for _, row in tmp.iterrows():
                            s = node_index[f"G: {row['gene']}"]
                            t = node_index.get(f"Rx: {row['drug_name']}")
                            if t is not None:
                                val = int(pd.to_numeric(row.get('ot_phase'), errors='coerce') or 0)
                                links.append((s, t, max(val, 1)))

                    if links:
                        sources = [s for s, _, _ in links]
                        targets = [t for _, t, _ in links]
                        values = [v for *_, v in links]
                        fig_sankey = go.Figure(data=[go.Sankey(
                            node=dict(pad=12, thickness=14, label=nodes),
                            link=dict(source=sources, target=targets, value=values),
                        )])
                        fig_sankey.update_layout(title_text="Gene → Disease (→ Drug) connections", height=700)
                        st.plotly_chart(fig_sankey, use_container_width=True)
            except Exception as e:
                st.warning(f"Sankey could not be drawn: {e}")

        # Network: Pathway ↔ Genes (top pathways)
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
st.caption("APIs: NCBI E-utilities, KEGG REST, OpenTargets GraphQL + ChEMBL (approval). Optional DrugCentral presence check. Data is fetched live and cached for 1h. Validate findings with primary sources.")
