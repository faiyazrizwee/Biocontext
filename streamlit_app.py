# streamlit_app.py
# -------------------------------------------------------------
# BioContext – Gene2Therapy
# Gene list → KEGG enrichment (counts-only) → Disease links (OpenTargets)
# → Drug suggestions (DGIdb + DrugBank vocab CSV)
# → Visualizations
# -------------------------------------------------------------

import time
import re
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
    page_icon="logo.png",
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
    .hint { color:#9bbcff; }
    </style>
    """,
    unsafe_allow_html=True,
)

# Top header (branding)
left, right = st.columns([1, 9])
with left:
    st.image("logo.png", width=110)
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
                target_col = df.columns[0]
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
    st.caption("Gene metadata, pathway enrichment, disease links & drug suggestions")
    st.markdown("---")
    st.markdown("**Tips**")
    st.markdown("- Keep gene lists modest (≤100) to avoid API throttling.\n- Re-run if APIs rate-limit (we cache results for 1h).")

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
# OpenTargets (no API key) – disease links
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

# ----------------------------
# DrugBank vocabulary loader & normalizer (CSV you upload)
# ----------------------------
@st.cache_data(ttl=3600)
def load_drugbank_vocab(file) -> pd.DataFrame:
    try:
        df = pd.read_csv(file)
        df.columns = [c.strip() for c in df.columns]
        return df
    except Exception:
        return pd.DataFrame()

@st.cache_data(ttl=3600)
def build_drug_normalizer(vocab_df: pd.DataFrame):
    alias_to_canon, canon_info = {}, {}
    if vocab_df.empty:
        return alias_to_canon, canon_info
    for _, r in vocab_df.iterrows():
        drugbank_id = str(r.get("DrugBankID") or "").strip()
        name = str(r.get("Name") or "").strip()
        if not drugbank_id or not name:
            continue
        canon_key = f"{name} [{drugbank_id}]"
        canon_info[canon_key] = {
            "drugbank_id": drugbank_id,
            "name": name,
            "atc": r.get("ATC Level 5 Code"),
            "groups": r.get("Groups"),
        }
        aliases = set()
        for col in ["Name", "Synonyms", "Brand Names"]:
            val = r.get(col)
            if pd.isna(val):
                continue
            parts = re.split(r";|\||,|\t", str(val))
            for p in parts:
                p = p.strip()
                if p:
                    aliases.add(p.lower())
        for a in aliases:
            alias_to_canon.setdefault(a, canon_key)
    return alias_to_canon, canon_info

def normalize_drug_names(df: pd.DataFrame, name_col: str, alias_to_canon: dict, canon_info: dict) -> pd.DataFrame:
    if df.empty:
        return df
    s = df[name_col].astype(str).str.lower()
    df = df.copy()
    df["canon_key"] = s.map(alias_to_canon)
    df["drugbank_id"] = df["canon_key"].map(lambda k: canon_info.get(k, {}).get("drugbank_id") if k else None)
    df["canon_name"]  = df["canon_key"].map(lambda k: canon_info.get(k, {}).get("name") if k else None)
    df["atc_l5"]      = df["canon_key"].map(lambda k: canon_info.get(k, {}).get("atc") if k else None)
    df["groups"]      = df["canon_key"].map(lambda k: canon_info.get(k, {}).get("groups") if k else None)
    return df

# ----------------------------
# DGIdb – Gene↔Drug interactions (no API key)
# ----------------------------
DGIDB = "https://dgidb.org/api/v2"

@st.cache_data(ttl=3600)
def dgidb_interactions_for_genes(genes: list[str]) -> pd.DataFrame:
    """
    Query DGIdb interactions for a list of gene symbols.
    Returns columns: gene, drug_name, interaction_types, sources
    """
    if not genes:
        return pd.DataFrame(columns=["gene","drug_name","interaction_types","sources"])

    q = ",".join(sorted(set(g for g in genes if g)))
    try:
        r = requests.get(f"{DGIDB}/interactions.json", params={"genes": q}, timeout=40)
        r.raise_for_status()
        js = r.json()
    except Exception:
        return pd.DataFrame(columns=["gene","drug_name","interaction_types","sources"])

    rows = []
    for match in (js or {}).get("matchedTerms", []):
        gene = match.get("geneName")
        for itx in match.get("interactions", []) or []:
            drug_name = itx.get("drugName") or itx.get("drugNameSynonyms") or None
            if isinstance(drug_name, list):
                drug_name = drug_name[0] if drug_name else None
            interaction_types = "; ".join(sorted(set(itx.get("interactionTypes") or [])))
            sources = "; ".join(sorted({s.get("source") for s in itx.get("sources") or [] if s.get("source")}))
            if drug_name:
                rows.append({
                    "gene": gene,
                    "drug_name": str(drug_name),
                    "interaction_types": interaction_types,
                    "sources": sources
                })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = (df
              .groupby(["gene","drug_name"], as_index=False)
              .agg(interaction_types=("interaction_types", lambda s: "; ".join(sorted({x for x in "; ".join(s).split("; ") if x}))),
                   sources=("sources", lambda s: "; ".join(sorted({x for x in "; ".join(s).split("; ") if x})))))
    return df

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
# Cohort-level: mapping, diseases (OT)
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

# Optional: DrugBank vocabulary CSV for normalization/filters
st.markdown("#### Optional: DrugBank vocabulary CSV (for name/ID normalization)")
drugbank_csv = st.file_uploader(
    "Upload DrugBank vocabulary CSV",
    type=["csv"],
    key="db_vocab"
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
            show_meta.insert(0, "#", range(1, len(show_meta) + 1))
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

    # -------- Step 3: Disease Links (OpenTargets) --------
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
            showD.insert(0, "#", range(1, len(showD) + 1))
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

    # -------- Step 4: Drug Suggestions (DGIdb + DrugBank vocab) --------
    with drug_tab:
        st.markdown('<div class="section-title">Step 4 — Repurposable Drug Suggestions (DGIdb + DrugBank vocab)</div>', unsafe_allow_html=True)

        # Load DrugBank vocab if provided
        vocab_df = pd.DataFrame()
        alias_to_canon, canon_info = {}, {}
        if drugbank_csv is not None:
            vocab_df = load_drugbank_vocab(drugbank_csv)
            alias_to_canon, canon_info = build_drug_normalizer(vocab_df)
            if vocab_df.empty:
                st.warning("Could not read DrugBank vocabulary CSV. Proceeding without normalization.")
            else:
                st.success(f"Loaded DrugBank vocab with {len(vocab_df):,} rows for normalization.")

        with st.spinner("Fetching DGIdb gene–drug interactions..."):
            df_dgidb = dgidb_interactions_for_genes(genes)

        if df_dgidb.empty:
            st.info("No DGIdb interactions found for the provided genes.")
        else:
            # Normalize drug names (if vocab present)
            df_norm = normalize_drug_names(df_dgidb, "drug_name", alias_to_canon, canon_info) if alias_to_canon else df_dgidb.copy()

            def _agg_unique(series):
                return "; ".join(sorted({x for x in series if x}))

            drug_sum = (
                df_norm.groupby(["drug_name","canon_name","drugbank_id","atc_l5","groups"], dropna=False)
                       .agg(genes=("gene", lambda s: ";".join(sorted(set(s)))),
                            interaction_types=("interaction_types", _agg_unique),
                            sources=("sources", _agg_unique))
                       .reset_index()
            )

            only_approved = st.checkbox("Show only drugs marked 'approved' in DrugBank CSV", value=False,
                                        help="Uses the 'Groups' column in your DrugBank vocabulary file.")
            if only_approved and "groups" in drug_sum.columns:
                drug_sum = drug_sum[drug_sum["groups"].astype(str).str.contains("approved", case=False, na=False)]

            showRx = drug_sum.copy()
            showRx.insert(0, "#", range(1, len(showRx) + 1))
            cols = [c for c in ["#", "drugbank_id", "canon_name", "drug_name", "genes", "interaction_types", "sources", "atc_l5", "groups"] if c in showRx.columns]
            st.dataframe(showRx[cols], use_container_width=True, hide_index=True)

            st.download_button(
                "⬇️ Download DGIdb interactions (per gene)",
                data=df_norm.to_csv(index=False).encode("utf-8"),
                file_name="drug_suggestions_per_gene_dgidb.csv",
                mime="text/csv"
            )
            st.download_button(
                "⬇️ Download drug summary (normalized with DrugBank)",
                data=drug_sum.to_csv(index=False).encode("utf-8"),
                file_name="drug_suggestions_aggregated_dgidb_drugbank.csv",
                mime="text/csv"
            )

    # -------- Step 5: Visualizations --------
    with viz_tab:
        st.markdown('<div class="section-title">Step 5 — Visualize the landscape</div>', unsafe_allow_html=True)
        colA, colB = st.columns(2)

        # Sankey: Genes → Diseases (top 10) → Drugs (top 15 by gene count)
        with colA:
            try:
                links = []
                nodes = []
                node_index = {}

                if 'df_dis' in locals() and not df_dis.empty:
                    aggD = (
                        df_dis.groupby("disease_name").agg(n_genes=("gene", lambda s: len(set(s)))).reset_index()
                        .sort_values("n_genes", ascending=False).head(10)
                    )
                    top_dis = set(aggD["disease_name"].tolist())
                    genes_set = sorted(set(df_dis[df_dis["disease_name"].isin(top_dis)]["gene"]))
                    dis_list = sorted(top_dis)

                    drugs_set = []
                    if 'df_dgidb' in locals() and isinstance(df_dgidb, pd.DataFrame) and not df_dgidb.empty:
                        freq = df_dgidb.groupby("drug_name")["gene"].nunique().sort_values(ascending=False)
                        drugs_set = freq.head(15).index.tolist()

                    nodes = [f"G: {g}" for g in genes_set] + [f"D: {d}" for d in dis_list] + [f"Rx: {d}" for d in drugs_set]
                    node_index = {n: i for i, n in enumerate(nodes)}

                    for d in dis_list:
                        sub = df_dis[df_dis["disease_name"] == d]
                        for g, cnt in Counter(sub["gene"]).items():
                            links.append((node_index[f"G: {g}"], node_index[f"D: {d}"], max(cnt, 1)))

                    if drugs_set and 'df_dgidb' in locals() and not df_dgidb.empty:
                        tmp = df_dgidb[df_dgidb['drug_name'].isin(drugs_set) & df_dgidb['gene'].isin(genes_set)]
                        for _, row in tmp.iterrows():
                            s = node_index[f"G: {row['gene']}"]
                            t = node_index.get(f"Rx: {row['drug_name']}")
                            if t is not None:
                                links.append((s, t, 1))

                    if links:
                        sources = [s for s, _, _ in links]
                        targets = [t for _, t, _ in links]
                        values = [v for *_, v in links]
                        fig_sankey = go.Figure(data=[go.Sankey(
                            node=dict(pad=12, thickness=14, label=nodes),
                            link=dict(source=sources, target=targets, value=values),
                        )])
                        fig_sankey.update_layout(title_text="Gene → Disease (→ Drug, DGIdb) connections", height=700)
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
st.caption("APIs: NCBI E-utilities, KEGG REST, OpenTargets GraphQL (disease links), DGIdb (gene–drug). Optional DrugBank CSV for normalization/filters. Data is fetched live and cached for 1h. Validate findings with primary sources.")
