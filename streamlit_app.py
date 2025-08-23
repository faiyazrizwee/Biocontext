# streamlit_app.py
# -------------------------------------------------------------
# BioContext – Gene2Therapy
# Gene list → KEGG enrichment → Disease links (OpenTargets)
# → Drug repurposing → Visualizations
# -------------------------------------------------------------

import time
import math
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
# App Config / Theming (MUST be the first Streamlit call)
# ----------------------------
st.set_page_config(
    page_title="Gene2Therapy – BioContext",
    page_icon="logoo.png",   # ensure logoo.png sits next to this file
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
    </style>
    """,
    unsafe_allow_html=True,
)

# Top header (branding)
left, right = st.columns([1, 9])
with left:
    st.image("logoo.png", width=120)
with right:
    st.markdown("# Gene2Therapy")
    st.caption("Gene → Enrichment → Disease → Drug repurposing")

# ---------- UX helpers & stats ----------
def init_state():
    for k, v in {
        "genes": [],
        "parsed_from": None,  # "paste" | "upload" | "demo"
        "df_meta": None,
        "pathway_to_genes": None,
        "df_enrich": None,
        "g2t": None,
        "df_dis": None,
        "df_drugs": None,
    }.items():
        st.session_state.setdefault(k, v)

def show_gene_chips(genes: list[str], max_show: int = 40):
    if not genes:
        return
    shown = genes[:max_show]
    chips = " ".join([
        f"<span style='background:#0f172a;border:1px solid #1f2a44;border-radius:999px;padding:4px 10px;margin:4px;display:inline-block'>{g}</span>"
        for g in shown
    ])
    extra = f" +{len(genes)-len(shown)} more" if len(genes) > len(shown) else ""
    st.markdown(chips + (f"<span style='opacity:.7;margin-left:8px'>{extra}</span>" if extra else ""),
                unsafe_allow_html=True)

def benjamini_hochberg(pvals: pd.Series) -> pd.Series:
    """Return BH FDR (q-values) for a pandas Series of p-values."""
    s = pvals.copy().astype(float)
    n = s.notna().sum()
    if n == 0:
        return s
    order = s.sort_values().index
    ranks = pd.Series(range(1, n+1), index=order, dtype=float)
    q = (s.loc[order] * n / ranks).cummin().clip(upper=1.0)
    out = pd.Series(index=s.index, dtype=float)
    out.loc[order] = q
    return out

init_state()

# ----------------------------
# Helper: load genes from any supported input (define BEFORE using)
# ----------------------------
def load_genes_from_any(uploaded_file) -> list[str]:
    """
    Read genes from CSV/TSV/XLSX/TXT.
    Prefer columns named 'Gene.symbol' or 'Symbol' (case-insensitive).
    Returns up to the first 200 unique symbols (uppercased).
    """
    name = (uploaded_file.name or "").lower()

    def _clean(series: pd.Series) -> list[str]:
        vals = (
            series.dropna()
                  .astype(str)
                  .str.strip()
                  .str.upper()
        )
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
                target_col = df.columns[0]  # fall back to first column
            return _clean(df[target_col])
    except Exception:
        pass  # fall back to plain text parsing

    # Plain-text list (one gene per line)
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
    r = requests.get(f"https://rest.kegg.jp{path}", timeout=30)
    r.raise_for_status()
    return r.text

@st.cache_data(ttl=3600)
def ncbi_esearch_gene_ids(gene_symbol: str, organism_entrez: str) -> list[str]:
    handle = Entrez.esearch(
        db="gene",
        term=f"{gene_symbol}[Gene] AND {organism_entrez}[Organism]",
        retmode="xml"
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
# OpenTargets helpers (no API key)
# ----------------------------
OT_GQL = "https://api.platform.opentargets.org/api/v4/graphql"

@st.cache_data(ttl=3600)
def ot_query(query: str, variables: dict | None = None) -> dict:
    """Return {} on HTTP/GraphQL errors so the UI keeps running."""
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
            "phase": r.get("phase"),
            "moa": r.get("mechanismOfAction"),
            "diseases": "; ".join(filter(None, [disease_obj.get("name")])),
        })
    return pd.DataFrame(out)

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

def hypergeom_pval(M: int, K: int, n: int, x: int) -> float:
    denom = math.comb(M, n) if 0 <= n <= M else 1
    s = 0.0
    upper = min(K, n)
    for i in range(x, upper + 1):
        s += math.comb(K, i) * math.comb(M - K, n - i)
    return s / denom if denom else 1.0

def compute_enrichment(pathway_to_genes: dict, gene_list: list[str], kegg_org_prefix: str, universe_size: int = 20000):
    K_cache = {}
    for pid in pathway_to_genes.keys():
        pid_clean = pid.replace("path:", "")
        try:
            txt = kegg_get(f"/link/genes/{pid_clean}")
            genes_on_pathway = {line.split("\t")[1].strip() for line in txt.strip().split("\n") if "\t" in line}
            K_cache[pid] = len(genes_on_pathway)
        except Exception:
            K_cache[pid] = None

    n = len(gene_list)
    rows = []
    for pid, genes in sorted(pathway_to_genes.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        count = len(genes)
        pname = kegg_pathway_name(pid) or ""
        K = K_cache.get(pid, None)
        pval = None
        if K is not None and universe_size is not None:
            try:
                pval = hypergeom_pval(universe_size, K, n, count)
            except Exception:
                pval = None
        rows.append({
            "Pathway_ID": pid.replace("path:", ""),
            "Pathway_Name": pname,
            "Count": count,
            "Genes": ";".join(sorted(genes)),
            "PValue": pval
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["Count", "PValue"], ascending=[False, True], na_position="last")
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
def collect_disease_links(gene_to_target: dict, size: int) -> pd.DataFrame:
    frames = []
    for g, tgt in gene_to_target.items():
        tid = tgt.get("id")
        if not tid:
            continue
        df = ot_diseases_for_target(tid, size=size)
        if not df.empty:
            df.insert(0, "gene", g)
            frames.append(df)
        time.sleep(0.05)
    if frames:
        return pd.concat(frames, ignore_index=True)
    return pd.DataFrame(columns=["gene", "target", "disease_id", "disease_name", "association_score"])

@st.cache_data(ttl=3600)
def collect_drug_suggestions(gene_to_target: dict, size: int) -> pd.DataFrame:
    frames = []
    for g, tgt in gene_to_target.items():
        tid = tgt.get("id")
        if not tid:
            continue
        df = ot_drugs_for_target(tid, size=size)
        if not df.empty:
            df.insert(0, "gene", g)
            frames.append(df)
        time.sleep(0.05)
    if frames:
        return pd.concat(frames, ignore_index=True)
    return pd.DataFrame(columns=["gene", "target", "drug_id", "drug_name", "phase", "moa", "diseases"])

# ----------------------------
# UI – Inputs (guided panel)
# ----------------------------
st.markdown("Upload a gene list (CSV/TSV/XLSX/TXT) **or** paste genes. I’ll cap at 200 unique symbols.")

colA, colB = st.columns([3, 2])
with colA:
    email = st.text_input("NCBI Entrez email (required)", value="", help="NCBI requests a contact email for E-Utilities.")
    if email:
        Entrez.email = email
with colB:
    st.write("")
    if st.button("Clear cache"):
        st.cache_data.clear()
        st.toast("Cache cleared.", icon="🧹")

org_label = st.selectbox("Organism", ["Homo sapiens (human)", "Mus musculus (mouse)", "Rattus norvegicus (rat)"], index=0)
organisms = {
    "Homo sapiens (human)": {"entrez": "Homo sapiens", "kegg": "hsa"},
    "Mus musculus (mouse)": {"entrez": "Mus musculus", "kegg": "mmu"},
    "Rattus norvegicus (rat)": {"entrez": "Rattus norvegicus", "kegg": "rno"},
}
organism_entrez = organisms[org_label]["entrez"]
kegg_org_prefix = organisms[org_label]["kegg"]

with st.expander("Advanced options", expanded=False):
    universe_size = st.number_input(
        "Gene universe size for enrichment",
        min_value=1000, max_value=100000, value=20000, step=1000,
        help="Used for hypergeometric p-values. ~20k is common for human."
    )
    ot_size_diseases = st.slider("Max diseases per target (OpenTargets)", 10, 200, 50, 10)
    ot_size_drugs    = st.slider("Max drugs per target (OpenTargets)", 10, 200, 100, 10)

st.markdown("### Input Options")
src = st.radio("Choose input source", ["Upload a file", "Paste genes", "Use demo list"], horizontal=True)

genes_from_input: list[str] = []
if src == "Upload a file":
    uploaded = st.file_uploader(
        "Upload (.csv, .tsv, .txt, .xlsx). If a table, I'll use the 'Gene.symbol' or 'Symbol' column.",
        type=["csv", "tsv", "txt", "xlsx"]
    )
    if uploaded is not None:
        genes_from_input = load_genes_from_any(uploaded)
        st.session_state.parsed_from = "upload"
elif src == "Paste genes":
    manual_input = st.text_area(
        "Paste gene symbols (comma/space/newline separated):",
        placeholder="TP53, BRCA1, EGFR, MYC, ...",
        height=120
    )
    if manual_input.strip():
        raw = manual_input.replace(",", "\n").replace(" ", "\n")
        seen, cleaned = set(), []
        for g in raw.splitlines():
            gg = g.strip().upper()
            if gg and gg not in seen:
                seen.add(gg); cleaned.append(gg)
            if len(cleaned) >= 200:
                break
        genes_from_input = cleaned
        st.session_state.parsed_from = "paste"
else:
    demo = ["TP53","BRCA1","EGFR","MYC","PIK3CA","PTEN","KRAS","BRAF","CDK4","CDK6","CDKN2A","ERBB2","NRAS","ALK","ROS1"]
    genes_from_input = demo
    st.session_state.parsed_from = "demo"
    st.info("Loaded a small demo list so you can try the app.")

# Preview chips + persist
if genes_from_input:
    st.session_state.genes = genes_from_input
    st.toast(f"Parsed {len(genes_from_input)} gene(s).", icon="✅")
    show_gene_chips(genes_from_input)

# Final go button
run_btn = st.button(
    "Analyze",
    type="primary",
    disabled=(not st.session_state.genes or not email)
)

# ----------------------------
# Tabs
# ----------------------------
meta_tab, enrich_tab, disease_tab, drug_tab, viz_tab = st.tabs([
    "1) Metadata", "2) Enrichment", "3) Disease Links", "4) Drug Suggestions", "5) Visualize"
])

if run_btn:
    # -------- Load genes --------
    try:
        genes = st.session_state.genes
        if not genes:
            st.error("No genes parsed. Upload a file or paste symbols.")
            st.stop()
        st.success(f"Loaded {len(genes)} genes.")
    except Exception as e:
        st.error(f"Could not read input: {e}")
        st.stop()

    # -------- Step 1: Metadata --------
    with meta_tab:
        st.subheader("Step 1 — NCBI + KEGG annotations")
        progress = st.progress(0.0)
        with st.spinner("Querying NCBI and KEGG..."):
            df_meta, pathway_to_genes = fetch_gene_metadata_and_kegg(
                genes, organism_entrez, kegg_org_prefix, progress=progress
            )
        st.toast("Metadata retrieval complete.", icon="📘")
        st.dataframe(df_meta, use_container_width=True)
        st.download_button(
            "Download metadata CSV",
            data=df_meta.to_csv(index=False).encode("utf-8"),
            file_name="gene_metadata_with_kegg.csv",
            mime="text/csv"
        )

    # -------- Step 2: Enrichment --------
    with enrich_tab:
        st.subheader("Step 2 — Pathway Enrichment (KEGG)")
        with st.spinner("Computing enrichment..."):
            df_enrich = compute_enrichment(pathway_to_genes, genes, kegg_org_prefix, universe_size=universe_size)

        if df_enrich.empty:
            st.info("No pathways found for enrichment with the current gene list.")
        else:
            # Add BH FDR (q-value)
            if "PValue" in df_enrich.columns:
                df_enrich["q_value"] = benjamini_hochberg(df_enrich["PValue"])

            # Filters
            colf1, colf2, colf3 = st.columns([1,1,2])
            with colf1:
                min_count = st.number_input("Min genes in pathway", 1, 100, 2)
            with colf2:
                q_cut = st.selectbox("q-value max", ["none", 0.05, 0.10, 0.20], index=1)
            with colf3:
                search = st.text_input("Search pathway name (contains)", "")

            dfv = df_enrich.copy()
            dfv = dfv[dfv["Count"] >= min_count]
            if q_cut != "none" and "q_value" in dfv:
                dfv = dfv[dfv["q_value"].le(float(q_cut))]
            if search:
                dfv = dfv[dfv["Pathway_Name"].str.contains(search, case=False, na=False)]

            st.dataframe(
                dfv,
                use_container_width=True,
                column_config={
                    "PValue": st.column_config.NumberColumn(format="%.3e"),
                    "q_value": st.column_config.NumberColumn(format="%.3e", help="Benjamini–Hochberg FDR"),
                    "Count": st.column_config.NumberColumn(help="# of your genes in this pathway")
                }
            )
            st.download_button(
                "Download enrichment CSV",
                data=dfv.to_csv(index=False).encode("utf-8"),
                file_name="pathway_enrichment_filtered.csv",
                mime="text/csv"
            )
            try:
                topN = dfv.head(15).copy()
                fig = px.bar(topN, x="Count", y="Pathway_Name", orientation="h",
                             title="Top enriched pathways (filtered)")
                fig.update_layout(height=600)
                st.plotly_chart(fig, use_container_width=True)
            except Exception:
                pass

    # -------- Step 3: Disease Links --------
    with disease_tab:
        st.subheader("Step 3 — Disease Associations (OpenTargets)")
        with st.spinner("Mapping symbols to Ensembl IDs and fetching disease links..."):
            g2t = build_gene_to_ot_target_map(genes, species="Homo sapiens")
            df_dis = collect_disease_links(g2t, size=ot_size_diseases)

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

            # Filters
            topN_dis = st.slider("Show top N diseases", 10, 200, 50, 10)
            min_genes_dis = st.number_input("Min #genes per disease", 1, 100, 2)
            search_dis = st.text_input("Search disease name", "")

            aggv = agg.copy()
            aggv = aggv[aggv["n_genes"] >= min_genes_dis]
            if search_dis:
                aggv = aggv[aggv["disease_name"].str.contains(search_dis, case=False, na=False)]

            st.markdown("**Top diseases hit by your gene list**")
            st.dataframe(aggv.head(topN_dis), use_container_width=True)

            st.download_button(
                "Download disease associations (per gene)",
                data=df_dis.to_csv(index=False).encode("utf-8"),
                file_name="gene_disease_links_opentargets.csv",
                mime="text/csv"
            )
            st.download_button(
                "Download disease summary (aggregated)",
                data=agg.to_csv(index=False).encode("utf-8"),
                file_name="disease_summary_aggregated.csv",
                mime="text/csv"
            )
            try:
                topD = aggv.head(min(topN_dis, 20))
                figd = px.bar(topD, x="n_genes", y="disease_name", orientation="h",
                              title="Top disease associations (by #genes)")
                figd.update_layout(height=650)
                st.plotly_chart(figd, use_container_width=True)
            except Exception:
                pass

    # -------- Step 4: Drug Suggestions --------
    with drug_tab:
        st.subheader("Step 4 — Repurposable Drug Suggestions (targets from your list)")
        with st.spinner("Fetching known drugs targeting your genes..."):
            df_drugs = collect_drug_suggestions(g2t, size=ot_size_drugs)

        if df_drugs.empty:
            st.info("No drugs found for the mapped targets. Try different genes or check human mapping.")
        else:
            phase_rank = {None: 0, 1: 1, 2: 2, 3: 3, 4: 4}
            df_drugs["phase_rank"] = df_drugs["phase"].map(phase_rank).fillna(0)
            drug_sum = (
                df_drugs.groupby(["drug_id", "drug_name"]).agg(
                    targets=("target", lambda s: ";".join(sorted(set(s)))),
                    genes=("gene", lambda s: ";".join(sorted(set(s)))),
                    max_phase=("phase", "max"),
                    max_phase_rank=("phase_rank", "max"),
                    indications=("diseases", lambda s: "; ".join(sorted({x for x in "; ".join(s).split("; ") if x}))),
                    moa=("moa", lambda s: "; ".join(sorted({x for x in s if x})))
                ).reset_index()
                 .sort_values(["max_phase_rank", "drug_name"], ascending=[False, True])
            )

            # Quick filters
            min_phase = st.selectbox("Minimum clinical phase", [0,1,2,3,4], index=3,
                                     help="Filter by highest phase per drug")
            search_drug = st.text_input("Search drug/target/gene", "")

            dv = drug_sum[drug_sum["max_phase"].fillna(0).astype(int) >= int(min_phase)].copy()
            if search_drug:
                mask = (
                    dv["drug_name"].str.contains(search_drug, case=False, na=False) |
                    dv["targets"].str.contains(search_drug, case=False, na=False) |
                    dv["genes"].str.contains(search_drug, case=False, na=False)
                )
                dv = dv[mask]

            st.markdown("**Drugs that hit at least one of your targets** (higher phase first)")
            st.dataframe(dv, use_container_width=True)
            st.download_button(
                "Download drug suggestions (per target)",
                data=df_drugs.to_csv(index=False).encode("utf-8"),
                file_name="drug_suggestions_per_target.csv",
                mime="text/csv"
            )
            st.download_button(
                "Download drug suggestions (aggregated)",
                data=drug_sum.to_csv(index=False).encode("utf-8"),
                file_name="drug_suggestions_aggregated.csv",
                mime="text/csv"
            )

    # -------- Step 5: Visualizations --------
    with viz_tab:
        st.subheader("Step 5 — Visualize the landscape")
        colA, colB = st.columns(2)

        # Sankey: Genes → Diseases (top 10) → (optional) Drugs
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
                        tmp['phase_rank'] = tmp['phase'].map({None:0,1:1,2:2,3:3,4:4}).fillna(0)
                        tmp = tmp.sort_values('phase_rank', ascending=False)
                        drugs_set = sorted(set(tmp.head(100)['drug_name']))[:15]
                    nodes = (
                        [f"G: {g}" for g in genes_set] +
                        [f"D: {d}" for d in dis_list] +
                        ([f"Rx: {d}" for d in drugs_set] if drugs_set else [])
                    )
                    node_index = {n:i for i,n in enumerate(nodes)}
                    links_s = []
                    for d in dis_list:
                        sub = df_dis[df_dis["disease_name"] == d]
                        for g, cnt in Counter(sub["gene"]).items():
                            s = node_index[f"G: {g}"]
                            t = node_index[f"D: {d}"]
                            links_s.append((s, t, max(cnt, 1)))
                    links_r = []
                    if drugs_set and 'df_drugs' in locals() and not df_drugs.empty:
                        tmp = df_drugs[df_drugs['drug_name'].isin(drugs_set) & df_drugs['gene'].isin(genes_set)]
                        for _, row in tmp.iterrows():
                            s = node_index[f"G: {row['gene']}"]
                            t = node_index.get(f"Rx: {row['drug_name']}")
                            if t is not None:
                                links_r.append((s, t, max((row['phase'] or 0), 1)))
                    sources = [s for s,_,_ in links_s + links_r]
                    targets = [t for _,t,_ in links_s + links_r]
                    values  = [v for *_,v in links_s + links_r]
                    fig_sankey = go.Figure(data=[go.Sankey(
                        node=dict(pad=12, thickness=14, label=nodes),
                        link=dict(source=sources, target=targets, value=values)
                    )])
                    fig_sankey.update_layout(title_text="Gene → Disease (→ Drug) connections", height=700)
                    st.plotly_chart(fig_sankey, use_container_width=True)
            except Exception as e:
                st.warning(f"Sankey could not be drawn: {e}")

        # Network: Pathway ↔ Genes (top enriched)
        with colB:
            try:
                if 'df_enrich' in locals() and not df_enrich.empty:
                    top_paths = df_enrich.head(8).copy()
                    edges = []
                    for _, r in top_paths.iterrows():
                        p = r["Pathway_Name"] or r["Pathway_ID"]
                        glist = (r["Genes"] or "").split(";")
                        for g in glist:
                            if g:
                                edges.append((p, g))
                    if edges:
                        G = nx.Graph()
                        G.add_edges_from(edges)
                        pos = nx.spring_layout(G, seed=42, k=0.7)
                        x_nodes = [pos[n][0] for n in G.nodes()]
                        y_nodes = [pos[n][1] for n in G.nodes()]
                        node_text = list(G.nodes())
                        xe, ye = [], []
                        for u, v in G.edges():
                            xe += [pos[u][0], pos[v][0], None]
                            ye += [pos[u][1], pos[v][1], None]
                        fig_net = go.Figure()
                        fig_net.add_trace(go.Scatter(x=xe, y=ye, mode='lines', opacity=0.5))
                        fig_net.add_trace(go.Scatter(x=x_nodes, y=y_nodes, mode='markers+text',
                                                     text=node_text, textposition='top center'))
                        fig_net.update_layout(title="Gene–Pathway network (top enriched)", height=700, showlegend=False)
                        st.plotly_chart(fig_net, use_container_width=True)
            except Exception as e:
                st.warning(f"Network could not be drawn: {e}")

st.markdown("---")
st.caption("APIs used: NCBI E-utilities, KEGG REST, OpenTargets GraphQL. Data is fetched live and cached for 1h. For clinical use, validate findings with primary sources.")
