# streamlit_app.py
# -------------------------------------------------------------
# BioContext – Gene2Therapy (OpenTargets drugs)
# Gene list → KEGG enrichment (counts-only) → Disease links (OpenTargets)
# → Drug repurposing (Phase-4 filter from OpenTargets)
# → Visualizations
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
from pathlib import Path

# ----------------------------
# App Config / Theming (MUST be first Streamlit call)
# ----------------------------
st.set_page_config(
    page_title="Gene2Therapy – BioContext",
    page_icon="💊",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ---------- Global CSS (light = default; dark uses your 3 colors) ----------
st.markdown(
    """
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700;800&display=swap');

    /* brand tokens provided */
    :root{
      --silver:#C0C5C1;
      --battleship-gray:#999999;
      --night:#0b0b0b;
    }

    /* ===== LIGHT Theme (unchanged) ===== */
    :root{
      --bg:#f7f8fb; --panel:#ffffff; --glass:#ffffffea;
      --text:#000000; --muted:#4b5563; --sub:#475569;
      --border:#e6eaf2; --border-strong:#cbd5e1;
      --input-bg:#ffffff; --placeholder:#ffffff;   /* kept as you had it */
      --accent:#2563eb; --accent2:#06b6d4;
      --hero1:#1f2937; --hero2:#2563eb;
    }

    /* ===== DARK Theme ===== */
    @media (prefers-color-scheme: dark) {
      :root{
        --bg:var(--night);
        --panel:var(--night);
        --glass:var(--night);
        --text:#ffffff;
        --muted:var(--battleship-gray);
        --sub:var(--battleship-gray);
        /* borders derived from battleship-gray with alpha */
        --border:#99999933;
        --border-strong:#99999966;
        --input-bg:var(--night);
        --placeholder:var(--battleship-gray);
        /* use silver for accents/icons; button gradient set explicitly below */
        --accent:var(--silver);
        --accent2:var(--silver);
        /* title gradient hues for hero text */
        --hero1:#ffffff;
        --hero2:var(--battleship-gray);
      }
    }

    .stApp {
      background: var(--bg);
      color: var(--text);
      font-family: 'Inter', system-ui, -apple-system, 'Segoe UI', Roboto, Ubuntu, 'Helvetica Neue', Arial;
    }

    /* prevent hero clipping under Streamlit header */
    .block-container { padding-top: 3.6rem !important; padding-bottom: 2rem; }

    /* Sidebar */
    section[data-testid="stSidebar"]{
      background: var(--panel);
      border-right:1px solid var(--border);
    }
    .sidebar-title{font-weight:700; font-size:1.05rem; margin-bottom:.25rem;}
    .sidebar-tip{color:var(--muted);}

    /* Hero */
    .hero{
      margin-top:.25rem; margin-bottom:.9rem; padding:18px 20px;
      border:1px solid var(--border); border-radius:18px;
      background: linear-gradient(135deg, rgba(37,99,235,.10) 0%, rgba(34,211,238,.06) 100%);
    }
    @media (prefers-color-scheme: dark){
      .hero{
        /* soft wash using your colors */
        background: linear-gradient(135deg, rgba(59,13,17,.30) 0%, rgba(153,153,153,.12) 100%);
      }
    }
    .hero h1{
      margin:0; font-weight:800; letter-spacing:.2px; font-size:1.7rem;
      background: linear-gradient(90deg, var(--hero1), var(--hero2));
      -webkit-background-clip:text; background-clip:text; color:transparent;
    }
    .hero p{ margin:.25rem 0 0 0; color:var(--sub); }

    /* Section titles */
    .section-title{ font-weight:700; margin:0 0 .5rem 0; display:flex; align-items:center; gap:.5rem; font-size:1.05rem; }
    .section-title .icon{ color:var(--accent); }

    /* Inputs */
    label{ font-weight:600; color:var(--text); }
    .stTextInput>div>div>input,
    .stTextArea>div>textarea,
    .stSelectbox>div>div>div{
      background:var(--input-bg)!important;
      border:1.5px solid var(--border-strong)!important;
      color:var(--text)!important; border-radius:12px;
    }
    .stTextInput input::placeholder, .stTextArea textarea::placeholder{
      color:var(--placeholder)!important; opacity:1;
    }

    /* File uploader + textarea: visible borders */
    div[data-testid="stFileUploaderDropzone"]{
      min-height:140px; display:flex; align-items:center; border-radius:14px;
      background:var(--input-bg)!important;
      border:1.5px dashed var(--border-strong)!important;
    }
    .stTextArea textarea{ min-height:80px; max-height:80px; }

    /* Buttons — keep Analyze identical to light mode across themes */
    .stButton>button{
      border-radius:12px; font-weight:700; padding:.6rem 1rem;
      background: linear-gradient(90deg, #2563eb, #06b6d4) !important;  /* fixed gradient */
      color:#ffffff; border:none;
    }
    .stButton>button:disabled{ background:linear-gradient(90deg,#94a3b8,#64748b) !important; color:#e5e7eb; }

    /* Tabs & tables */
    .stTabs [data-baseweb="tab"] { font-weight:700; color:var(--sub); }
    .stTabs [aria-selected="true"] { color:var(--text); border-bottom:2px solid var(--accent); }
    .stDataFrame{ border:1px solid var(--border-strong); border-radius:12px; overflow:hidden; }

    /* Plotly readability */
    .js-plotly-plot .plotly .xtick text,
    .js-plotly-plot .plotly .ytick text,
    .js-plotly-plot .plotly .legend text,
    .js-plotly-plot .plotly .gtitle,
    .js-plotly-plot .plotly .sankey text,
    .js-plotly-plot .plotly .sankey .node text{
      fill: var(--text) !important; font-weight:700 !important;
    }
        /* --- Fix: Drug filters header getting overlapped in light theme --- */
    .drug-filters { 
      position: relative; 
      z-index: 5;          /* sit above uploader/textarea */
      clear: both;         /* start on a new row cleanly */
      margin-top: 8px;
    }

    /* keep inputs below the header in stacking order */
    div[data-testid="stFileUploaderDropzone"],
    .stTextArea, .stTextArea > div, .stTextArea textarea {
      position: relative;
      z-index: 1;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# ---------- Top header / hero ----------
def render_logo():
    for p in ("logoo.png", "logo.png", "assets/logo.png"):
        if Path(p).exists():
            st.image(p, width=96)
            return
    st.markdown("### 💊")

left, right = st.columns([1, 9])
with left:
    render_logo()
with right:
    st.markdown(
        """
        <div class="hero">
          <h1>Gene2Therapy</h1>
          <p>Fast annotations → enrichment → diseases → drug repurposing (OpenTargets)</p>
        </div>
        """,
        unsafe_allow_html=True,
    )

# ----------------------------
# Helpers
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
        pass

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
    st.markdown('<div class="sidebar-title">🧪 BioContext</div>', unsafe_allow_html=True)
    st.markdown('<div class="sidebar-tip">Gene metadata, pathway enrichment, disease links & drug repurposing</div>', unsafe_allow_html=True)
    st.markdown("---")
    st.markdown("**Tips**")
    st.markdown("- Keep gene lists modest (≤100) to avoid API throttling.\n- Re-run if APIs rate-limit (cache = 1h).")

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
        return pd.concat(frames, ignore_index=True)
    return pd.DataFrame(columns=["gene", "target", "drug_id", "drug_name", "phase", "moa", "diseases"])

# ----------------------------
# UI – Inputs
# ----------------------------
with st.container():
    st.markdown('<div class="section-title"><span class="icon">🔧</span>Input</div>', unsafe_allow_html=True)
    st.markdown('<div class="hint">Upload a gene list (CSV/TSV/XLSX/TXT) or paste genes, then explore annotations, enrichment, diseases, and drugs.</div>', unsafe_allow_html=True)

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

    col_u, col_p = st.columns(2)
    with col_u:
        uploaded = st.file_uploader("Upload gene list - csv, tsv, txt, xlsx", type=["csv", "tsv", "txt", "xlsx"])
    with col_p:
        manual_input = st.text_area(
            "Or paste gene symbols here (comma, space, or newline separated):",
            placeholder="e.g. TP53, BRCA1, EGFR, MYC",
            height=80,
        )

    st.markdown('<h4 class="drug-filters">Drug filters (applied in the Drug Suggestions tab)</h4>', unsafe_allow_html=True)
    opt_only_phase4 = st.checkbox(
        "Show only approved drugs (Phase 4)",
        value=True,
        help="Filters to max_phase ≥ 4."
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
    # Load genes
    try:
        genes = genes_from_input
        if not genes:
            st.error("Could not parse any genes. Provide a 'Gene.symbol' / 'Symbol' column or a plain list.")
            st.stop()
        st.success(f"Loaded {len(genes)} genes.")
    except Exception as e:
        st.error(f"Could not read input: {e}")
        st.stop()

    # Step 1
    with meta_tab:
        st.markdown('<div class="section-title"><span class="icon">📇</span>Step 1 — NCBI + KEGG annotations</div>', unsafe_allow_html=True)
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

    # Step 2
    with enrich_tab:
        st.markdown('<div class="section-title"><span class="icon">📊</span>Step 2 — Pathway Enrichment (counts-only)</div>', unsafe_allow_html=True)
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

    # Step 3
    with disease_tab:
        st.markdown('<div class="section-title"><span class="icon">🧬</span>Step 3 — Disease Associations (OpenTargets)</div>', unsafe_allow_html=True)
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

    # Step 4
    with drug_tab:
        st.markdown('<div class="section-title"><span class="icon">💊</span>Step 4 — Repurposable Drug Suggestions</div>', unsafe_allow_html=True)
        with st.spinner("Fetching known drugs targeting your genes..."):
            df_drugs = collect_drug_suggestions(g2t)

        if df_drugs.empty:
            st.info("No drugs found for the mapped targets.")
        else:
            df_drugs["phase_rank"] = pd.to_numeric(df_drugs["phase"], errors="coerce").fillna(0).astype(int)

            drug_sum = (
                df_drugs.groupby(["drug_id", "drug_name"]).agg(
                    targets=("target", lambda s: ";".join(sorted(set(s)))),
                    genes=("gene", lambda s: ";".join(sorted(set(s)))),
                    indications=("diseases", lambda s: "; ".join(sorted({x for x in "; ".join(s).split("; ") if x}))),
                    moa=("moa", lambda s: "; ".join(sorted({x for x in s if x}))),
                    max_phase=("phase", "max"),
                ).reset_index()
            )

            drug_sum["max_phase"] = pd.to_numeric(drug_sum["max_phase"], errors="coerce").fillna(0).astype(int)
            drug_sum["approved"] = drug_sum["max_phase"] >= 4

            if opt_only_phase4:
                drug_sum = drug_sum[drug_sum["approved"] == True]

            drug_sum = drug_sum.sort_values(
                ["approved", "max_phase", "drug_name"],
                ascending=[False, False, True]
            )

            if drug_sum.empty:
                st.info("No drugs met the selected filters. Tip: uncheck 'Show only approved (Phase 4)' to see investigational candidates.")
            else:
                showRx = drug_sum.copy()
                showRx.insert(0, "#", range(1, len(showRx) + 1))
                cols_order = [c for c in [
                    "#", "drug_id", "drug_name", "targets", "genes", "indications", "moa",
                    "max_phase", "approved"
                ] if c in showRx.columns]
                other_cols = [c for c in showRx.columns if c not in cols_order]
                st.dataframe(showRx[cols_order + other_cols], use_container_width=True, hide_index=True)

            st.download_button(
                "⬇️ Download drug suggestions (per target)",
                data=df_drugs.to_csv(index=False).encode("utf-8"),
                file_name="drug_suggestions_per_target.csv",
                mime="text/csv"
            )
            st.download_button(
                "⬇️ Download drug suggestions (aggregated + filters)",
                data=drug_sum.to_csv(index=False).encode("utf-8"),
                file_name="drug_suggestions_aggregated_filtered.csv",
                mime="text/csv"
            )

    # Step 5
    with viz_tab:
        st.markdown('<div class="section-title"><span class="icon">🌐</span>Step 5 — Visualize the landscape</div>', unsafe_allow_html=True)
        colA, colB = st.columns(2)

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
                    node_index = {n: i for i, n in enumerate(nodes)}

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
st.caption("APIs: NCBI E-utilities, KEGG REST, OpenTargets GraphQL. Data is fetched live and cached for 1h. Validate findings with primary sources.")
