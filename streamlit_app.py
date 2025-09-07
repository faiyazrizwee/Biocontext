# streamlit_app.py — Gene2Therapy (DrugCentral, dual-mode UI)
# -------------------------------------------------------------
# Gene list → KEGG enrichment (counts-only)
# → Disease links (OpenTargets)
# → Drug suggestions (DrugCentral via MyChem)
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

# =========================
# Page config
# =========================
st.set_page_config(
    page_title="Gene2Therapy – BioContext",
    page_icon="logo.png",
    layout="wide",
    initial_sidebar_state="expanded",
)

# =========================
# THEME (dark+light, equal input heights)
# =========================
st.markdown(
    """
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700;800&display=swap');

    /* ---------- Color tokens (default = DARK) ---------- */
    :root{
      --bg:#0b1220; --panel:#0e1630; --glass:#0f172aCC; --border:#1f2a44;
      --text:#e6edf3; --muted:#b8c7ef; --sub:#9bbcff;
      --accent:#2563eb; --accent2:#22d3ee;
      --success:#10b981; --warn:#f59e0b; --error:#ef4444;
      --input-bg:#0b1328; --placeholder:#9aa8c0;
    }
    /* Light scheme overrides */
    @media (prefers-color-scheme: light) {
      :root{
        --bg:#f7f8fb; --panel:#ffffff; --glass:#ffffffE6; --border:#e6eaf2;
        --text:#0b1220; --muted:#4b5563; --sub:#475569;
        --accent:#2563eb; --accent2:#06b6d4;
        --input-bg:#ffffff; --placeholder:#6b7280;
      }
    }

    .stApp { background: radial-gradient(1400px 700px at 10% -10%, rgba(37,99,235,.10) 0%, var(--bg) 45%) fixed; color: var(--text); font-family:'Inter',system-ui,-apple-system,'Segoe UI',Roboto,Ubuntu,Cantarell,'Helvetica Neue',Arial,'Noto Sans'; }
    .block-container { padding-top: 1.2rem; padding-bottom: 2rem; }

    /* Sidebar */
    section[data-testid="stSidebar"] { background: linear-gradient(180deg, var(--panel) 0%, var(--bg) 100%); border-right: 1px solid var(--border); }
    .sidebar-title{font-weight:700; font-size: 1.1rem; margin-bottom:.25rem;}
    .sidebar-tip{color:var(--muted);}

    /* Hero */
    .hero{
      margin-top:.2rem; margin-bottom:.8rem; padding:18px 20px;
      border:1px solid var(--border); border-radius:18px;
      background: linear-gradient(135deg, rgba(37,99,235,.12) 0%, rgba(34,211,238,.08) 100%);
    }
    .hero h1{
      margin:0; font-weight:800; letter-spacing:.2px; font-size:1.65rem;
      background: linear-gradient(90deg, var(--accent), var(--accent2));
      -webkit-background-clip:text; background-clip:text; color:transparent;
    }
    .hero p{ margin:.25rem 0 0 0; color:var(--sub); }

    /* Cards / sections */
    .card{ background: var(--glass); border:1px solid var(--border); border-radius:18px; padding:18px; margin-bottom:14px; box-shadow:0 10px 26px rgba(0,0,0,.18); }
    .section-title{ font-weight:700; margin:0 0 .5rem 0; display:flex; align-items:center; gap:.5rem; font-size:1.05rem; }
    .section-title .icon{ color:var(--accent); }

    /* Inputs */
    label{ font-weight:600; color:var(--text); }
    .stTextInput>div>div>input,
    .stTextArea>div>textarea,
    .stSelectbox>div>div>div{ background:var(--input-bg)!important; border:1px solid var(--border)!important; color:var(--text)!important; border-radius:12px; }
    .stTextInput input::placeholder, .stTextArea textarea::placeholder{ color:var(--placeholder) !important; opacity:1; }

    /* File uploader — make dropzone SAME HEIGHT as text area */
    .stFileUploader div[data-testid="stFileUploaderDropzone"]{
      min-height: 180px; display:flex; align-items:center; border-radius:14px;
      background:var(--input-bg)!important; border:1px dashed var(--border)!important;
    }
    .stTextArea textarea{ min-height:180px; }

    /* Buttons */
    .stButton>button{
      border-radius:12px; font-weight:700; padding:.6rem 1rem;
      background: linear-gradient(90deg, var(--accent), var(--accent2)); color:#ffffff; border:none;
      box-shadow:0 10px 20px rgba(37,99,235,.25);
    }
    .stButton>button:disabled{ background:linear-gradient(90deg,#94a3b8,#64748b); color:#e5e7eb; box-shadow:none; }

    /* Tabs */
    .stTabs [data-baseweb="tab"] { font-weight:700; color:var(--sub); }
    .stTabs [aria-selected="true"] { color:var(--text); border-bottom:2px solid var(--accent); }

    /* DataFrame */
    .stDataFrame{ border:1px solid var(--border); border-radius:12px; overflow:hidden; }

    /* Plotly: make tick labels readable in both modes */
    .js-plotly-plot .plotly .xtick text,
    .js-plotly-plot .plotly .ytick text,
    .js-plotly-plot .plotly .legend text,
    .js-plotly-plot .plotly .gtitle{ fill: var(--text) !important; }

    /* Footer */
    .footer{ color:var(--muted); font-size:.86rem; }
    </style>
    """,
    unsafe_allow_html=True,
)

# =========================
# HERO header
# =========================
left, right = st.columns([1, 9], vertical_alignment="center")
with left:
    st.image("logo.png", width=96)
with right:
    st.markdown(
        """
        <div class="hero">
          <h1>Gene2Therapy</h1>
          <p>Fast annotations → enrichment → diseases → <b>DrugCentral</b> suggestions</p>
        </div>
        """,
        unsafe_allow_html=True,
    )

# Helper for Plotly (transparent bg; axis labels get CSS override above)
def themed_bar(fig, height=560, title=None):
    fig.update_layout(
        template="plotly",  # neutral (works under both), colors overridden via CSS
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        height=height,
        title=title or fig.layout.title.text,
        margin=dict(l=40, r=20, t=60, b=40),
    )
    return fig

# =========================
# Utils
# =========================
def _extract_uniprot_id(u):
    """Return a UniProt accession from DrugCentral 'uniprot' field (dict/list/str/None)."""
    if not u:
        return None
    if isinstance(u, dict):
        return u.get("uniprot_id") or u.get("id") or u.get("accession")
    if isinstance(u, list):
        for item in u:
            val = _extract_uniprot_id(item)
            if val:
                return val
        return None
    if isinstance(u, str):
        return u.strip()
    return None

# =========================
# Load genes from any supported input
# =========================
def load_genes_from_any(uploaded_file) -> list[str]:
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
                seen.add(v); out.append(v)
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
                    target_col = lower_map[key]; break
            if target_col is None:
                target_col = df.columns[0]
            return _clean_series_to_genes(df[target_col])
    except Exception:
        pass

    try:
        try: uploaded_file.seek(0)
        except Exception: pass
        raw = uploaded_file.read()
        text = raw.decode("utf-8", errors="ignore") if isinstance(raw, (bytes, bytearray)) else str(raw)
        return _clean_series_to_genes(pd.Series([text]))
    except Exception:
        return []

# =========================
# Sidebar
# =========================
with st.sidebar:
    st.markdown('<div class="sidebar-title">🧪 BioContext</div>', unsafe_allow_html=True)
    st.markdown('<div class="sidebar-tip">Gene metadata, enrichment, disease links & drug repurposing</div>', unsafe_allow_html=True)
    st.markdown("---")
    st.markdown("**Tips**")
    st.markdown("• Keep gene lists modest (≤100)\n\n• Re-run if APIs rate-limit (cache = 1h)")

# =========================
# KEGG / NCBI
# =========================
@st.cache_data(ttl=3600)
def kegg_get(path: str) -> str:
    r = requests.get(f"https://rest.kegg.jp{path}", timeout=30); r.raise_for_status()
    return r.text

@st.cache_data(ttl=3600)
def ncbi_esearch_gene_ids(gene_symbol: str, organism_entrez: str) -> list[str]:
    handle = Entrez.esearch(db="gene", term=f"{gene_symbol}[Gene] AND {organism_entrez}[Organism]", retmode="xml")
    record = Entrez.read(handle); handle.close()
    return record.get("IdList", [])

@st.cache_data(ttl=3600)
def ncbi_esummary_description(gene_id: str) -> str:
    handle = Entrez.esummary(db="gene", id=gene_id, retmode="xml")
    raw_xml = handle.read(); handle.close()
    root = ET.fromstring(raw_xml)
    docsum = root.find(".//DocumentSummary")
    return (docsum.findtext("Description", default="") or "").strip()

@st.cache_data(ttl=3600)
def kegg_ncbi_to_kegg_gene_id(ncbi_gene_id: str, kegg_org_prefix: str) -> str | None:
    txt = kegg_get(f"/conv/genes/ncbi-geneid:{ncbi_gene_id}")
    if not txt.strip(): return None
    for line in txt.strip().split("\n"):
        parts = line.split("\t")
        if len(parts)==2 and parts[0].endswith(f"{ncbi_gene_id}") and parts[1].startswith(f"{kegg_org_prefix}:"):
            return parts[1].strip()
    return None

@st.cache_data(ttl=3600)
def kegg_gene_pathways(kegg_gene_id: str) -> list[str]:
    txt = kegg_get(f"/link/pathway/{kegg_gene_id}")
    if not txt.strip(): return []
    pids=[]
    for line in txt.strip().split("\n"):
        parts = line.split("\t")
        if len(parts)==2 and parts[1].startswith("path:"):
            pids.append(parts[1])
    return pids

@st.cache_data(ttl=3600)
def kegg_pathway_name(pathway_id: str) -> str | None:
    pid = pathway_id.replace("path:", "")
    txt = kegg_get(f"/get/{pid}")
    for line in txt.split("\n"):
        if line.startswith("NAME"):
            return line.replace("NAME","").strip()
    return None

# =========================
# OpenTargets (diseases only)
# =========================
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
    }"""
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
    }"""
    data = ot_query(q, {"id": ensembl_id, "size": size})
    rows = (((data or {}).get("data", {})).get("target", {}) or {}).get("associatedDiseases", {}).get("rows", [])
    out=[]
    for r in rows:
        d = r.get("disease", {})
        out.append({"target": ensembl_id, "disease_id": d.get("id"), "disease_name": d.get("name"), "association_score": r.get("score")})
    return pd.DataFrame(out)

# =========================
# DrugCentral via MyChem
# =========================
MYGENE_BASE = "https://mygene.info"; MYCHEM_BASE = "https://mychem.info"
MYGENE_URL = f"{MYGENE_BASE}/v3/query"; MYCHEM_URL = f"{MYCHEM_BASE}/v1/query"

@st.cache_data(ttl=300)
def _net_status() -> dict:
    out = {"mygene": None, "mychem": None, "example_hits": None}
    try: out["mygene"] = requests.get(f"{MYGENE_BASE}/status", timeout=10).status_code
    except Exception as e: out["mygene"] = f"error: {e}"
    try: out["mychem"] = requests.get(f"{MYCHEM_BASE}/status", timeout=10).status_code
    except Exception as e: out["mychem"] = f"error: {e}"
    try:
        r = requests.get(MYCHEM_URL, params={"q":"imatinib","size":1}, timeout=10)
        out["example_hits"] = r.json().get("total", None)
    except Exception as e:
        out["example_hits"] = f"error: {e}"
    return out

@st.cache_data(ttl=3600)
def _to_mygene_species_label(entrez_organism: str) -> str:
    return {"Homo sapiens":"human","Mus musculus":"mouse","Rattus norvegicus":"rat"}.get(entrez_organism,"human")

@st.cache_data(ttl=3600)
def mygene_symbol_to_uniprot(symbol: str, species_label: str = "human") -> str | None:
    params = {"q": f"symbol:{symbol} AND species:{species_label}","fields": "uniprot.Swiss-Prot,entrezgene,symbol,name","size": 1}
    try:
        r = requests.get(MYGENE_URL, params=params, timeout=20); r.raise_for_status()
        hits = r.json().get("hits", [])
        if not hits: return None
        swiss = hits[0].get("uniprot", {}).get("Swiss-Prot")
        return swiss[0] if isinstance(swiss, list) else swiss
    except Exception:
        return None

@st.cache_data(ttl=3600)
def mychem_drugcentral_query(uniprot_id: str | None = None, gene_symbol: str | None = None, size: int = 1000) -> list[dict]:
    clauses=[]
    if uniprot_id:
        clauses += [
            f"drugcentral.bioactivity.uniprot.uniprot_id:{uniprot_id}",
            f"drugcentral.targets.uniprot:{uniprot_id}",
            f"drugcentral.xref.uniprot:{uniprot_id}"
        ]
    if gene_symbol:
        clauses += [
            f"drugcentral.targets.gene_symbol:{gene_symbol}",
            f"drugcentral.bioactivity.gene_symbol:{gene_symbol}",
            f"drugcentral.bioactivity.target_name:{gene_symbol}",
            f"drugcentral.targets.target_name:{gene_symbol}",
        ]
    if not clauses: return []
    q = " OR ".join(clauses)
    fields = "drugcentral,chembl.pref_name,drugbank.name,inchikey"
    try:
        r = requests.get(MYCHEM_URL, params={"q":q,"fields":fields,"size":size}, timeout=30); r.raise_for_status()
        return r.json().get("hits", [])
    except Exception:
        return []

@st.cache_data(ttl=3600)
def collect_drug_suggestions_dc(genes: list[str], species_label: str = "human") -> pd.DataFrame:
    rows=[]
    for g in genes:
        up = mygene_symbol_to_uniprot(g, species_label=species_label)
        time.sleep(0.05)
        hits=[]
        if up: hits.extend(mychem_drugcentral_query(uniprot_id=up))
        hits.extend(mychem_drugcentral_query(gene_symbol=g))
        time.sleep(0.05)

        seen=set()
        for h in hits:
            dc = h.get("drugcentral") or {}
            if not dc: continue
            name = dc.get("drug_name") or (h.get("chembl") or {}).get("pref_name") or (h.get("drugbank") or {}).get("name") or "Unknown"
            xref = dc.get("xref") or {}
            pid = xref.get("drugbank_id") or xref.get("inchikey") or name
            key=(g,pid)
            if key in seen: continue
            seen.add(key)

            bio = dc.get("bioactivity") or []
            if isinstance(bio, dict): bio=[bio]

            actions, act_types, pots = set(), set(), []
            for b in bio:
                if not isinstance(b, dict): continue
                uni = _extract_uniprot_id(b.get("uniprot"))
                keep = (up and uni == up) or (not up)
                if not keep:
                    gs = str(b.get("gene_symbol","")).upper()
                    tn = str(b.get("target_name","")).upper()
                    keep = (gs == g.upper()) or (tn == g.upper())
                if not keep: continue

                act = b.get("action_type") or b.get("moa")
                if act: actions.add(act)
                if b.get("act_type"): act_types.add(b.get("act_type"))
                try:
                    pots.append(float(b.get("act_value")))
                except (TypeError, ValueError):
                    pass

            best_pot = min(pots) if pots else None

            atc = (dc.get("atc") or {}).get("level5") or []
            if isinstance(atc, dict): atc_list = list(atc.values())
            elif isinstance(atc, list):
                atc_list=[]
                for item in atc:
                    if isinstance(item, dict): atc_list.extend([str(v) for v in item.values()])
                    else: atc_list.append(str(item))
            else:
                atc_list=[str(atc)] if atc else []

            routes = (dc.get("pharmacology") or {}).get("routes") or []
            if isinstance(routes, str): routes=[routes]

            rows.append({
                "gene": g, "uniprot": up, "drug_name": name,
                "actions": "; ".join(sorted(actions)) or None,
                "act_types": "; ".join(sorted(act_types)) or None,
                "best_potency_nM": best_pot,
                "atc_level5": "; ".join(sorted({str(x) for x in atc_list})) or None,
                "routes": "; ".join(sorted({str(x) for x in routes})) or None,
                "approved_like": bool(atc_list),
            })

    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(columns=[
            "gene","uniprot","drug_name","actions","act_types",
            "best_potency_nM","atc_level5","routes","approved_like","potency_rank"
        ])

    def bucket(v):
        try: v=float(v)
        except (TypeError, ValueError): return 4
        return 1 if v<=10 else 2 if v<=100 else 3 if v<=1000 else 4

    df["potency_rank"] = df["best_potency_nM"].apply(bucket)
    return df

# =========================
# Core functions
# =========================
def fetch_gene_metadata_and_kegg(gene_list: list[str], organism_entrez: str, kegg_org_prefix: str, progress=None):
    results=[]; pathway_to_genes=defaultdict(set)
    for i, gene in enumerate(gene_list, start=1):
        if progress: progress.progress(min(i / max(len(gene_list), 1), 1.0))
        try:
            ids = ncbi_esearch_gene_ids(gene, organism_entrez)
            if not ids:
                results.append({"Gene": gene, "NCBI_ID": None, "Description": "No match found", "KEGG_Pathways": None}); continue
            gene_id = ids[0]
            description = ncbi_esummary_description(gene_id)
            kegg_id = kegg_ncbi_to_kegg_gene_id(gene_id, kegg_org_prefix)
            if not kegg_id:
                results.append({"Gene": gene, "NCBI_ID": gene_id, "Description": description, "KEGG_Pathways": None}); continue
            pids = kegg_gene_pathways(kegg_id)
            pairs=[]
            for pid in pids:
                name = kegg_pathway_name(pid) or ""
                pairs.append(f"{pid} - {name}")
                pathway_to_genes[pid].add(gene)
            results.append({"Gene": gene, "NCBI_ID": gene_id, "Description": description, "KEGG_Pathways": "; ".join(pairs) if pairs else None})
            time.sleep(0.20)
        except Exception as e:
            results.append({"Gene": gene, "NCBI_ID": None, "Description": f"Error: {e}", "KEGG_Pathways": None})
    return pd.DataFrame(results), pathway_to_genes

def compute_enrichment_counts_only(pathway_to_genes: dict) -> pd.DataFrame:
    rows=[]
    for pid, genes in sorted(pathway_to_genes.items(), key=lambda kv:(-len(kv[1]), kv[0])):
        rows.append({"Pathway_ID": pid.replace("path:",""),"Pathway_Name": kegg_pathway_name(pid) or "","Count": len(genes),"Genes": ";".join(sorted(genes))})
    df=pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["Count","Pathway_Name"], ascending=[False, True]).reset_index(drop=True)
    return df

# =========================
# OpenTargets cohort helpers
# =========================
@st.cache_data(ttl=3600)
def build_gene_to_ot_target_map(genes: list[str], species: str = "Homo sapiens") -> dict:
    g2t={}
    for g in genes:
        hit = ot_target_from_symbol(g, species)
        if hit: g2t[g]=hit
        time.sleep(0.05)
    return g2t

@st.cache_data(ttl=3600)
def collect_disease_links(gene_to_target: dict) -> pd.DataFrame:
    frames=[]
    for g,tgt in gene_to_target.items():
        tid=tgt.get("id")
        if not tid: continue
        df=ot_diseases_for_target(tid)
        if not df.empty:
            df.insert(0,"gene",g); frames.append(df)
        time.sleep(0.05)
    if frames: return pd.concat(frames, ignore_index=True)
    return pd.DataFrame(columns=["gene","target","disease_id","disease_name","association_score"])

# =========================
# UI — Inputs
# =========================
with st.container():
    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.markdown('<div class="section-title"><span class="icon">🔧</span>Input</div>', unsafe_allow_html=True)
    st.markdown('<div class="hint">Upload a list or paste genes; then explore annotations, enrichment, diseases, and drugs.</div>', unsafe_allow_html=True)

    email = st.text_input("NCBI Entrez email (required)", value="", help="Required by NCBI E-utilities.")
    if email: Entrez.email = email

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
        uploaded = st.file_uploader("Upload gene list - csv, tsv, txt, xlsx", type=["csv","tsv","txt","xlsx"])
    with col_p:
        manual_input = st.text_area("Or paste gene symbols here (comma, space, or newline separated)", placeholder="e.g. TP53, BRCA1, EGFR, MYC", height=120)

    st.markdown("#### Drug filters")
    opt_only_phase4 = st.checkbox("Show only approved drugs", value=True, help="Heuristic: keep drugs with ATC codes (proxy for marketed/approved).")

    genes_from_input=[]
    if manual_input.strip():
        genes_from_input = (
            pd.Series([manual_input]).str.replace(r"[,;|\\t ]+", "\\n", regex=True)
            .str.split("\\n").explode().str.strip().str.upper()
        )
        genes_from_input = [g for g in genes_from_input.tolist() if g][:200]
    elif uploaded is not None:
        genes_from_input = load_genes_from_any(uploaded)

    run_btn = st.button("▶️ Analyze", type="primary", disabled=(not genes_from_input or not email))
    st.markdown("</div>", unsafe_allow_html=True)

st.markdown("")

# =========================
# Tabs
# =========================
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
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<div class="section-title"><span class="icon">📇</span>Step 1 — NCBI + KEGG annotations</div>', unsafe_allow_html=True)
        progress = st.progress(0.0)
        with st.spinner("Querying NCBI and KEGG..."):
            df_meta, pathway_to_genes = fetch_gene_metadata_and_kegg(genes, organism_entrez, kegg_org_prefix, progress=progress)
        st.success("Metadata retrieval complete.")
        if not df_meta.empty:
            show_meta = df_meta.copy(); show_meta.insert(0, "#", range(1, len(show_meta) + 1))
            st.dataframe(show_meta, use_container_width=True, hide_index=True)
            st.download_button("⬇️ Download metadata CSV", data=df_meta.to_csv(index=False).encode("utf-8"),
                               file_name="gene_metadata_with_kegg.csv", mime="text/csv")
        else:
            st.info("No metadata found for the provided gene list.")
        st.markdown("</div>", unsafe_allow_html=True)

    # -------- Step 2: Enrichment --------
    with enrich_tab:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<div class="section-title"><span class="icon">📊</span>Step 2 — Pathway Enrichment (counts-only)</div>', unsafe_allow_html=True)
        with st.spinner("Summarizing pathway hits..."):
            df_enrich = compute_enrichment_counts_only(pathway_to_genes)
        if df_enrich.empty:
            st.info("No pathways found for enrichment with the current gene list.")
        else:
            show = df_enrich.copy(); show.insert(0, "#", range(1, len(show) + 1))
            st.dataframe(show, use_container_width=True, hide_index=True)
            st.download_button("⬇️ Download enrichment CSV (counts-only)",
                               data=df_enrich.to_csv(index=False).encode("utf-8"),
                               file_name="pathway_enrichment_counts_only.csv", mime="text/csv")
            try:
                topN = df_enrich.head(15).copy()
                fig = px.bar(topN, x="Count", y="Pathway_Name", orientation="h", title="Top pathways by gene hits")
                st.plotly_chart(themed_bar(fig, height=560), use_container_width=True)
            except Exception:
                pass
        st.markdown("</div>", unsafe_allow_html=True)

    # -------- Step 3: Disease Links --------
    with disease_tab:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<div class="section-title"><span class="icon">🧬</span>Step 3 — Disease Associations (OpenTargets)</div>', unsafe_allow_html=True)
        with st.spinner("Mapping symbols to Ensembl IDs and fetching disease links..."):
            g2t = build_gene_to_ot_target_map(genes, species="Homo sapiens")
            df_dis = collect_disease_links(g2t)

        if df_dis.empty:
            st.info("No disease associations retrieved (try human genes or a smaller list).")
        else:
            agg = (
                df_dis.groupby(["disease_id", "disease_name"])
                .agg(n_genes=("gene", lambda s: len(set(s))), max_score=("association_score", "max"))
                .reset_index().sort_values(["n_genes","max_score"], ascending=[False, False])
            )
            showD = agg.copy(); showD.insert(0, "#", range(1, len(showD) + 1))
            st.dataframe(showD, use_container_width=True, hide_index=True)

            st.download_button("⬇️ Download disease associations (per gene)",
                               data=df_dis.to_csv(index=False).encode("utf-8"),
                               file_name="gene_disease_links_opentargets.csv", mime="text/csv")
            st.download_button("⬇️ Download disease summary (aggregated)",
                               data=agg.to_csv(index=False).encode("utf-8"),
                               file_name="disease_summary_aggregated.csv", mime="text/csv")

            try:
                topD = agg.head(20)
                figd = px.bar(topD, x="n_genes", y="disease_name", orientation="h", title="Top disease associations (by # genes)")
                st.plotly_chart(themed_bar(figd, height=640), use_container_width=True)
            except Exception:
                pass
        st.markdown("</div>", unsafe_allow_html=True)

    # -------- Step 4: Drug Suggestions --------
    with drug_tab:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<div class="section-title"><span class="icon">💊</span>Step 4 — Drug Suggestions (DrugCentral)</div>', unsafe_allow_html=True)
        with st.spinner("Fetching DrugCentral suggestions via MyChem for your genes..."):
            species_label = _to_mygene_species_label(organism_entrez)
            map_preview = [{"gene": g, "uniprot": mygene_symbol_to_uniprot(g, species_label=species_label)} for g in genes[:10]]
            df_drugs = collect_drug_suggestions_dc(genes, species_label=species_label)

        with st.expander("🔎 Diagnostics (DrugCentral/MyChem)", expanded=False):
            status = _net_status()
            st.write("**Connectivity**:", status)
            st.write(pd.DataFrame(map_preview))
            st.caption("If UniProt IDs are missing we still query by gene symbol; if `mychem` != 200 or `example_hits` is 0/error, your environment cannot reach MyChem.")

        if df_drugs.empty:
            st.info("No DrugCentral hits found for these genes. Try unchecking the 'approved' filter or confirm MyChem connectivity.")
        else:
            if opt_only_phase4:
                df_drugs = df_drugs[df_drugs["approved_like"] == True]

            if df_drugs.empty:
                st.info("All hits were investigational (no ATC). Uncheck the filter to see them.")
            else:
                def agg_unique(series):
                    vals=[v for v in series if v]; parts=set()
                    for v in vals: parts.update(str(v).split("; "))
                    return "; ".join(sorted({p for p in parts if p})) if parts else None

                drug_sum = (
                    df_drugs.groupby(["drug_name"]).agg(
                        genes=("gene", lambda s: ";".join(sorted(set(s)))),
                        uniprots=("uniprot", lambda s: ";".join(sorted(set(s)))),
                        actions=("actions", agg_unique),
                        act_types=("act_types", agg_unique),
                        best_potency_nM=("best_potency_nM", "min"),
                        atc_level5=("atc_level5", agg_unique),
                        routes=("routes", agg_unique),
                        approved=("approved_like", "max"),
                    ).reset_index()
                ).sort_values(by=["approved","best_potency_nM","drug_name"], ascending=[False,True,True])

                showRx = drug_sum.copy(); showRx.insert(0, "#", range(1, len(showRx) + 1))
                cols_order = [c for c in ["#","drug_name","approved","best_potency_nM","actions","act_types","atc_level5","routes","genes","uniprots"] if c in showRx.columns]
                st.dataframe(showRx[cols_order], use_container_width=True, hide_index=True)

                st.download_button("⬇️ Download drug suggestions (per gene)",
                                   data=df_drugs.to_csv(index=False).encode("utf-8"),
                                   file_name="drug_suggestions_drugcentral_per_gene.csv", mime="text/csv")
                st.download_button("⬇️ Download drug suggestions (aggregated + filters)",
                                   data=drug_sum.to_csv(index=False).encode("utf-8"),
                                   file_name="drug_suggestions_drugcentral_aggregated.csv", mime="text/csv")
        st.markdown("</div>", unsafe_allow_html=True)

    # -------- Step 5: Visualizations --------
    with viz_tab:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<div class="section-title"><span class="icon">🌐</span>Step 5 — Visualize the landscape</div>', unsafe_allow_html=True)
        colA, colB = st.columns(2)

        with colA:
            try:
                if 'df_dis' in locals() and not df_dis.empty:
                    aggD = (df_dis.groupby("disease_name").agg(n_genes=("gene", lambda s: len(set(s)))).reset_index()
                            .sort_values("n_genes", ascending=False).head(10))
                    top_dis = set(aggD["disease_name"].tolist())
                    genes_set = sorted(set(df_dis[df_dis["disease_name"].isin(top_dis)]["gene"]))
                    dis_list = sorted(top_dis)

                    drugs_set=[]
                    if 'df_drugs' in locals() and not df_drugs.empty:
                        tmp = df_drugs[df_drugs['gene'].isin(genes_set)].copy()
                        tmp = tmp.sort_values('potency_rank', ascending=True)
                        drugs_set = sorted(set(tmp.head(100)['drug_name']))[:15]

                    nodes = [f"G: {g}" for g in genes_set] + [f"D: {d}" for d in dis_list] + [f"Rx: {d}" for d in drugs_set]
                    node_index = {n:i for i,n in enumerate(nodes)}
                    links=[]

                    for d in dis_list:
                        sub = df_dis[df_dis["disease_name"] == d]
                        for g, cnt in Counter(sub["gene"]).items():
                            links.append((node_index[f"G: {g}"], node_index[f"D: {d}"], max(cnt,1)))

                    if drugs_set and 'df_drugs' in locals() and not df_drugs.empty:
                        tmp = df_drugs[df_drugs['drug_name'].isin(drugs_set) & df_drugs['gene'].isin(genes_set)]
                        for _, row in tmp.iterrows():
                            s = node_index[f"G: {row['gene']}"]
                            t = node_index.get(f"Rx: {row['drug_name']}")
                            if t is not None:
                                pr = int(row.get('potency_rank', 4))
                                val = max(1, 5 - pr)
                                links.append((s, t, val))

                    if links:
                        sources=[s for s,_,_ in links]; targets=[t for _,t,_ in links]; values=[v for *_,v in links]
                        fig_sankey = go.Figure(data=[go.Sankey(node=dict(pad=12, thickness=14, label=nodes),
                                                               link=dict(source=sources, target=targets, value=values))])
                        st.plotly_chart(themed_bar(fig_sankey, height=680, title="Gene → Disease (→ Drug) connections"), use_container_width=True)
            except Exception as e:
                st.warning(f"Sankey could not be drawn: {e}")

        with colB:
            try:
                if 'df_enrich' in locals() and not df_enrich.empty:
                    top_paths = df_enrich.head(8).copy()
                    edges=[]
                    for _, r in top_paths.iterrows():
                        p = r["Pathway_Name"] or r["Pathway_ID"]
                        for g in (r["Genes"] or "").split(";"):
                            if g: edges.append((p,g))
                    if edges:
                        G = nx.Graph(); G.add_edges_from(edges)
                        pos = nx.spring_layout(G, seed=42, k=0.7)
                        xe, ye = [], []
                        for u,v in G.edges():
                            xe += [pos[u][0], pos[v][0], None]
                            ye += [pos[u][1], pos[v][1], None]
                        fig_net = go.Figure()
                        fig_net.add_trace(go.Scatter(x=xe, y=ye, mode='lines', opacity=0.45))
                        fig_net.add_trace(go.Scatter(
                            x=[pos[n][0] for n in G.nodes()],
                            y=[pos[n][1] for n in G.nodes()],
                            mode='markers+text', text=list(G.nodes()), textposition='top center'
                        ))
                        st.plotly_chart(themed_bar(fig_net, height=680, title="Gene–Pathway network (top pathways)"),
                                        use_container_width=True)
            except Exception as e:
                st.warning(f"Network could not be drawn: {e}")
        st.markdown("</div>", unsafe_allow_html=True)

# Footer
st.markdown("---")
st.markdown(
    '<div class="footer">APIs: NCBI E-utilities, KEGG REST, OpenTargets GraphQL (diseases), '
    'MyGene.info (mapping), MyChem.info → DrugCentral (drugs). Cache = 1h.</div>',
    unsafe_allow_html=True
)
