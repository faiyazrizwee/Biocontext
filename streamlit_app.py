# streamlit_app.py
# -------------------------------------------------------------
# Gene2Therapy
# Gene list → KEGG enrichment → Disease links (OpenTargets)
# → Drug repurposing → Visualizations
# -------------------------------------------------------------

import io
import time
import math
import json
import requests
import pandas as pd
import streamlit as st
from collections import defaultdict, Counter
from xml.etree import ElementTree as ET
from Bio import Entrez
import plotly.express as px
import plotly.graph_objects as go
import networkx as nx

# NEW: robust HTTP session with retries/backoff
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# ----------------------------
# App Config / Theming (MUST be first Streamlit call)
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

# ----------------------------
# HTTP session with retries/backoff (faster & more reliable)
# ----------------------------
RETRY = Retry(
    total=5, connect=3, read=3,
    backoff_factor=0.5,
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=["GET", "POST"]
)
SESSION = requests.Session()
SESSION.mount("https://", HTTPAdapter(max_retries=RETRY))
SESSION.mount("http://",  HTTPAdapter(max_retries=RETRY))

# ----------------------------
# Helper: load genes from any supported input (define BEFORE using)
# ----------------------------
def load_genes_from_any(uploaded_file) -> list[str]:
    """
    Read genes from CSV/TSV/XLSX/TXT.
    Prefer columns named 'Gene.symbol' or 'Symbol' (case-insensitive).
    Returns up to the first 200 unique, cleaned symbols (uppercased).
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
    r = SESSION.get(f"https://rest.kegg.jp{path}", timeout=30)
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

# NEW: robust count of genes in a KEGG pathway (used for K)
@st.cache_data(ttl=3600)
def kegg_link_gene_count(pid_clean: str, retries: int = 3, backoff: float = 0.7) -> int | None:
    """
    Return number of genes on a KEGG pathway (e.g. 'hsa04310').
    Retries a few times to avoid transient empty responses / rate limits.
    """
    url = f"/link/genes/{pid_clean}"
    for i in range(retries):
        try:
            txt = kegg_get(url)
            if txt.strip():
                lines = [ln for ln in txt.splitlines() if "\t" in ln]
                return len(lines) if lines else 0
        except Exception:
            pass
        time.sleep(backoff * (i + 1))
    return None

# NEW: total KEGG genes for an organism (universe M)
@st.cache_data(ttl=3600)
def kegg_total_genes_for_org(kegg_org_prefix: str, retries: int = 3, backoff: float = 0.7) -> int | None:
    """
    Count all KEGG genes for an organism (e.g. '/list/genes/hsa').
    Returns None if KEGG is temporarily unavailable.
    """
    url = f"/list/genes/{kegg_org_prefix}"
    for i in range(retries):
        try:
            txt = kegg_get(url)
            if txt.strip():
                lines = [ln for ln in txt.splitlines() if "\t" in ln]
                return len(lines) if lines else 0
        except Exception:
            pass
        time.sleep(backoff * (i + 1))
    return None

# ----------------------------
# OpenTargets helpers (no API key)
# ----------------------------
OT_GQL = "https://api.platform.opentargets.org/api/v4/graphql"

@st.cache_data(ttl=3600)
def ot_query(query: str, variables: dict | None = None) -> dict:
    """Return {} on HTTP/GraphQL errors so the UI keeps running."""
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

def hypergeom_pval(M: int, K: int, n: int, x: int) -> float | None:
    """
    Right-tail hypergeometric P(X >= x).
    Guards and returns None only on hard errors.
    """
    try:
        M = max(int(M), 1)
        K = max(min(int(K), M), 0)
        n = max(min(int(n), M), 0)
        x = max(min(int(x), min(K, n)), 0)

        denom = math.comb(M, n)
        if denom == 0:
            return None

        if x == 0:
            return 1.0

        s = 0.0
        upper = min(K, n)
        for i in range(x, upper + 1):
            s += math.comb(K, i) * math.comb(M - K, n - i)
        return s / denom
    except Exception:
        return None

def bh_fdr(pvals: list[float]) -> list[float]:
    """Benjamini–Hochberg FDR; returns q-values in the original order."""
    m = len(pvals)
    pairs = sorted([(1.0 if p is None or p < 0 else float(p), i) for i, p in enumerate(pvals)], key=lambda x: x[0])
    q = [0.0] * m
    prev = 1.0
    for rank, (p, idx) in enumerate(pairs, start=1):
        val = min(p * m / rank, 1.0)
        prev = min(prev, val)
        q[idx] = prev
    return q

def compute_enrichment(pathway_to_genes: dict, gene_list: list[str], kegg_org_prefix: str, universe_size: int = 20000):
    """
    KEGG pathway enrichment with a hypergeometric right-tail test.

    M (universe): prefer live KEGG total for organism; fallback to user 'universe_size'.
    K: number of genes on each pathway (from /link/genes/<pathway> with retries).
    n: number of unique input genes that mapped to any KEGG pathway.
    x: number of your genes on the specific pathway.
    """
    # Universe M
    M_live = kegg_total_genes_for_org(kegg_org_prefix)
    M = M_live if (M_live is not None and M_live > 0) else int(universe_size or 20000)

    # n = unique input genes that mapped anywhere
    mapped_genes = sorted({g for s in pathway_to_genes.values() for g in s})
    n = len(mapped_genes)

    # Pre-compute K for each pathway
    K_cache: dict[str, int | None] = {}
    for pid in pathway_to_genes.keys():
        pid_clean = pid.replace("path:", "")
        K_cache[pid] = kegg_link_gene_count(pid_clean)

    rows, pvals = [], []
    for pid, genes in sorted(pathway_to_genes.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        x = len(genes)
        pname = kegg_pathway_name(pid) or ""
        K = K_cache.get(pid)

        pval = None
        if K is not None and K > 0 and n > 0 and M >= max(K, n):
            pval = hypergeom_pval(M, K, n, x)

        rows.append({
            "Pathway_ID": pid.replace("path:", ""),
            "Pathway_Name": pname,
            "Count": x,
            "Genes": ";".join(sorted(genes)),
            "PValue": pval
        })
        pvals.append(1.0 if pval is None else pval)

    df = pd.DataFrame(rows)
    if not df.empty:
        df["QValue"] = bh_fdr(pvals)
        df = df.sort_values(["PValue", "Count"], ascending=[True, False], na_position="last").reset_index(drop=True)
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
            continue# streamlit_app.py
# -------------------------------------------------------------
# Gene2Therapy
# Gene list → KEGG enrichment → Disease links (OpenTargets)
# → Drug repurposing → Visualizations
# -------------------------------------------------------------

import io
import time
import math
import json
import requests
import pandas as pd
import streamlit as st
from collections import defaultdict, Counter
from xml.etree import ElementTree as ET
from Bio import Entrez
import plotly.express as px
import plotly.graph_objects as go
import networkx as nx

# NEW: robust HTTP session with retries/backoff
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# ----------------------------
# App Config / Theming (MUST be first Streamlit call)
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

# ----------------------------
# HTTP session with retries/backoff (faster & more reliable)
# ----------------------------
RETRY = Retry(
    total=5, connect=3, read=3,
    backoff_factor=0.5,
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=["GET", "POST"]
)
SESSION = requests.Session()
SESSION.mount("https://", HTTPAdapter(max_retries=RETRY))
SESSION.mount("http://",  HTTPAdapter(max_retries=RETRY))

# ----------------------------
# Helper: load genes from any supported input (define BEFORE using)
# ----------------------------
def load_genes_from_any(uploaded_file) -> list[str]:
    """
    Read genes from CSV/TSV/XLSX/TXT.
    Prefer columns named 'Gene.symbol' or 'Symbol' (case-insensitive).
    Returns up to the first 200 unique, cleaned symbols (uppercased).
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
    r = SESSION.get(f"https://rest.kegg.jp{path}", timeout=30)
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

# NEW: robust count of genes in a KEGG pathway (used for K)
@st.cache_data(ttl=3600)
def kegg_link_gene_count(pid_clean: str, retries: int = 3, backoff: float = 0.7) -> int | None:
    """
    Return number of genes on a KEGG pathway (e.g. 'hsa04310').
    Retries a few times to avoid transient empty responses / rate limits.
    """
    url = f"/link/genes/{pid_clean}"
    for i in range(retries):
        try:
            txt = kegg_get(url)
            if txt.strip():
                lines = [ln for ln in txt.splitlines() if "\t" in ln]
                return len(lines) if lines else 0
        except Exception:
            pass
        time.sleep(backoff * (i + 1))
    return None

# NEW: total KEGG genes for an organism (universe M)
@st.cache_data(ttl=3600)
def kegg_total_genes_for_org(kegg_org_prefix: str, retries: int = 3, backoff: float = 0.7) -> int | None:
    """
    Count all KEGG genes for an organism (e.g. '/list/genes/hsa').
    Returns None if KEGG is temporarily unavailable.
    """
    url = f"/list/genes/{kegg_org_prefix}"
    for i in range(retries):
        try:
            txt = kegg_get(url)
            if txt.strip():
                lines = [ln for ln in txt.splitlines() if "\t" in ln]
                return len(lines) if lines else 0
        except Exception:
            pass
        time.sleep(backoff * (i + 1))
    return None

# ----------------------------
# OpenTargets helpers (no API key)
# ----------------------------
OT_GQL = "https://api.platform.opentargets.org/api/v4/graphql"

@st.cache_data(ttl=3600)
def ot_query(query: str, variables: dict | None = None) -> dict:
    """Return {} on HTTP/GraphQL errors so the UI keeps running."""
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

def hypergeom_pval(M: int, K: int, n: int, x: int) -> float | None:
    """
    Right-tail hypergeometric P(X >= x).
    Guards and returns None only on hard errors.
    """
    try:
        M = max(int(M), 1)
        K = max(min(int(K), M), 0)
        n = max(min(int(n), M), 0)
        x = max(min(int(x), min(K, n)), 0)

        denom = math.comb(M, n)
        if denom == 0:
            return None

        if x == 0:
            return 1.0

        s = 0.0
        upper = min(K, n)
        for i in range(x, upper + 1):
            s += math.comb(K, i) * math.comb(M - K, n - i)
        return s / denom
    except Exception:
        return None

def bh_fdr(pvals: list[float]) -> list[float]:
    """Benjamini–Hochberg FDR; returns q-values in the original order."""
    m = len(pvals)
    pairs = sorted([(1.0 if p is None or p < 0 else float(p), i) for i, p in enumerate(pvals)], key=lambda x: x[0])
    q = [0.0] * m
    prev = 1.0
    for rank, (p, idx) in enumerate(pairs, start=1):
        val = min(p * m / rank, 1.0)
        prev = min(prev, val)
        q[idx] = prev
    return q

def compute_enrichment(pathway_to_genes: dict, gene_list: list[str], kegg_org_prefix: str, universe_size: int = 20000):
    """
    KEGG pathway enrichment with a hypergeometric right-tail test.

    M (universe): prefer live KEGG total for organism; fallback to user 'universe_size'.
    K: number of genes on each pathway (from /link/genes/<pathway> with retries).
    n: number of unique input genes that mapped to any KEGG pathway.
    x: number of your genes on the specific pathway.
    """
    # Universe M
    M_live = kegg_total_genes_for_org(kegg_org_prefix)
    M = M_live if (M_live is not None and M_live > 0) else int(universe_size or 20000)

    # n = unique input genes that mapped anywhere
    mapped_genes = sorted({g for s in pathway_to_genes.values() for g in s})
    n = len(mapped_genes)

    # Pre-compute K for each pathway
    K_cache: dict[str, int | None] = {}
    for pid in pathway_to_genes.keys():
        pid_clean = pid.replace("path:", "")
        K_cache[pid] = kegg_link_gene_count(pid_clean)

    rows, pvals = [], []
    for pid, genes in sorted(pathway_to_genes.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        x = len(genes)
        pname = kegg_pathway_name(pid) or ""
        K = K_cache.get(pid)

        pval = None
        if K is not None and K > 0 and n > 0 and M >= max(K, n):
            pval = hypergeom_pval(M, K, n, x)

        rows.append({
            "Pathway_ID": pid.replace("path:", ""),
            "Pathway_Name": pname,
            "Count": x,
            "Genes": ";".join(sorted(genes)),
            "PValue": pval
        })
        pvals.append(1.0 if pval is None else pval)

    df = pd.DataFrame(rows)
    if not df.empty:
        df["QValue"] = bh_fdr(pvals)
        df = df.sort_values(["PValue", "Count"], ascending=[True, False], na_position="last").reset_index(drop=True)
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
# ----------------------------
# Cohort-level OpenTargets: mapping, diseases, drugs
# ----------------------------
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

    # fallback: empty dataframe with expected columns
    return pd.DataFrame(
        columns=["gene", "target", "disease_id", "disease_name", "association_score"]
    )


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

    # fallback: empty dataframe with expected columns
    return pd.DataFrame(
        columns=["gene", "target", "drug_id", "drug_name", "phase", "moa", "diseases"]
    )


# ----------------------------
# UI – Inputs
# ----------------------------
st.markdown("Upload a gene list (CSV/TSV/XLSX/TXT) or paste genes, then download annotations and enrichment. Explore disease links and repurposable drugs.")

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

universe_size = st.number_input(
    "Gene universe size for enrichment (approx.)",
    min_value=1000, max_value=100000, value=20000, step=1000,
    help="Used for hypergeometric p-values. ~20,000 is a common default for human protein-coding genes."
)

st.markdown("### Input Options")

# Option 1: File upload
uploaded = st.file_uploader(
    "Upload gene list (.csv, .tsv, .txt, .xlsx). If a table, I'll use the 'Gene.symbol' or 'Symbol' column.",
    type=["csv", "tsv", "txt", "xlsx"]
)

# Option 2: Manual input
manual_input = st.text_area(
    "Or paste gene symbols here (comma, space, or newline separated):",
    placeholder="e.g. TP53, BRCA1, EGFR, MYC"
)

# Combine input sources (prefer manual when provided)
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

# ----------------------------
# UI – Inputs
# ----------------------------
st.markdown("Upload a gene list (CSV/TSV/XLSX/TXT) or paste genes, then download annotations and enrichment. Explore disease links and repurposable drugs.")

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

universe_size = st.number_input(
    "Gene universe size for enrichment (approx.)",
    min_value=1000, max_value=100000, value=20000, step=1000,
    help="Used for hypergeometric p-values. ~20,000 is a common default for human protein-coding genes."
)

st.markdown("### Input Options")

# Option 1: File upload
uploaded = st.file_uploader(
    "Upload gene list (.csv, .tsv, .txt, .xlsx). If a table, I'll use the 'Gene.symbol' or 'Symbol' column.",
    type=["csv", "tsv", "txt", "xlsx"]
)

# Option 2: Manual input
manual_input = st.text_area(
    "Or paste gene symbols here (comma, space, or newline separated):",
    placeholder="e.g. TP53, BRCA1, EGFR, MYC"
)

# Combine input sources (prefer manual when provided)
genes_from_input: list[str] = []
if manual_input.strip():
    raw = manual_input.replace(",", "\n").replace(" ", "\n")
    seen = set()
    cleaned = []
    for g in raw.splitlines():
        gg = g.strip().upper()
        if gg and gg not in seen:
            seen.add(gg)
           
