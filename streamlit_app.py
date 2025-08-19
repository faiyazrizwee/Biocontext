# streamlit_app.py
import io
import time
import math
import requests
import pandas as pd
import streamlit as st
from collections import defaultdict
from xml.etree import ElementTree as ET
from Bio import Entrez

# ----------------------------
# App Config
# ----------------------------
st.set_page_config(page_title="BioContext – Gene Metadata & Pathway Enrichment", layout="wide")

# ----------------------------
# Caching helpers
# ----------------------------
@st.cache_data(ttl=3600)
def kegg_get(path: str) -> str:
    """Cached GET to KEGG REST to avoid re-downloading identical results."""
    r = requests.get(f"https://rest.kegg.jp{path}", timeout=30)
    r.raise_for_status()
    return r.text

@st.cache_data(ttl=3600)
def ncbi_esearch_gene_ids(gene_symbol: str, organism_entrez: str) -> list[str]:
    """Search NCBI Gene for symbol in organism; return a list of Gene IDs (strings)."""
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
    """Get the Description from NCBI Gene esummary."""
    handle = Entrez.esummary(db="gene", id=gene_id, retmode="xml")
    raw_xml = handle.read()
    handle.close()
    root = ET.fromstring(raw_xml)
    docsum = root.find(".//DocumentSummary")
    return (docsum.findtext("Description", default="") or "").strip()

@st.cache_data(ttl=3600)
def kegg_ncbi_to_kegg_gene_id(ncbi_gene_id: str, kegg_org_prefix: str) -> str | None:
    """
    Convert an NCBI Gene ID to a KEGG gene ID for the organism.
    Example: ncbi-geneid:7157 -> hsa:7157
    """
    txt = kegg_get(f"/conv/genes/ncbi-geneid:{ncbi_gene_id}")
    # Lines like: ncbi-geneid:7157\thsa:7157
    if not txt.strip():
        return None
    for line in txt.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) == 2 and parts[0].endswith(f"{ncbi_gene_id}") and parts[1].startswith(f"{kegg_org_prefix}:"):
            return parts[1].strip()
    return None

@st.cache_data(ttl=3600)
def kegg_gene_pathways(kegg_gene_id: str) -> list[str]:
    """
    Return pathway IDs linked to a KEGG gene ID.
    Example: ['path:hsa04110', 'path:hsa05200', ...]
    """
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
    """Fetch the pathway NAME line from KEGG for a given pathway id."""
    # pathway_id like 'path:hsa04110'
    pid = pathway_id.replace("path:", "")
    txt = kegg_get(f"/get/{pid}")
    for line in txt.split("\n"):
        if line.startswith("NAME"):
            return line.replace("NAME", "").strip()
    return None

# ----------------------------
# Core functions (vectorized-ish)
# ----------------------------
def fetch_gene_metadata_and_kegg(gene_list: list[str], organism_entrez: str, kegg_org_prefix: str, progress=None):
    """
    For each gene symbol -> (Gene, NCBI_ID, Description, KEGG_Pathways (pid - name; ...))
    Also returns a mapping: pathway_id -> set(genes) used later for enrichment.
    """
    results = []
    pathway_to_genes = defaultdict(set)

    for i, gene in enumerate(gene_list, start=1):
        # Show progress in UI
        if progress:
            progress.progress(min(i / max(len(gene_list), 1), 1.0))
        try:
            # 1) NCBI search
            ids = ncbi_esearch_gene_ids(gene, organism_entrez)
            if not ids:
                results.append({"Gene": gene, "NCBI_ID": None, "Description": "No match found", "KEGG_Pathways": None})
                continue
            gene_id = ids[0]

            # 2) NCBI description
            description = ncbi_esummary_description(gene_id)

            # 3) KEGG mapping
            kegg_id = kegg_ncbi_to_kegg_gene_id(gene_id, kegg_org_prefix)
            if not kegg_id:
                results.append({"Gene": gene, "NCBI_ID": gene_id, "Description": description, "KEGG_Pathways": None})
                continue

            # 4) KEGG pathways
            pids = kegg_gene_pathways(kegg_id)
            pairs = []
            for pid in pids:
                name = kegg_pathway_name(pid) or ""
                pairs.append(f"{pid} - {name}")
                pathway_to_genes[pid].add(gene)

            pathways_str = "; ".join(pairs) if pairs else None
            results.append({"Gene": gene, "NCBI_ID": gene_id, "Description": description, "KEGG_Pathways": pathways_str})

            # Be gentle to NCBI/KEGG servers
            time.sleep(0.20)
        except Exception as e:
            results.append({"Gene": gene, "NCBI_ID": None, "Description": f"Error: {e}", "KEGG_Pathways": None})

    return pd.DataFrame(results), pathway_to_genes

def hypergeom_pval(M: int, K: int, n: int, x: int) -> float:
    """
    Upper-tail (P[X >= x]) hypergeometric p-value.
    M: Population size (gene universe)
    K: # of population successes (genes associated with pathway)
    n: sample size (your input list)
    x: overlap (hits in pathway)
    """
    # sum_{i=x}^{min(K, n)} [ C(K,i) * C(M-K, n-i) ] / C(M, n)
    denom = math.comb(M, n) if 0 <= n <= M else 1
    s = 0.0
    upper = min(K, n)
    for i in range(x, upper + 1):
        s += math.comb(K, i) * math.comb(M - K, n - i)
    return s / denom if denom else 1.0

def compute_enrichment(pathway_to_genes: dict, gene_list: list[str], kegg_org_prefix: str, universe_size: int = 20000):
    """
    Build enrichment table:
    - Pathway_ID
    - Pathway_Name
    - Count (how many of your genes hit it)
    - Genes (semicolon-separated)
    - PValue (hypergeometric)
    Notes:
      * K (genes-in-pathway in universe) is approximated using KEGG /link/genes per pathway.
      * Universe size defaults to 20,000 (approximate number of protein-coding genes in human).
    """
    # Precompute KEGG pathway -> all KEGG genes size K (approximate coverage)
    K_cache = {}
    for pid in pathway_to_genes.keys():
        pid_clean = pid.replace("path:", "")
        try:
            txt = kegg_get(f"/link/genes/{pid_clean}")
            # lines like: path:hsa04110  hsa:1234
            genes_on_pathway = {line.split("\t")[1].strip() for line in txt.strip().split("\n") if "\t" in line}
            K_cache[pid] = len(genes_on_pathway)
        except Exception:
            K_cache[pid] = None

    n = len(gene_list)
    rows = []
    for pid, genes in sorted(pathway_to_genes.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        count = len(genes)
        pname = kegg_pathway_name(pid) or ""
        # K: pathway size in universe
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
# UI
# ----------------------------
st.title("🧬 BioContext – Contextual Gene Metadata & Pathway Enrichment")
st.markdown("Upload a gene list (CSV or TXT), select organism, and download annotated results + enrichment.")

# Entrez email (required by NCBI)
email = st.text_input("NCBI Entrez email (required)", value="", help="NCBI asks for a contact email for E-Utilities.")
if email:
    Entrez.email = email

# Organism selection
organisms = {
    "Homo sapiens (human)": {"entrez": "Homo sapiens", "kegg": "hsa"},
    "Mus musculus (mouse)": {"entrez": "Mus musculus", "kegg": "mmu"},
    "Rattus norvegicus (rat)": {"entrez": "Rattus norvegicus", "kegg": "rno"},
}
org_label = st.selectbox("Organism", list(organisms.keys()), index=0)
organism_entrez = organisms[org_label]["entrez"]
kegg_org_prefix = organisms[org_label]["kegg"]

# Universe size for enrichment
universe_size = st.number_input(
    "Gene universe size for enrichment (approx.)",
    min_value=1000, max_value=100000, value=20000, step=1000,
    help="Used for hypergeometric p-values. ~20,000 is a common default for human protein-coding genes."
)

uploaded = st.file_uploader("Upload gene list (.csv or .txt). If CSV, gene symbols must be in the first column.", type=["csv", "txt"])

run_btn = st.button("Run Annotation + Enrichment", type="primary", disabled=not uploaded or not email)

if run_btn:
    # Load genes
    try:
        ext = uploaded.name.lower().split(".")[-1]
        if ext == "csv":
            df_in = pd.read_csv(uploaded)
            genes = df_in.iloc[:, 0].dropna().astype(str).str.strip().str.upper().unique().tolist()
        else:
            text = uploaded.read().decode("utf-8", errors="ignore")
            genes = [line.strip().upper() for line in text.splitlines() if line.strip()]
        st.success(f"Loaded {len(genes)} genes.")
    except Exception as e:
        st.error(f"Could not read file: {e}")
        st.stop()

    # Progress UI
    st.subheader("Step 1 — Fetch NCBI + KEGG annotations")
    progress = st.progress(0.0)
    with st.spinner("Querying NCBI and KEGG..."):
        df_meta, pathway_to_genes = fetch_gene_metadata_and_kegg(
            genes, organism_entrez, kegg_org_prefix, progress=progress
        )
    st.success("Metadata retrieval complete.")
    st.dataframe(df_meta, use_container_width=True)

    # Download metadata CSV
    csv_bytes_meta = df_meta.to_csv(index=False).encode("utf-8")
    st.download_button(
        "Download metadata CSV",
        data=csv_bytes_meta,
        file_name="gene_metadata_with_kegg.csv",
        mime="text/csv"
    )

    # Enrichment
    st.subheader("Step 2 — Pathway Enrichment (KEGG)")
    with st.spinner("Computing enrichment..."):
        df_enrich = compute_enrichment(pathway_to_genes, genes, kegg_org_prefix, universe_size=universe_size)
    if df_enrich.empty:
        st.info("No pathways found for enrichment with the current gene list.")
    else:
        st.dataframe(df_enrich, use_container_width=True)
        csv_bytes_enrich = df_enrich.to_csv(index=False).encode("utf-8")
        st.download_button(
            "Download enrichment CSV",
            data=csv_bytes_enrich,
            file_name="pathway_enrichment.csv",
            mime="text/csv"
        )

st.markdown("---")
st.caption("Tip: If you run into rate limits, try again after a short pause. Requests are lightly throttled to be kind to the APIs.")
