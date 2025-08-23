# streamlit_app.py
# -------------------------------------------------------------
# BioContext – Gene2Therapy (UI-polished)
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
# App Config / Theming
# ----------------------------
st.set_page_config(
    page_title="Gene2Therapy – BioContext",
    page_icon="logoo.png",
    layout="wide",
)

# Global CSS – colors, cards, hero, chips, sticky bar
st.markdown(
    """
    <style>
      :root{
        --bg:#0b1220; --panel:#0f172a; --panel-2:#111a2c;
        --muted:#9fb0c1; --text:#e6edf3; --brand:#7dd3fc; --line:#1e293b;
        --ok:#10b981; --warn:#f59e0b; --bad:#ef4444;
      }
      .stApp{ background:var(--bg); color:var(--text); }
      .hero{
        background: linear-gradient(120deg,#0e1a33 0%, #0b1220 40%, #06233a 100%);
        border:1px solid var(--line); border-radius:20px; padding:22px 20px; margin-top:8px;
      }
      .hero h1{ margin:0; font-size:30px; letter-spacing:.4px }
      .hero .sub{ color:var(--muted); margin-top:6px }
      .pill{
        display:inline-flex; align-items:center; gap:8px;
        background:rgba(125,211,252,.12); border:1px solid rgba(125,211,252,.35);
        color:var(--brand); padding:4px 10px; border-radius:999px; font-size:12px;
      }
      .card{
        background:var(--panel); border:1px solid var(--line);
        border-radius:18px; padding:16px 16px 12px 16px; margin-bottom:14px;
      }
      .metric{
        background:var(--panel-2); border:1px solid var(--line);
        border-radius:16px; padding:14px; text-align:center;
      }
      .metric .t{ color:var(--muted); font-size:12px }
      .metric .v{ font-size:22px; font-weight:700; margin-top:4px }
      .chips span{
        background:var(--panel-2); border:1px solid var(--line);
        border-radius:999px; padding:5px 10px; margin:4px; display:inline-block;
        transition: transform .08s ease-in-out;
      }
      .chips span:hover{ transform: translateY(-1px); }
      .divider{ height:1px; background:var(--line); margin:18px 0 }
      .sticky{
        position: sticky; bottom: 8px; z-index: 9; padding:8px;
        background: linear-gradient(180deg, rgba(11,18,32,0) 0%, rgba(11,18,32,.9) 30%, rgba(11,18,32,1) 100%);
        backdrop-filter: blur(4px);
      }
      .help{ color:var(--muted); font-size:13px }
      .side-step{
        padding:8px 10px; border-radius:10px; border:1px solid var(--line); background:var(--panel);
        margin-bottom:6px; display:flex; gap:8px; align-items:center; font-size:13px;
      }
      .side-step .dot{ width:10px; height:10px; border-radius:50%; background:#334155; }
      .side-step.done .dot{ background:var(--ok); }
      .side-step span{ color:var(--muted) }
      .side-step.done span{ color:var(--text) }
    </style>
    """, unsafe_allow_html=True
)

# ---------------------------------
# Session helpers & small utilities
# ---------------------------------
def sset(k,v): st.session_state[k]=v
def sget(k,default=None): return st.session_state.get(k,default)

def init_state():
    sset("genes", sget("genes", []))
    sset("parsed_from", sget("parsed_from", None))
    sset("df_meta", None); sset("pathway_to_genes", None)
    sset("df_enrich", None); sset("g2t", None)
    sset("df_dis", None); sset("df_drugs", None)
    sset("step_done", sget("step_done", {"input":False,"meta":False,"enrich":False,"disease":False,"drugs":False}))

def show_gene_chips(genes, max_show=48):
    if not genes: return
    shown = genes[:max_show]
    chips = " ".join([f"<span>{g}</span>" for g in shown])
    extra = f"<span style='color:#9fb0c1'>&nbsp;+{len(genes)-len(shown)} more</span>" if len(genes)>len(shown) else ""
    st.markdown(f"<div class='chips'>{chips}{extra}</div>", unsafe_allow_html=True)

def bh_fdr(pvals: pd.Series)->pd.Series:
    s=pvals.copy().astype(float); n=s.notna().sum()
    if n==0: return s
    order=s.sort_values().index
    ranks=pd.Series(range(1,n+1), index=order, dtype=float)
    q=(s.loc[order]*n/ranks).cummin().clip(upper=1.0)
    out=pd.Series(index=s.index, dtype=float); out.loc[order]=q; return out

init_state()

# ----------------------------
# Parse genes from uploads/paste
# ----------------------------
def load_genes_from_any(uploaded_file)->list[str]:
    name=(uploaded_file.name or "").lower()
    def _clean(series: pd.Series)->list[str]:
        seen,out=set(),[]
        for v in (series.dropna().astype(str).str.strip().str.upper()):
            if v and v not in seen:
                seen.add(v); out.append(v)
            if len(out)>=200: break
        return out
    try:
        if name.endswith((".csv",".csv.gz")):
            df=pd.read_csv(uploaded_file, compression="infer")
        elif name.endswith((".tsv",".tsv.gz")):
            df=pd.read_csv(uploaded_file, sep="\t", compression="infer")
        elif name.endswith((".xlsx",".xls")):
            df=pd.read_excel(uploaded_file)
        else:
            df=None
        if isinstance(df,pd.DataFrame) and not df.empty:
            lower={str(c).lower():c for c in df.columns}
            col=lower.get("gene.symbol") or lower.get("symbol") or df.columns[0]
            return _clean(df[col])
    except Exception:
        pass
    try:
        raw=uploaded_file.read()
        text=raw.decode("utf-8","ignore") if isinstance(raw,(bytes,bytearray)) else str(raw)
        return _clean(pd.Series([ln for ln in (t.strip() for t in text.splitlines()) if ln]))
    except Exception:
        return []

# ----------------------------
# KEGG / NCBI (cached)
# ----------------------------
@st.cache_data(ttl=3600)
def kegg_get(path:str)->str:
    r=requests.get(f"https://rest.kegg.jp{path}", timeout=30)
    r.raise_for_status(); return r.text

@st.cache_data(ttl=3600)
def ncbi_esearch_gene_ids(gene_symbol:str, organism_entrez:str)->list[str]:
    handle=Entrez.esearch(db="gene", term=f"{gene_symbol}[Gene] AND {organism_entrez}[Organism]", retmode="xml")
    record=Entrez.read(handle); handle.close(); return record.get("IdList",[])

@st.cache_data(ttl=3600)
def ncbi_esummary_description(gene_id:str)->str:
    handle=Entrez.esummary(db="gene", id=gene_id, retmode="xml")
    raw=handle.read(); handle.close()
    root=ET.fromstring(raw); docsum=root.find(".//DocumentSummary")
    return (docsum.findtext("Description", default="") or "").strip()

@st.cache_data(ttl=3600)
def kegg_ncbi_to_kegg_gene_id(ncbi_gene_id:str, kegg_org_prefix:str)->str|None:
    txt=kegg_get(f"/conv/genes/ncbi-geneid:{ncbi_gene_id}")
    if not txt.strip(): return None
    for line in txt.strip().split("\n"):
        parts=line.split("\t")
        if len(parts)==2 and parts[0].endswith(f"{ncbi_gene_id}") and parts[1].startswith(f"{kegg_org_prefix}:"):
            return parts[1].strip()
    return None

@st.cache_data(ttl=3600)
def kegg_gene_pathways(kegg_gene_id:str)->list[str]:
    txt=kegg_get(f"/link/pathway/{kegg_gene_id}")
    if not txt.strip(): return []
    out=[]
    for line in txt.strip().split("\n"):
        parts=line.split("\t")
        if len(parts)==2 and parts[1].startswith("path:"):
            out.append(parts[1])
    return out

@st.cache_data(ttl=3600)
def kegg_pathway_name(pathway_id:str)->str|None:
    pid=pathway_id.replace("path:","")
    txt=kegg_get(f"/get/{pid}")
    for line in txt.split("\n"):
        if line.startswith("NAME"):
            return line.replace("NAME","").strip()
    return None

# ----------------------------
# OpenTargets (cached)
# ----------------------------
OT_GQL="https://api.platform.opentargets.org/api/v4/graphql"

@st.cache_data(ttl=3600)
def ot_query(query:str, variables:dict|None=None)->dict:
    try:
        r=requests.post(OT_GQL, json={"query":query,"variables":variables or {}}, timeout=40)
        data=r.json()
        if r.status_code>=400 or (isinstance(data,dict) and data.get("errors") and not data.get("data")):
            return {}
        return data if isinstance(data,dict) else {}
    except Exception:
        return {}

@st.cache_data(ttl=3600)
def ot_target_from_symbol(symbol:str, species:str="Homo sapiens")->dict|None:
    q="""
    query FindTarget($q: String!) {
      search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 5}) {
        hits { id name entity }
      }
    }"""
    data=ot_query(q,{"q":symbol})
    hits=(((data or {}).get("data",{})).get("search",{}) or {}).get("hits",[])
    for h in hits:
        if h.get("entity")=="target":
            return {"id":h.get("id"), "approvedSymbol":h.get("name")}
    return None

@st.cache_data(ttl=3600)
def ot_diseases_for_target(ensembl_id:str, size:int=25)->pd.DataFrame:
    q="""
    query Associations($id: String!, $size: Int!) {
      target(ensemblId: $id) {
        id
        associatedDiseases(page: {size: $size, index: 0}) {
          rows { disease { id name } score }
        }
      }
    }"""
    data=ot_query(q,{"id":ensembl_id,"size":size})
    rows=(((data or {}).get("data",{})).get("target",{}) or {}).get("associatedDiseases",{}).get("rows",[])
    out=[]
    for r in rows:
        d=r.get("disease",{})
        out.append({"target":ensembl_id,"disease_id":d.get("id"),"disease_name":d.get("name"),
                    "association_score":r.get("score")})
    return pd.DataFrame(out)

@st.cache_data(ttl=3600)
def ot_drugs_for_target(ensembl_id:str, size:int=50)->pd.DataFrame:
    q="""
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
    }"""
    data=ot_query(q,{"id":ensembl_id,"size":size})
    rows=(((data or {}).get("data",{})).get("target",{}) or {}).get("knownDrugs",{}).get("rows",[])
    out=[]
    for r in rows:
        drug=r.get("drug") or {}; disease=r.get("disease") or {}
        out.append({"target":ensembl_id,"drug_id":drug.get("id"),"drug_name":drug.get("name"),
                    "phase":r.get("phase"),"moa":r.get("mechanismOfAction"),
                    "diseases":"; ".join(filter(None,[disease.get("name")]))})
    return pd.DataFrame(out)

# ----------------------------
# Core logic
# ----------------------------
def fetch_gene_metadata_and_kegg(gene_list, organism_entrez, kegg_org_prefix, progress=None):
    results=[]; pathway_to_genes=defaultdict(set)
    for i,gene in enumerate(gene_list, start=1):
        if progress: progress.progress(min(i/max(len(gene_list),1),1.0), text=f"Fetching {gene}")
        try:
            ids=ncbi_esearch_gene_ids(gene, organism_entrez)
            if not ids:
                results.append({"Gene":gene,"NCBI_ID":None,"Description":"No match found","KEGG_Pathways":None}); continue
            gid=ids[0]; desc=ncbi_esummary_description(gid)
            kegg_id=kegg_ncbi_to_kegg_gene_id(gid, kegg_org_prefix)
            if not kegg_id:
                results.append({"Gene":gene,"NCBI_ID":gid,"Description":desc,"KEGG_Pathways":None}); continue
            pids=kegg_gene_pathways(kegg_id)
            pairs=[]
            for pid in pids:
                name=kegg_pathway_name(pid) or ""
                pairs.append(f"{pid} - {name}")
                pathway_to_genes[pid].add(gene)
            results.append({"Gene":gene,"NCBI_ID":gid,"Description":desc,
                            "KEGG_Pathways":"; ".join(pairs) if pairs else None})
            time.sleep(0.18)
        except Exception as e:
            results.append({"Gene":gene,"NCBI_ID":None,"Description":f"Error: {e}","KEGG_Pathways":None})
    return pd.DataFrame(results), pathway_to_genes

def hypergeom_pval(M,K,n,x)->float:
    denom=math.comb(M,n) if 0<=n<=M else 1
    s=0.0; upper=min(K,n)
    for i in range(x, upper+1):
        s += math.comb(K,i) * math.comb(M-K, n-i)
    return s/denom if denom else 1.0

def compute_enrichment(pathway_to_genes, gene_list, kegg_org_prefix, universe_size=20000):
    K_cache={}
    for pid in pathway_to_genes.keys():
        pidc=pid.replace("path:","")
        try:
            txt=kegg_get(f"/link/genes/{pidc}")
            genes_on_pathway={line.split("\t")[1].strip() for line in txt.strip().split("\n") if "\t" in line}
            K_cache[pid]=len(genes_on_pathway)
        except Exception:
            K_cache[pid]=None
    n=len(gene_list); rows=[]
    for pid, genes in sorted(pathway_to_genes.items(), key=lambda kv:(-len(kv[1]), kv[0])):
        cnt=len(genes); pname=kegg_pathway_name(pid) or ""; K=K_cache.get(pid,None); pval=None
        if K is not None and universe_size is not None:
            try: pval=hypergeom_pval(universe_size, K, n, cnt)
            except Exception: pval=None
        rows.append({"Pathway_ID":pid.replace("path:",""),"Pathway_Name":pname,"Count":cnt,
                     "Genes":";".join(sorted(genes)),"PValue":pval})
    df=pd.DataFrame(rows)
    if not df.empty: df=df.sort_values(["Count","PValue"], ascending=[False,True], na_position="last")
    return df

@st.cache_data(ttl=3600)
def build_gene_to_ot_target_map(genes, species="Homo sapiens")->dict:
    g2t={}
    for g in genes:
        hit=ot_target_from_symbol(g, species)
        if hit: g2t[g]=hit
        time.sleep(0.05)
    return g2t

@st.cache_data(ttl=3600)
def collect_disease_links(gene_to_target, size:int)->pd.DataFrame:
    frames=[]
    for g,tgt in gene_to_target.items():
        tid=tgt.get("id"); if not tid: continue
        df=ot_diseases_for_target(tid, size=size)
        if not df.empty: df.insert(0,"gene",g); frames.append(df)
        time.sleep(0.05)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
        columns=["gene","target","disease_id","disease_name","association_score"]
    )

@st.cache_data(ttl=3600)
def collect_drug_suggestions(gene_to_target, size:int)->pd.DataFrame:
    frames=[]
    for g,tgt in gene_to_target.items():
        tid=tgt.get("id"); if not tid: continue
        df=ot_drugs_for_target(tid, size=size)
        if not df.empty: df.insert(0,"gene",g); frames.append(df)
        time.sleep(0.05)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
        columns=["gene","target","drug_id","drug_name","phase","moa","diseases"]
    )

# ----------------------------
# HERO + SIDEBAR PROGRESS
# ----------------------------
with st.container():
    colh1, colh2 = st.columns([7,3])
    with colh1:
        st.markdown(
            """
            <div class='hero'>
              <div class='pill'>🧬 BioContext • Gene2Therapy</div>
              <h1>From genes to pathways, diseases & repurposable drugs</h1>
              <div class='sub'>Upload or paste a gene list (≤200 symbols). Live annotation via NCBI, KEGG & OpenTargets with 1-hour caching.</div>
            </div>
            """, unsafe_allow_html=True
        )
    with colh2:
        st.image("logoo.png", width=120)

with st.sidebar:
    st.markdown("#### Progress")
    steps = [("Input", "input"), ("Metadata", "meta"), ("Enrichment", "enrich"),
             ("Diseases", "disease"), ("Drugs", "drugs")]
    for label,key in steps:
        cls="side-step done" if sget("step_done")[key] else "side-step"
        st.markdown(f"<div class='{cls}'><div class='dot'></div><span>{label}</span></div>", unsafe_allow_html=True)
    st.markdown("---")
    st.caption("Tips\n- Keep lists modest (≤300) to avoid API throttling.\n- Re-run if rate-limited; results are cached ~1h.")

# ----------------------------
# INPUT CARD
# ----------------------------
with st.container():
    st.markdown("### 1) Input")
    with st.container():
        c1,c2 = st.columns([3,2])
        with c1:
            email = st.text_input("NCBI Entrez email", value="", placeholder="name@institute.org",
                                  help="NCBI requests a contact email for E-Utilities.")
            if email: Entrez.email=email
        with c2:
            density = st.toggle("Compact tables", value=False)

    orgs = {
        "Homo sapiens (human)": {"entrez":"Homo sapiens","kegg":"hsa"},
        "Mus musculus (mouse)": {"entrez":"Mus musculus","kegg":"mmu"},
        "Rattus norvegicus (rat)": {"entrez":"Rattus norvegicus","kegg":"rno"},
    }
    org_label = st.selectbox("Organism", list(orgs.keys()), index=0)
    organism_entrez = orgs[org_label]["entrez"]; kegg_org_prefix = orgs[org_label]["kegg"]

    with st.expander("Advanced options", expanded=False):
        universe_size = st.number_input("Gene universe size for enrichment", 1000, 100000, 20000, 1000)
        ot_size_diseases = st.slider("Max diseases per target (OpenTargets)", 10, 200, 50, 10)
        ot_size_drugs    = st.slider("Max drugs per target (OpenTargets)", 10, 200, 100, 10)

    st.markdown("<div class='divider'></div>", unsafe_allow_html=True)

    src = st.radio("Choose input source", ["Upload a file","Paste genes","Use demo list"], horizontal=True)
    genes_from_input=[]
    if src=="Upload a file":
        uploaded = st.file_uploader("Upload CSV/TSV/XLSX/TXT (I’ll auto-detect the gene column)", type=["csv","tsv","txt","xlsx"])
        if uploaded: genes_from_input = load_genes_from_any(uploaded); sset("parsed_from","upload")
    elif src=="Paste genes":
        manual = st.text_area("Paste gene symbols (comma/space/newline separated)", height=120,
                              placeholder="TP53, BRCA1, EGFR, MYC, ...")
        if manual.strip():
            raw = manual.replace(",","\n").replace(" ","\n")
            seen=set(); cleaned=[]
            for g in raw.splitlines():
                gg=g.strip().upper()
                if gg and gg not in seen:
                    seen.add(gg); cleaned.append(gg)
                if len(cleaned)>=200: break
            genes_from_input=cleaned; sset("parsed_from","paste")
    else:
        demo=["TP53","BRCA1","EGFR","MYC","PIK3CA","PTEN","KRAS","BRAF","CDK4","CDK6","CDKN2A","ERBB2","NRAS","ALK","ROS1"]
        genes_from_input=demo; sset("parsed_from","demo"); st.info("Demo list loaded.")

    if genes_from_input:
        sset("genes", genes_from_input)
        sget("step_done")["input"]=True
        show_gene_chips(genes_from_input)

# Sticky action bar
st.markdown("<div class='sticky'>", unsafe_allow_html=True)
colact1,colact2,colact3 = st.columns([1,1,8])
with colact1:
    run_btn = st.button("▶ Analyze", type="primary", disabled=(not sget("genes") or not email))
with colact2:
    if st.button("🧹 Clear cache"):
        st.cache_data.clear(); st.toast("Cache cleared.", icon="🧹")
st.markdown("</div>", unsafe_allow_html=True)

# ----------------------------
# TABS
# ----------------------------
meta_tab, enrich_tab, disease_tab, drug_tab, viz_tab = st.tabs(
    ["🗂️ Metadata","📊 Enrichment","🩺 Disease Links","💊 Drug Suggestions","🧠 Visualize"]
)

# ----------------------------
# PIPELINE
# ----------------------------
if run_btn:
    genes = sget("genes")

    # Step 1: Metadata
    with meta_tab:
        st.subheader("NCBI + KEGG annotations")
        progress = st.progress(0.0)
        with st.spinner("Querying NCBI and KEGG..."):
            df_meta, pathway_to_genes = fetch_gene_metadata_and_kegg(genes, organism_entrez, kegg_org_prefix, progress)
        sset("df_meta", df_meta); sset("pathway_to_genes", pathway_to_genes)
        sget("step_done")["meta"]=True

        # Metrics
        n_genes = len(genes)
        n_kegg = sum(df_meta["KEGG_Pathways"].notna()) if not df_meta.empty else 0
        c1,c2,c3 = st.columns(3)
        c1.markdown(f"<div class='metric'><div class='t'>Genes</div><div class='v'>{n_genes}</div></div>", unsafe_allow_html=True)
        c2.markdown(f"<div class='metric'><div class='t'>Genes with KEGG</div><div class='v'>{n_kegg}</div></div>", unsafe_allow_html=True)
        c3.markdown(f"<div class='metric'><div class='t'>Unique pathways (raw)</div><div class='v'>{len(pathway_to_genes)}</div></div>", unsafe_allow_html=True)

        st.dataframe(df_meta, use_container_width=True, hide_index=True,
                     height=420 if not density else 300)
        st.download_button("⬇️ Download metadata CSV", df_meta.to_csv(index=False).encode("utf-8"),
                           file_name="gene_metadata_with_kegg.csv", mime="text/csv")

    # Step 2: Enrichment
    with enrich_tab:
        st.subheader("Pathway Enrichment (KEGG)")
        with st.spinner("Computing enrichment..."):
            df_enrich = compute_enrichment(sget("pathway_to_genes"), genes, kegg_org_prefix, universe_size=universe_size)
        sset("df_enrich", df_enrich); sget("step_done")["enrich"]=True

        if df_enrich.empty:
            st.info("No pathways found for enrichment with the current gene list.")
        else:
            if "PValue" in df_enrich.columns:
                df_enrich["q_value"] = bh_fdr(df_enrich["PValue"])

            f1,f2,f3 = st.columns([1,1,2])
            with f1: min_count = st.number_input("Min genes in pathway", 1, 100, 2)
            with f2: q_cut = st.selectbox("q-value max", ["none",0.05,0.10,0.20], index=1)
            with f3: search = st.text_input("Search pathway name", "")

            dfv = df_enrich.copy()
            dfv = dfv[dfv["Count"] >= min_count]
            if q_cut!="none" and "q_value" in dfv: dfv = dfv[dfv["q_value"].le(float(q_cut))]
            if search: dfv = dfv[dfv["Pathway_Name"].str.contains(search, case=False, na=False)]

            st.dataframe(dfv, use_container_width=True, hide_index=True,
                         height=460 if not density else 320,
                         column_config={"PValue":st.column_config.NumberColumn(format="%.3e"),
                                        "q_value":st.column_config.NumberColumn(format="%.3e")})
            st.download_button("⬇️ Download enrichment (filtered)",
                               dfv.to_csv(index=False).encode("utf-8"),
                               file_name="pathway_enrichment_filtered.csv", mime="text/csv")
            try:
                topN = dfv.head(15).copy()
                fig = px.bar(topN, x="Count", y="Pathway_Name", orientation="h",
                             title="Top enriched pathways (filtered)")
                fig.update_layout(height=560)
                st.plotly_chart(fig, use_container_width=True)
            except Exception:
                pass

    # Step 3: Diseases
    with disease_tab:
        st.subheader("Disease Associations (OpenTargets)")
        with st.spinner("Mapping symbols and fetching disease links..."):
            g2t = build_gene_to_ot_target_map(genes, species="Homo sapiens")
            df_dis = collect_disease_links(g2t, size=ot_size_diseases)
        sset("g2t", g2t); sset("df_dis", df_dis); sget("step_done")["disease"]=True

        if df_dis.empty:
            st.info("No disease associations retrieved (try human genes or a smaller list).")
        else:
            agg = (df_dis.groupby(["disease_id","disease_name"])
                        .agg(n_genes=("gene", lambda s: len(set(s))),
                             max_score=("association_score","max"))
                        .reset_index()
                        .sort_values(["n_genes","max_score"], ascending=[False,False]))
            topN_dis = st.slider("Show top N diseases", 10, 200, 50, 10)
            min_genes_dis = st.number_input("Min #genes per disease", 1, 100, 2)
            search_dis = st.text_input("Search disease name", "")

            aggv = agg[agg["n_genes"]>=min_genes_dis].copy()
            if search_dis: aggv = aggv[aggv["disease_name"].str.contains(search_dis, case=False, na=False)]

            # Metrics
            c1,c2 = st.columns(2)
            c1.markdown(f"<div class='metric'><div class='t'>Diseases found</div><div class='v'>{len(agg)}</div></div>", unsafe_allow_html=True)
            genes_w_hits = len(set(df_dis['gene'])) if not df_dis.empty else 0
            c2.markdown(f"<div class='metric'><div class='t'>Genes with disease hits</div><div class='v'>{genes_w_hits}</div></div>", unsafe_allow_html=True)

            st.dataframe(aggv.head(topN_dis), use_container_width=True, hide_index=True,
                         height=460 if not density else 320)
            cdl1, cdl2 = st.columns(2)
            cdl1.download_button("⬇️ Per-gene disease links",
                                 df_dis.to_csv(index=False).encode("utf-8"),
                                 file_name="gene_disease_links_opentargets.csv", mime="text/csv")
            cdl2.download_button("⬇️ Disease summary (aggregated)",
                                 agg.to_csv(index=False).encode("utf-8"),
                                 file_name="disease_summary_aggregated.csv", mime="text/csv")
            try:
                topD = aggv.head(min(topN_dis, 20))
                figd = px.bar(topD, x="n_genes", y="disease_name", orientation="h",
                              title="Top disease associations (by #genes)")
                figd.update_layout(height=560)
                st.plotly_chart(figd, use_container_width=True)
            except Exception:
                pass

    # Step 4: Drugs
    with drug_tab:
        st.subheader("Repurposable Drug Suggestions (OpenTargets known drugs)")
        with st.spinner("Fetching known drugs targeting your genes..."):
            df_drugs = collect_drug_suggestions(sget("g2t"), size=ot_size_drugs)
        sset("df_drugs", df_drugs); sget("step_done")["drugs"]=True

        if df_drugs.empty:
            st.info("No drugs found for the mapped targets.")
        else:
            phase_rank = {None:0,1:1,2:2,3:3,4:4}
            df_drugs["phase_rank"]=df_drugs["phase"].map(phase_rank).fillna(0)
            drug_sum = (df_drugs.groupby(["drug_id","drug_name"]).agg(
                targets=("target", lambda s:";".join(sorted(set(s)))),
                genes=("gene", lambda s:";".join(sorted(set(s)))),
                max_phase=("phase","max"),
                max_phase_rank=("phase_rank","max"),
                indications=("diseases", lambda s:"; ".join(sorted({x for x in "; ".join(s).split("; ") if x}))),
                moa=("moa", lambda s:"; ".join(sorted({x for x in s if x})))
            ).reset_index().sort_values(["max_phase_rank","drug_name"], ascending=[False,True]))

            # Metrics
            n_drugs=len(drug_sum)
            st.markdown(f"<div class='metric'><div class='t'>Unique drugs</div><div class='v'>{n_drugs}</div></div>", unsafe_allow_html=True)

            min_phase = st.selectbox("Minimum clinical phase", [0,1,2,3,4], index=3)
            search_drug = st.text_input("Search drug/target/gene", "")
            dv = drug_sum[drug_sum["max_phase"].fillna(0).astype(int)>=int(min_phase)].copy()
            if search_drug:
                mask=(dv["drug_name"].str.contains(search_drug, case=False, na=False) |
                      dv["targets"].str.contains(search_drug, case=False, na=False) |
                      dv["genes"].str.contains(search_drug, case=False, na=False))
                dv=dv[mask]

            st.dataframe(dv, use_container_width=True, hide_index=True,
                         height=460 if not density else 320)
            cdl1, cdl2 = st.columns(2)
            cdl1.download_button("⬇️ Drug suggestions (per target)",
                                 df_drugs.to_csv(index=False).encode("utf-8"),
                                 file_name="drug_suggestions_per_target.csv", mime="text/csv")
            cdl2.download_button("⬇️ Drug suggestions (aggregated)",
                                 drug_sum.to_csv(index=False).encode("utf-8"),
                                 file_name="drug_suggestions_aggregated.csv", mime="text/csv")

    # Step 5: Visualize
    with viz_tab:
        st.subheader("Landscape")
        colA, colB = st.columns(2)

        with colA:
            try:
                df_dis = sget("df_dis")
                df_drugs = sget("df_drugs")
                if df_dis is not None and not df_dis.empty:
                    aggD = (df_dis.groupby("disease_name")
                                  .agg(n_genes=("gene", lambda s: len(set(s)))).reset_index()
                                  .sort_values("n_genes", ascending=False).head(10))
                    top_dis=set(aggD["disease_name"].tolist())
                    genes_set=sorted(set(df_dis[df_dis["disease_name"].isin(top_dis)]["gene"]))
                    dis_list=sorted(top_dis)
                    drugs_set=[]
                    if df_drugs is not None and not df_drugs.empty:
                        tmp=df_drugs[df_drugs['gene'].isin(genes_set)].copy()
                        tmp['phase_rank']=tmp['phase'].map({None:0,1:1,2:2,3:3,4:4}).fillna(0)
                        tmp=tmp.sort_values('phase_rank', ascending=False)
                        drugs_set=sorted(set(tmp.head(100)['drug_name']))[:15]
                    nodes=[*(f"G: {g}" for g in genes_set), *(f"D: {d}" for d in dis_list), *(f"Rx: {d}" for d in drugs_set)]
                    node_index={n:i for i,n in enumerate(nodes)}
                    links=[]
                    for d in dis_list:
                        sub=df_dis[df_dis["disease_name"]==d]
                        for g, cnt in Counter(sub["gene"]).items():
                            links.append((node_index[f"G: {g}"], node_index[f"D: {d}"], max(cnt,1)))
                    if drugs_set and df_drugs is not None and not df_drugs.empty:
                        tmp=df_drugs[df_drugs['drug_name'].isin(drugs_set) & df_drugs['gene'].isin(genes_set)]
                        for _, row in tmp.iterrows():
                            s=node_index[f"G: {row['gene']}"]; t=node_index.get(f"Rx: {row['drug_name']}"); 
                            if t is not None: links.append((s,t,max((row['phase'] or 0),1)))
                    if links:
                        sources=[s for s,_,_ in links]; targets=[t for _,t,_ in links]; values=[v for *_,v in links]
                        fig=go.Figure(data=[go.Sankey(node=dict(pad=12, thickness=14, label=nodes),
                                                      link=dict(source=sources, target=targets, value=values))])
                        fig.update_layout(title_text="Gene → Disease (→ Drug) connections", height=700)
                        st.plotly_chart(fig, use_container_width=True)
            except Exception as e:
                st.warning(f"Sankey could not be drawn: {e}")

        with colB:
            try:
                df_enrich = sget("df_enrich")
                if df_enrich is not None and not df_enrich.empty:
                    top_paths=df_enrich.head(8).copy()
                    edges=[]
                    for _, r in top_paths.iterrows():
                        p=r["Pathway_Name"] or r["Pathway_ID"]
                        for g in (r["Genes"] or "").split(";"):
                            if g: edges.append((p,g))
                    if edges:
                        G=nx.Graph(); G.add_edges_from(edges)
                        pos=nx.spring_layout(G, seed=42, k=0.7)
                        x_nodes=[pos[n][0] for n in G.nodes()]; y_nodes=[pos[n][1] for n in G.nodes()]
                        node_text=list(G.nodes())
                        xe,ye=[],[]
                        for u,v in G.edges():
                            xe += [pos[u][0], pos[v][0], None]; ye += [pos[u][1], pos[v][1], None]
                        fig=go.Figure()
                        fig.add_trace(go.Scatter(x=xe,y=ye,mode='lines',opacity=0.5))
                        fig.add_trace(go.Scatter(x=x_nodes,y=y_nodes,mode='markers+text',
                                                 text=node_text,textposition='top center'))
                        fig.update_layout(title="Gene–Pathway network (top enriched)", height=700, showlegend=False)
                        st.plotly_chart(fig, use_container_width=True)
            except Exception as e:
                st.warning(f"Network could not be drawn: {e}")

# Footer
st.markdown("<div class='divider'></div>", unsafe_allow_html=True)
st.caption("APIs: NCBI E-utilities, KEGG REST, OpenTargets GraphQL • Results cached for 1h • Validate with primary sources for clinical decisions.")
