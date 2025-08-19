import os
import time
import argparse
import pandas as pd
import requests
from Bio import Entrez
from urllib.error import HTTPError
from xml.etree import ElementTree as ET
from tqdm import tqdm
from scipy.stats import fisher_exact
from collections import Counter

# ---------- CONFIG ----------
Entrez.email = "your_real_email@example.com"  # replace with your real email
DATA_DIR = "data"
os.makedirs(DATA_DIR, exist_ok=True)


def load_gene_list(filepath):
    """Load gene list from CSV or TXT."""
    ext = os.path.splitext(filepath)[-1].lower()
    if ext == ".csv":
        df = pd.read_csv(filepath)
        genes = df.iloc[:, 0].dropna().unique().tolist()
    elif ext == ".txt":
        with open(filepath) as f:
            genes = [line.strip() for line in f if line.strip()]
    else:
        raise ValueError("File format not supported. Use CSV or TXT.")
    return [g.upper() for g in genes]


def fetch_gene_metadata(gene, organism="Homo sapiens"):
    """Fetch gene metadata from NCBI Entrez."""
    try:
        search_handle = Entrez.esearch(
            db="gene",
            term=f"{gene}[Gene] AND {organism}[Organism]",
            retmode="xml"
        )
        search_results = Entrez.read(search_handle)
        search_handle.close()

        if not search_results["IdList"]:
            return {"Gene": gene, "NCBI_ID": None, "Description": "No match found", "KEGG_Pathways": None}

        gene_id = search_results["IdList"][0]

        summary_handle = Entrez.esummary(db="gene", id=gene_id, retmode="xml")
        raw_xml = summary_handle.read()
        summary_handle.close()

        root = ET.fromstring(raw_xml)
        docsum = root.find(".//DocumentSummary")
        description = docsum.findtext("Description", default="").strip()

        # Fetch KEGG pathways
        pathways = fetch_kegg_pathways(gene_id)

        return {"Gene": gene, "NCBI_ID": gene_id, "Description": description, "KEGG_Pathways": pathways}

    except HTTPError as e:
        return {"Gene": gene, "NCBI_ID": None, "Description": f"HTTPError: {e.code}", "KEGG_Pathways": None}
    except Exception as e:
        return {"Gene": gene, "NCBI_ID": None, "Description": f"Error: {str(e)}", "KEGG_Pathways": None}


def fetch_kegg_pathways(ncbi_gene_id):
    """Fetch KEGG pathways using NCBI Gene ID (human-specific)."""
    try:
        # Step 1: Convert NCBI Gene ID → KEGG human gene ID
        r = requests.get(f"https://rest.kegg.jp/conv/genes/ncbi-geneid:{ncbi_gene_id}")
        if r.status_code != 200 or not r.text.strip():
            return None
        kegg_gene_id = r.text.split("\t")[1].strip()  # e.g., hsa:7157

        # Step 2: Get linked pathways
        r = requests.get(f"https://rest.kegg.jp/link/pathway/{kegg_gene_id}")
        if r.status_code != 200 or not r.text.strip():
            return None
        pathway_ids = [line.split("\t")[1] for line in r.text.strip().split("\n")]

        # Step 3: Get pathway names
        pathway_list = []
        for pid in pathway_ids:
            r = requests.get(f"https://rest.kegg.jp/get/{pid}")
            if r.status_code == 200:
                for line in r.text.split("\n"):
                    if line.startswith("NAME"):
                        pathway_list.append(f"{pid} - {line.replace('NAME', '').strip()}")
                        break

        return "; ".join(pathway_list) if pathway_list else None
    except:
        return None


def get_kegg_pathway_name(pathway_id):
    """Fetch KEGG pathway name given its ID (e.g., hsa04151)."""
    try:
        r = requests.get(f"https://rest.kegg.jp/get/{pathway_id}")
        if r.status_code == 200:
            for line in r.text.split("\n"):
                if line.startswith("NAME"):
                    return line.replace("NAME", "").strip()
    except:
        return None
    return None


def pathway_enrichment(results, background_size=20000):
    """Perform simple pathway enrichment using Fisher's exact test, with pathway names."""
    all_pathways = []
    for r in results:
        if r["KEGG_Pathways"]:
            all_pathways.extend([p.split(" - ")[0] for p in r["KEGG_Pathways"].split("; ")])

    pathway_counts = Counter(all_pathways)
    enrichment_results = []

    for pathway, count in pathway_counts.items():
        in_list_in_pathway = count
        in_list_not_in_pathway = len(results) - count
        background_in_pathway = 200  # placeholder
        background_not_in_pathway = background_size - background_in_pathway

        table = [[in_list_in_pathway, in_list_not_in_pathway],
                 [background_in_pathway, background_not_in_pathway]]
        odds, pval = fisher_exact(table, alternative="greater")

        # Fetch pathway name
        pathway_name = get_kegg_pathway_name(pathway)

        enrichment_results.append({
            "Pathway_ID": pathway,
            "Pathway_Name": pathway_name,
            "Count": count,
            "PValue": pval
        })

    df = pd.DataFrame(enrichment_results).sort_values("PValue")
    out_file = os.path.join(DATA_DIR, "pathway_enrichment.csv")
    df.to_csv(out_file, index=False)
    print(f"Pathway enrichment saved to {out_file}")


def main():
    parser = argparse.ArgumentParser(description="BioContext - Gene metadata retriever with KEGG enrichment")
    parser.add_argument("input", help="Input gene list (.txt or .csv)")
    parser.add_argument("--organism", default="Homo sapiens", help="Organism name (default: Homo sapiens)")
    parser.add_argument("--output", default="gene_metadata_with_kegg.csv", help="Output filename")
    args = parser.parse_args()

    genes = load_gene_list(args.input)

    print(f"Loaded {len(genes)} genes. Fetching metadata + KEGG pathways...")
    results = []
    for g in tqdm(genes, desc="Fetching genes"):
        results.append(fetch_gene_metadata(g, organism=args.organism))
        time.sleep(0.34)  # NCBI rate limit

    df = pd.DataFrame(results)
    output_path = os.path.join(DATA_DIR, args.output)
    df.to_csv(output_path, index=False)
    print(f"Metadata + KEGG pathways saved to {output_path}")

    # Run enrichment
    pathway_enrichment(results)


if __name__ == "__main__":
    main()
