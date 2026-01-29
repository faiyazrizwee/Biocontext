"""
Annotation App Views
Handles gene mapping, pathway enrichment, GO analysis, and disease/drug associations.
"""
import json
import requests
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
import time

class GeneMapView(APIView):
    """
    Map gene symbols to official IDs using NCBI Gene API.
    POST: {"genes": ["TP53", "BRCA1", "EGFR"]}
    """
    def post(self, request):
        genes = request.data.get('genes', [])
        if not genes:
            return Response({'error': 'No genes provided'}, status=status.HTTP_400_BAD_REQUEST)
        
        results = []
        for gene in genes[:50]:  # Limit to 50 genes per request
            try:
                # NCBI Entrez E-utilities
                url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                params = {
                    'db': 'gene',
                    'term': f'{gene}[Gene Name] AND Homo sapiens[Organism]',
                    'retmode': 'json'
                }
                response = requests.get(url, params=params, timeout=10)
                data = response.json()
                
                gene_ids = data.get('esearchresult', {}).get('idlist', [])
                results.append({
                    'symbol': gene,
                    'ncbi_gene_ids': gene_ids[:5],  # Top 5 matches
                    'status': 'found' if gene_ids else 'not_found'
                })
                time.sleep(0.35)  # NCBI rate limit: 3 requests/second
            except Exception as e:
                results.append({
                    'symbol': gene,
                    'ncbi_gene_ids': [],
                    'status': 'error',
                    'error': str(e)
                })
        
        return Response({'results': results})


class KEGGEnrichmentView(APIView):
    """
    Perform KEGG pathway enrichment analysis with pathway names.
    POST: {"genes": ["TP53", "BRCA1", "EGFR"]}
    """
    def post(self, request):
        genes = request.data.get('genes', [])
        if not genes:
            return Response({'error': 'No genes provided'}, status=status.HTTP_400_BAD_REQUEST)
        
        pathways = []
        pathway_names_cache = {}  # Cache pathway names
        
        try:
            # Get KEGG gene IDs for human genes
            for gene in genes[:15]:  # Limit to avoid timeout
                url = f"https://rest.kegg.jp/find/genes/{gene}"
                response = requests.get(url, timeout=10)
                if response.status_code == 200:
                    lines = response.text.strip().split('\n')
                    for line in lines:
                        if 'hsa:' in line:  # Human genes
                            kegg_id = line.split('\t')[0]
                            # Get pathways for this gene
                            pathway_url = f"https://rest.kegg.jp/link/pathway/{kegg_id}"
                            pw_response = requests.get(pathway_url, timeout=10)
                            if pw_response.status_code == 200:
                                for pw_line in pw_response.text.strip().split('\n'):
                                    if pw_line:
                                        parts = pw_line.split('\t')
                                        if len(parts) >= 2:
                                            pathway_id = parts[1]
                                            # Only include human pathways (hsa)
                                            if pathway_id.startswith('path:hsa'):
                                                pathways.append({
                                                    'gene': gene,
                                                    'kegg_id': kegg_id,
                                                    'pathway_id': pathway_id
                                                })
                            break
                time.sleep(0.4)  # Rate limiting
            
            # Get pathway names for unique pathways
            unique_pathways = list(set(p['pathway_id'] for p in pathways))
            for pw_id in unique_pathways[:30]:  # Limit name lookups
                if pw_id not in pathway_names_cache:
                    try:
                        name_url = f"https://rest.kegg.jp/get/{pw_id.replace('path:', '')}"
                        name_response = requests.get(name_url, timeout=8)
                        if name_response.status_code == 200:
                            # Parse the first NAME line
                            for line in name_response.text.split('\n'):
                                if line.startswith('NAME'):
                                    pathway_names_cache[pw_id] = line.replace('NAME', '').strip()
                                    break
                        time.sleep(0.3)
                    except:
                        pathway_names_cache[pw_id] = pw_id
                        
        except Exception as e:
            return Response({
                'pathways': [],
                'error': str(e)
            })
        
        # Aggregate pathways with names
        pathway_counts = {}
        for p in pathways:
            pw_id = p['pathway_id']
            pw_name = pathway_names_cache.get(pw_id, pw_id)
            
            if pw_id not in pathway_counts:
                pathway_counts[pw_id] = {
                    'pathway_id': pw_id,
                    'pathway': pw_name,
                    'genes': [],
                    'count': 0
                }
            if p['gene'] not in pathway_counts[pw_id]['genes']:
                pathway_counts[pw_id]['genes'].append(p['gene'])
                pathway_counts[pw_id]['count'] += 1
        
        enriched = sorted(pathway_counts.values(), key=lambda x: x['count'], reverse=True)
        return Response({'pathways': enriched[:25]})


class GOEnrichmentView(APIView):
    """
    Perform Gene Ontology enrichment using UniProt API.
    POST: {"genes": ["TP53", "BRCA1"]}
    """
    def post(self, request):
        genes = request.data.get('genes', [])
        if not genes:
            return Response({'error': 'No genes provided'}, status=status.HTTP_400_BAD_REQUEST)
        
        go_terms = []
        for gene in genes[:15]:
            try:
                # Use UniProt to map gene to protein, then get GO terms
                url = f"https://rest.uniprot.org/uniprotkb/search"
                params = {
                    'query': f'gene:{gene} AND organism_id:9606',
                    'format': 'json',
                    'size': 1
                }
                response = requests.get(url, params=params, timeout=10)
                if response.status_code == 200:
                    data = response.json()
                    results = data.get('results', [])
                    if results:
                        entry = results[0]
                        # Extract GO annotations
                        xrefs = entry.get('uniProtKBCrossReferences', [])
                        for xref in xrefs:
                            if xref.get('database') == 'GO':
                                go_id = xref.get('id', '')
                                props = xref.get('properties', [])
                                term = next((p['value'] for p in props if p.get('key') == 'GoTerm'), '')
                                go_terms.append({
                                    'gene': gene,
                                    'go_id': go_id,
                                    'term': term
                                })
                time.sleep(0.4)
            except Exception as e:
                go_terms.append({'gene': gene, 'error': str(e)})
        
        return Response({'go_terms': go_terms})


class DiseaseAssociationView(APIView):
    """
    Get disease associations using Open Targets Platform API.
    POST: {"genes": ["TP53", "BRCA1"]}
    """
    def post(self, request):
        genes = request.data.get('genes', [])
        if not genes:
            return Response({'error': 'No genes provided'}, status=status.HTTP_400_BAD_REQUEST)
        
        associations = []
        url = "https://api.platform.opentargets.org/api/v4/graphql"
        headers = {'Content-Type': 'application/json'}
        
        for gene in genes[:10]:
            try:
                # Step 1: Search for Ensembl ID
                search_query = """
                query searchTarget($q: String!) {
                    search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 1}) {
                        hits { id }
                    }
                }
                """
                search_res = requests.post(url, json={'query': search_query, 'variables': {'q': gene}}, timeout=10, headers=headers)
                if search_res.status_code != 200: 
                    print(f"Search failed for {gene}: {search_res.status_code}")
                    continue
                
                search_data = search_res.json()
                hits = search_data.get('data', {}).get('search', {}).get('hits', [])
                if not hits: 
                    print(f"No hits for gene: {gene}")
                    continue
                target_id = hits[0]['id']

                # Step 2: Get associations by Ensembl ID
                assoc_query = """
                query targetAssoc($id: String!) {
                    target(ensemblId: $id) {
                        associatedDiseases(page: {index: 0, size: 10}) {
                            rows {
                                disease { id name }
                                score
                            }
                        }
                    }
                }
                """
                assoc_res = requests.post(url, json={'query': assoc_query, 'variables': {'id': target_id}}, timeout=10, headers=headers)
                
                if assoc_res.status_code == 200:
                    data = assoc_res.json()
                    target_data = data.get('data', {}).get('target', {})
                    if not target_data: 
                        print(f"No target data for ID: {target_id}")
                        continue
                    
                    diseases = target_data.get('associatedDiseases', {}).get('rows', [])
                    print(f"Found {len(diseases)} diseases for {gene}")
                    for d in diseases[:5]:
                        disease_info = d.get('disease', {})
                        associations.append({
                            'gene': gene,
                            'disease': disease_info.get('name', 'Unknown'),
                            'disease_id': disease_info.get('id', ''),
                            'score': round(d.get('score', 0), 4)
                        })
                else:
                    print(f"Assoc query failed for {target_id}: {assoc_res.status_code}")
                time.sleep(0.2)
                
            except Exception as e:
                print(f"Error fetching diseases for {gene}: {str(e)}")
        
        # Remove duplicates
        seen = set()
        unique_associations = []
        for a in associations:
            key = f"{a['gene']}_{a['disease']}"
            if key not in seen:
                seen.add(key)
                unique_associations.append(a)
        
        unique_associations.sort(key=lambda x: x['score'], reverse=True)
        return Response({'associations': unique_associations})


class DrugAssociationView(APIView):
    """
    Get drug-gene interactions using DGIdb API and Open Targets.
    POST: {"genes": ["TP53", "BRCA1"]}
    """
    def post(self, request):
        genes = request.data.get('genes', [])
        if not genes:
            return Response({'error': 'No genes provided'}, status=status.HTTP_400_BAD_REQUEST)
        
        drugs = []
        
        # Try DGIdb first with better error handling
        try:
            gene_list = ','.join(genes[:20])
            url = "https://www.dgidb.org/api/v2/interactions.json"
            params = {'genes': gene_list}
            
            response = requests.get(url, params=params, timeout=15)
            
            if response.status_code == 200:
                try:
                    data = response.json()
                    matched_terms = data.get('matchedTerms', [])
                    for match in matched_terms:
                        gene_name = match.get('geneName', match.get('searchTerm', ''))
                        interactions = match.get('interactions', [])
                        for interaction in interactions:
                            drug_name = interaction.get('drugName', interaction.get('drugConceptName', ''))
                            drugs.append({
                                'gene': gene_name,
                                'drug_name': drug_name,
                                'interaction_types': ['Interaction'],
                                'source': 'DGIdb'
                            })
                except json.JSONDecodeError:
                    print("DGIdb returned invalid JSON")
        except Exception as e:
            print(f"DGIdb error: {str(e)}")
        
        # Fallback/Additional data from Open Targets
        ot_url = "https://api.platform.opentargets.org/api/v4/graphql"
        headers = {'Content-Type': 'application/json'}
        for gene in genes[:5]:
            try:
                # Get Ensembl ID
                search_query = """
                query searchTarget($q: String!) {
                    search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 1}) {
                        hits { id }
                    }
                }
                """
                search_res = requests.post(ot_url, json={'query': search_query, 'variables': {'q': gene}}, timeout=10, headers=headers)
                if search_res.status_code != 200: 
                    print(f"OT Search failed for {gene}: {search_res.status_code}")
                    continue
                
                hits = search_res.json().get('data', {}).get('search', {}).get('hits', [])
                if not hits: continue
                target_id = hits[0]['id']

                # Get Known Drugs
                drug_query = """
                query getKnownDrugs($id: String!) {
                    target(ensemblId: $id) {
                        knownDrugs {
                            rows {
                                drug { name }
                                phase
                            }
                        }
                    }
                }
                """
                drug_res = requests.post(ot_url, json={'query': drug_query, 'variables': {'id': target_id}}, timeout=10, headers=headers)
                
                if drug_res.status_code == 200:
                    ot_data = drug_res.json().get('data', {}).get('target', {})
                    if not ot_data: continue
                    
                    known_drugs = ot_data.get('knownDrugs', {}).get('rows', [])
                    print(f"Found {len(known_drugs)} drugs for {gene} from OT")
                    for d in known_drugs[:10]:
                        drug_info = d.get('drug', {})
                        if drug_info.get('name'):
                            drugs.append({
                                'gene': gene,
                                'drug_name': drug_info['name'],
                                'interaction_types': [f"Phase {d.get('phase', '?')}"],
                                'source': 'Open Targets'
                            })
                else:
                    print(f"OT Drug query failed for {target_id}: {drug_res.status_code}")
                time.sleep(0.2)
            except Exception as e:
                print(f"Open Targets drug error for {gene}: {str(e)}")
        
        # Remove duplicates
        seen = set()
        unique_drugs = []
        for d in drugs:
            key = f"{d['gene']}_{d['drug_name']}".lower()
            if key not in seen:
                seen.add(key)
                unique_drugs.append(d)
        
        return Response({'drugs': unique_drugs})
