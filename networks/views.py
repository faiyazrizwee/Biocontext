"""
Networks App Views
Generate network data structures for frontend visualization.
"""
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status


class GenePathwayDiseaseNetworkView(APIView):
    """
    Generate Gene -> Pathway -> Disease network.
    POST: {
        "genes": [...],
        "pathways": [...],
        "diseases": [...]
    }
    """
    def post(self, request):
        genes = request.data.get('genes', [])
        pathways = request.data.get('pathways', [])
        diseases = request.data.get('diseases', [])
        
        nodes = []
        edges = []
        
        # Add gene nodes
        for gene in genes:
            nodes.append({
                'data': {
                    'id': f"gene_{gene}",
                    'label': gene,
                    'type': 'gene'
                }
            })
        
        # Add pathway nodes and edges
        for pathway in pathways:
            pw_id = pathway.get('pathway', pathway.get('id', ''))
            pw_genes = pathway.get('genes', [])
            
            nodes.append({
                'data': {
                    'id': f"pathway_{pw_id}",
                    'label': pw_id,
                    'type': 'pathway'
                }
            })
            
            for gene in pw_genes:
                if gene in genes:
                    edges.append({
                        'data': {
                            'id': f"e_gene_{gene}_pathway_{pw_id}",
                            'source': f"gene_{gene}",
                            'target': f"pathway_{pw_id}",
                            'type': 'gene_pathway'
                        }
                    })
        
        # Add disease nodes and edges
        for disease in diseases:
            d_id = disease.get('disease_id', disease.get('disease', ''))
            d_name = disease.get('disease', d_id)
            d_gene = disease.get('gene', '')
            
            node_id = f"disease_{d_id}"
            if not any(n['data']['id'] == node_id for n in nodes):
                nodes.append({
                    'data': {
                        'id': node_id,
                        'label': d_name,
                        'type': 'disease'
                    }
                })
            
            if d_gene and d_gene in genes:
                edges.append({
                    'data': {
                        'id': f"e_gene_{d_gene}_disease_{d_id}",
                        'source': f"gene_{d_gene}",
                        'target': node_id,
                        'type': 'gene_disease'
                    }
                })
        
        return Response({
            'elements': {
                'nodes': nodes,
                'edges': edges
            }
        })


class GeneDrugNetworkView(APIView):
    """
    Generate Gene -> Drug network.
    POST: {
        "genes": [...],
        "drugs": [...]
    }
    """
    def post(self, request):
        genes = request.data.get('genes', [])
        drugs = request.data.get('drugs', [])
        
        nodes = []
        edges = []
        
        # Add gene nodes
        for gene in genes:
            nodes.append({
                'data': {
                    'id': f"gene_{gene}",
                    'label': gene,
                    'type': 'gene'
                }
            })
        
        # Add drug nodes and edges
        added_drugs = set()
        for drug in drugs:
            drug_name = drug.get('drug_name', '')
            gene = drug.get('gene', '')
            interaction = drug.get('interaction_types', [])
            
            if drug_name and drug_name not in added_drugs:
                nodes.append({
                    'data': {
                        'id': f"drug_{drug_name}",
                        'label': drug_name,
                        'type': 'drug'
                    }
                })
                added_drugs.add(drug_name)
            
            if gene and gene in genes and drug_name:
                edges.append({
                    'data': {
                        'id': f"e_gene_{gene}_drug_{drug_name}",
                        'source': f"gene_{gene}",
                        'target': f"drug_{drug_name}",
                        'type': 'gene_drug',
                        'interaction': ', '.join(interaction) if interaction else 'unknown'
                    }
                })
        
        return Response({
            'elements': {
                'nodes': nodes,
                'edges': edges
            }
        })


class FullTherapyNetworkView(APIView):
    """
    Generate comprehensive Gene -> Disease -> Drug network.
    POST: {
        "genes": [...],
        "diseases": [...],
        "drugs": [...]
    }
    """
    def post(self, request):
        genes = request.data.get('genes', [])
        diseases = request.data.get('diseases', [])
        drugs = request.data.get('drugs', [])
        
        nodes = []
        edges = []
        added_ids = set()
        
        # Gene nodes
        for gene in genes:
            node_id = f"gene_{gene}"
            if node_id not in added_ids:
                nodes.append({'data': {'id': node_id, 'label': gene, 'type': 'gene'}})
                added_ids.add(node_id)
        
        # Disease nodes
        for d in diseases:
            d_id = d.get('disease_id', d.get('disease', ''))
            d_name = d.get('disease', d_id)
            d_gene = d.get('gene', '')
            node_id = f"disease_{d_id}"
            
            if node_id not in added_ids:
                nodes.append({'data': {'id': node_id, 'label': d_name, 'type': 'disease'}})
                added_ids.add(node_id)
            
            if d_gene and f"gene_{d_gene}" in added_ids:
                edges.append({
                    'data': {
                        'id': f"e_gene_{d_gene}_disease_{d_id}",
                        'source': f"gene_{d_gene}",
                        'target': node_id
                    }
                })
        
        # Drug nodes
        for drug in drugs:
            drug_name = drug.get('drug_name', '')
            gene = drug.get('gene', '')
            node_id = f"drug_{drug_name}"
            
            if drug_name and node_id not in added_ids:
                nodes.append({'data': {'id': node_id, 'label': drug_name, 'type': 'drug'}})
                added_ids.add(node_id)
            
            if gene and f"gene_{gene}" in added_ids and drug_name:
                edges.append({
                    'data': {
                        'id': f"e_gene_{gene}_drug_{drug_name}",
                        'source': f"gene_{gene}",
                        'target': node_id
                    }
                })
        
        return Response({
            'elements': {
                'nodes': nodes,
                'edges': edges
            }
        })
