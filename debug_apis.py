import requests
import json

def test_opentargets(gene):
    url = "https://api.platform.opentargets.org/api/v4/graphql"
    # Step 1: Search for Ensembl ID
    search_query = """
    query searchTarget($q: String!) {
        search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 1}) {
            hits {
                id
                name
            }
        }
    }
    """
    response = requests.post(url, json={'query': search_query, 'variables': {'q': gene}}, timeout=20)
    if response.status_code != 200:
        print(f"Search failed: {response.text}")
        return

    search_data = response.json()
    hits = search_data.get('data', {}).get('search', {}).get('hits', [])
    if not hits:
        print(f"No hits for {gene}")
        return

    target_id = hits[0]['id']
    print(f"Found ID for {gene}: {target_id}")

    # Step 2: Get associations
    assoc_query = """
    query targetAssoc($id: String!) {
        target(ensemblId: $id) {
            associatedDiseases(page: {index: 0, size: 10}) {
                rows {
                    disease {
                        id
                        name
                    }
                    score
                }
            }
            knownDrugs(page: {index: 0, size: 10}) {
                rows {
                    drug {
                        name
                    }
                    phase
                }
            }
        }
    }
    """
    response = requests.post(url, json={'query': assoc_query, 'variables': {'id': target_id}}, timeout=20)
    print(f"Open Targets Assoc ({gene}) status: {response.status_code}")
    if response.status_code == 200:
        data = response.json()
        print(json.dumps(data, indent=2))
    else:
        print(response.text)

def test_dgidb(genes):
    url = "https://www.dgidb.org/api/v2/interactions.json"
    params = {'genes': ','.join(genes)}
    response = requests.get(url, params=params, timeout=30)
    print(f"DGIdb status: {response.status_code}")
    if response.status_code == 200:
        data = response.json()
        # Just print summary
        print(f"Matched terms: {len(data.get('matchedTerms', []))}")
        if data.get('matchedTerms'):
            print(f"First match interactions: {len(data['matchedTerms'][0].get('interactions', []))}")
    else:
        print(response.text)

if __name__ == "__main__":
    test_opentargets("TP53")
    # test_dgidb(["TP53", "BRCA1"])
