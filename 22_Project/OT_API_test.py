import requests

# reference: https://community.opentargets.org/t/downloading-variant-information/668/2

url = 'https://api.genetics.opentargets.org/graphql'

search_query  = """
	query searchRsId($rsId: String!) {
		search(queryString: $rsId) {
	    variants {
	      id
	    }
	  }
	}
"""

variant_query = """
	query variantInfo($variantId: String!) {
	  variantInfo(variantId: $variantId) {
		id
	    rsId
	    chromosome
	    position
	    refAllele
	    altAllele
	    nearestGeneDistance
	    nearestCodingGeneDistance
	  }
	}
"""

variables = {'rsId': 'rs10469840'} # Your rsID of interest

search_result = requests.post(url, json={'query': search_query, 'variables':variables}).json()
# This returns: {'data': {'search': {'variants': [{'id': '2_102476784_T_C'}]}}}

variables.update({'variantId': search_result['data']['search']['variants'][0]['id']})

variant_result = requests.post(url, json={'query': variant_query, 'variables':variables}).json()

print(variant_result)