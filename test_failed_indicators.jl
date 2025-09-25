using Pioneer

# Run only up to MaxLFQ to check failed indicators
params_json = "./data/ecoli_test/ecoli_test_params.json"

# Run search - this will show our debug output
Pioneer.SearchDIA(params_json)