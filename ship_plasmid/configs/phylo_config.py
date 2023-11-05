PHYLO_CONFIG = {
  'agglomerative-clustering': {
    'linkage': 'average',
    'distance-threshold': .4,
    'adaptative': True
  },
  'learn-weights': True,
  'learning-method': 'max-likelihood', # ratio, precomputed or max-likelihood
  'weights': {
    'mutation': 1,
    'mobile element': 1,
    'recombination': 20,
    'reversal': 1,
    'transposition': 1,
    'duplicate': 1
  },
  'break-distance': {
    'k': 3,
    'indel-size-penalty': 0.5
  }    
}
