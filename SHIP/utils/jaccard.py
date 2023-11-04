'''
Implements Jaccard distance and similarity.
'''
from typing import Iterable
import warnings
from utils.base_cluster import check_type_set

def jaccard(
    a: Iterable,
    b: Iterable,
    warn: bool = True
) -> float:
    '''
    Computes the Jaccard similarity between two vectors.
    Raises a warning for empty sets if warn is True. Always
    returns zero for empty vectors.
    '''

    a, b = check_type_set(a, b)

    if len(a)+len(b) == 0:
        if warn:
            warnings.warn(
                'Empty imput vectors. Score will be returned as 0.'
            )
        else:
            return 0

    return len(a.intersection(b))/len(a.union(b))

def jaccard_distance(
    a: Iterable,
    b: Iterable,
    warn: bool = False
) -> float:
    '''
    Calculates the complementary of the Jaccard index. If two sets are equal,
    the distance is 0. If they are completely distinct, the distance is 1.
    Throws a warning for empty sets when \'warn\' is True.
    '''

    return 1 - jaccard(a, b, warn=warn)

