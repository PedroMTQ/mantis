#cython: language_level=3
from time import time

def add_possible_combo(dict possible_combos,set temp_combo):
    if len(temp_combo) not in possible_combos:
        possible_combos[len(temp_combo)]=set()
    if not already_added(possible_combos[len(temp_combo)],temp_combo):
        possible_combos[len(temp_combo)].add(frozenset(temp_combo))

def already_added(set all_combos,set new_combo):
    if not all_combos: return False
    if not new_combo: return True
    cdef frozenset old_combo
    cdef int c
    for old_combo in all_combos:
        if len(old_combo) == len(new_combo):
            c = len(new_combo)
            for k in old_combo:
                if k in new_combo: c -= 1
            if c == 0: return True
    return False

def is_overlap(set temp_queries, tuple current_query):
    if not temp_queries or not current_query: return False
    #'env_coord_from'=1
    #'env_coord_to'=2
    #'hmm_name'=3
    y = range(current_query[1], current_query[2] + 1)
    for t in temp_queries:
        if t[3]==current_query[3]: return True
        x = range(t[1], t[2] + 1)
        xs = set(x)
        res = xs.intersection(y)
        if res: return True
    return False

def check_children(set parent, set children):
    cdef set res
    cdef tuple c
    cdef set temp_children
    res = set()
    temp_children=set(children)
    for c in temp_children:
        if is_overlap(parent, c):
            children.remove(c)
        else: res.add(c)
    return res



def add_to_combo(set possible_children, set  combo, dict possible_combos, set query_hits, double start_time, double time_limit):
    cdef tuple child_hit
    cdef set temp_combo
    cdef set temp_query_hits
    cdef set temp_possible_children
    #when this algorithm takes too long, we just return an empty dictionary
    if time()-start_time>time_limit:
      raise TimeoutError
    if not possible_children:
      if combo:
        add_possible_combo(possible_combos,combo)
    for child_hit in possible_children:
        temp_combo = set(combo)
        temp_combo.add(child_hit)
        temp_query_hits=set(query_hits)
        #add_possible_combo(possible_combos,temp_combo)
        # in order to avoid adding the same combinations
        temp_possible_children = check_children(temp_combo, temp_query_hits)
        add_to_combo(temp_possible_children, temp_combo, possible_combos, temp_query_hits,start_time,time_limit)
    return possible_combos

def get_non_overlapping_hits(set query_hits, double time_limit):
    '''
    The immediate approach to this would be to just test all combinations, but this grows exponentially according to the number of hits
    2 hits: 2^2=4
    5 hits: 2^5=120
    15 hits : 2^15=32767
    So instead we only generate all POSSIBLE NON-REPEATED combinations, thus reducing search space
    '''
    cdef double start_time = time()
    cdef dict possible_combos = dict()
    cdef set combo = set()
    cdef set possible_children = set(query_hits)
    res = add_to_combo(possible_children, combo, possible_combos, query_hits, start_time, time_limit)
    return res
