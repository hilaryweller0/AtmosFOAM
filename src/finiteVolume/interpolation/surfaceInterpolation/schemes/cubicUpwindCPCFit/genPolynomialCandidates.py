#!/usr/bin/env python3
import itertools

def each_degree_le_max(e):
    for i, max_e in enumerate(max_degrees):
        if e[i] > max_e:
            return False

    return True

def format(term):
    s = ""
    algebra = ["x", "y", "z"]
    for i, t in enumerate(term):
        if t == 1:
            s += algebra[i]
        elif t > 1:
            s += algebra[i] + "^" + str(t)
    return s

def ordered(pair):
    pair = tuple(pair)
    a = pair[0]
    b = pair[1]

    for i, x in enumerate(a):
        if x > b[i]:
            return False

    return True

def between_bounds(e, a, b):
    for i, x in enumerate(e):
        if x <= a[i] or x >= b[i]:
            return False

    return True

max_degrees = (3, 2) # x^3, y^2

exponents = itertools.product(range(max(max_degrees)+1), repeat=len(max_degrees))
exponents = [e for e in exponents if sum(e) <= max(max_degrees) and each_degree_le_max(e)]

M = set(exponents)

pairs = itertools.combinations(M, 2)
dense_subsets = []

for pair in pairs:
    if ordered(pair):
        a, b = pair
    elif ordered(reversed(pair)):
        b, a = pair
    else:
        continue

    print("bounds", a, b)
    S = {e for e in M if (e in pair) or between_bounds(e, a, b)}
    print(S)
    print()
    # now, filter M using lower bound 'a' and upper bound 'b'

# TODO: remove yz^2 and y^2z ?
