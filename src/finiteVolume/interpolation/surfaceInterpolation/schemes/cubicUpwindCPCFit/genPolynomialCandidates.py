#!/usr/bin/env python3
import itertools
import sys

def each_degree_le_max(e):
    for i, max_e in enumerate(max_degrees):
        if e[i] > max_e:
            return False

    return True

def dense(s):
    s = set(s)

    for upper in s:
        lower_order_terms = set(itertools.product(*ranges(upper)))
        if not lower_order_terms.issubset(s):
            return False

    return True

def ranges(upper):
    return [tuple(range(x+1)) for x in upper]

def between_bounds(e, a, b):
    for i, x in enumerate(e):
        if x <= a[i] or x >= b[i]:
            return False

    return True

def format_polynomial(polynomial):
    return "candidates.append(" + " | ".join([format_term(t) for t in sorted(polynomial, key=lambda p: sum(p))]) + ");"

def format_term(term):
    s = ""
    algebra = ["x", "y", "z"]
    for i, t in enumerate(term):
        if t == 1:
            s += algebra[i]
        elif t > 1:
            s += algebra[i] * t

    if s == "":
        s = "constant"

    return s

if len(sys.argv) <= 1:
    print("""Usage: ./genPolynomialCandidates.py <max_degree,...>
  Generate OpenFOAM C++ source code of polynomial candidates.
  max_degree is the maximum polynomial degree in a particular dimension.

  Example:
    ./genPolynomialCandidates.py 3 2
  Generates the polynomial candidates for two-dimensional cubicFit""", file=sys.stderr)
    sys.exit(2)

if len(sys.argv) > 3+1:
    print("More than three dimensions are not supported", file=sys.stderr)
    sys.exit(2)

max_degrees = tuple([int(arg) for arg in sys.argv[1:]])

exponents = itertools.product(range(max(max_degrees)+1), repeat=len(max_degrees))
exponents = [e for e in exponents if sum(e) <= max(max_degrees) and each_degree_le_max(e)]

M = set(exponents)

for length in range(2,len(M)+1):
    subsets = itertools.combinations(M, length)
    dense_subsets = [s for s in subsets if dense(s)]
    print("\n".join([format_polynomial(s) for s in dense_subsets]))

# TODO: stricter constraints for including xy/xz/yz?
# TODO: remove yz^2 and y^2z?
