#!/usr/bin/env python3
with open("build/thermalAdvectionSlantedCellCubicUpwindCPCFit/old.log") as f:
    old = f.readlines()

with open("build/thermalAdvectionSlantedCellCubicUpwindCPCFit/new.log") as f:
    new = f.readlines()

assert len(old) == len(new)

epsilon = 0.01

larger = 0
smaller = 0
same = 0

moreTerms = 0
fewerTerms = 0
sameTerms = 0

for old_stencil, new_stencil in zip(old, new):
    if old_stencil.startswith("==="):
        continue
    old_tokens = old_stencil.split(" ")
    old_upwind = float(old_tokens[4])

    new_tokens = new_stencil.split(" ")
    new_upwind = float(new_tokens[4])

    if new_upwind > old_upwind+epsilon:
        larger += 1
    elif old_upwind > new_upwind+epsilon:
        smaller += 1
        old_polyTerms = float(old_tokens[6])
        new_polyTerms = float(new_tokens[6])
        if new_polyTerms > old_polyTerms:
            moreTerms += 1
        elif new_polyTerms < old_polyTerms:
            fewerTerms += 1
        else:
            sameTerms += 1
    else:
        same += 1

print("smaller, same, larger")
print(smaller, same, larger)

print("fewerTerms, sameTerms, moreTerms")
print(fewerTerms, sameTerms, moreTerms)
