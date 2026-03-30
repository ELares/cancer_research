#!/usr/bin/env python3
"""Quick parameter tuning: find params where all 4 validations pass."""
import subprocess, json, sys

# The key parameters to tune
# Problem: lp_propagation=0.012 is too high, causing baseline deaths
# Solution: reduce propagation AND increase repair to create tighter balance

configs = [
    # (lp_prop, fenton, gpx4_rate, fsp1_rate, lp_rate, death_thresh)
    (0.003, 0.05, 0.35, 0.08, 0.06, 10.0),
    (0.005, 0.04, 0.30, 0.07, 0.05, 10.0),
    (0.004, 0.04, 0.30, 0.07, 0.06, 12.0),
    (0.003, 0.04, 0.25, 0.06, 0.05, 10.0),
    (0.004, 0.05, 0.35, 0.08, 0.05, 12.0),
]

for i, (lp_prop, fenton, gpx4, fsp1, lp, thresh) in enumerate(configs):
    print(f"\nConfig {i}: prop={lp_prop} fen={fenton} gpx4={gpx4} fsp1={fsp1} lp={lp} thresh={thresh}")
    # We'd need to modify the rust code for each config — skip for now
    # Instead, print the analysis

print("\nThe core tension:")
print("- Propagation must be high enough that SDT pushes persisters over the edge")
print("- But low enough that OXPHOS basal ROS doesn't self-propagate")
print("- FSP1 loss must matter — it's the key differentiator")
print("- RSL3 must kill persisters (GPX4 inhibition alone)")
print("\nThe fix: make propagation dependent on GSH depletion")
print("When GSH is high, propagation is quenched. When GSH is depleted, propagation runs away.")
print("This creates the switch: SDT depletes GSH → unlocks propagation → bistable death")
