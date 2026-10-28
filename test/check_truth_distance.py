#!/usr/bin/env python3
"""Distance between the reconstructed vertex and the true decay vertex.

usage: check_truth_distance.py <k4run output.root> <vertex collection>

The true decay vertex is the production vertex of the first final-state muon
in MCParticles (for HNL -> mu mu X both muons come from the same point).
Only the first (best) vertex of each event is used.
"""
import math
import sys
from podio.root_io import Reader

file_name = sys.argv[1]
collection = sys.argv[2]

n_events = 0
n_with_vertex = 0
distances = []

for event in Reader(file_name).get("events"):
    n_events += 1
    vertices = list(event.get(collection))
    if len(vertices) == 0:
        continue
    n_with_vertex += 1

    # true decay vertex = production point of a final-state muon
    true_vertex = None
    for particle in event.get("MCParticles"):
        if abs(particle.getPDG()) == 13 and particle.getGeneratorStatus() == 1:
            true_vertex = particle.getVertex()
            break
    if true_vertex is None:
        continue

    reco = vertices[0].getPosition()
    distance = math.dist((reco.x, reco.y, reco.z), (true_vertex.x, true_vertex.y, true_vertex.z))
    distances.append(distance)

distances.sort()
print(f"{n_events} events, {n_with_vertex} with a vertex")
if distances:
    median = distances[len(distances) // 2]
    p68 = distances[int(0.68 * len(distances))]
    print(f"|reco - true vertex|: median {median:.3f} mm, 68% below {p68:.3f} mm")
