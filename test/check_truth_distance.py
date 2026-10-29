#!/usr/bin/env python3
"""Distance between the reconstructed vertex and the true decay vertex.

usage: check_truth_distance.py <k4run output.root> <vertex vertex_collection>

The true decay vertex is the production vertex of the first final-state muon
in MCParticles (for HNL -> mu mu X both muons come from the same point).
Only the first (best) vertex of each event is used.
"""
import math
import sys
from podio.root_io import Reader

output_file = sys.argv[1]
vertex_collection = sys.argv[2]

number_of_events = 0
number_with_vertex = 0
distances = []

for event in Reader(output_file).get("events"):
    number_of_events += 1
    vertices = list(event.get(vertex_collection))
    if len(vertices) == 0:
        continue
    number_with_vertex += 1

    # true decay vertex = production point of a final-state muon
    true_vertex = None
    for particle in event.get("MCParticles"):
        if abs(particle.getPDG()) == 13 and particle.getGeneratorStatus() == 1:
            true_vertex = particle.getVertex()
            break
    if true_vertex is None:
        continue

    reconstructed_position = vertices[0].getPosition()
    distance = math.dist((reconstructed_position.x, reconstructed_position.y, reconstructed_position.z),
                         (true_vertex.x, true_vertex.y, true_vertex.z))
    distances.append(distance)

distances.sort()
print(f"{number_of_events} events, {number_with_vertex} with a vertex")
if distances:
    median = distances[len(distances) // 2]
    quantile_68 = distances[int(0.68 * len(distances))]
    print(f"|reconstructed - true vertex|: median {median:.3f} mm, 68% below {quantile_68:.3f} mm")
