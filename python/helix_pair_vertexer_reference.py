"""
Python reference of the helix pair vertexer (HelixPairVertexFinder) for muon-system tracks.
It works the same way in the barrel and in the endcaps.

Each track
  - A helix is fitted to the track's own hits, in the field where the muon chambers are
    (Bz = -1.7 T everywhere outside the coil in IDEA_o1_v05, also in the forward volume
    |z| > 2.62 m):
      * in xy: a circle (algebraic fit, then a geometric refinement)
      * in z:  a straight line in the turning angle
  - The track is described by its 3D path length from the innermost hit (negative = towards
    the decay). So a track crossing endcap discs (small xy motion, large z motion) is
    described as well as a barrel track.
  - Tracks without a measurable curvature (xy span below MINIMUM_XY_SPAN, or radius above
    MAXIMUM_RADIUS) are described by a straight 3D line through the hits.

Each pair of tracks
  - The vertex is the point of closest approach of the two trajectories in 3D (midpoint of
    the two closest points); the distance of closest approach (DCA) is their distance.
  - No region is assumed: the caller decides barrel / endcap / outside.

Momentum
  - pT = 0.3 |B| R with |B| in the chambers, p = pT * sqrt(1 + tanLambda^2), plus the
    ionisation loss between the vertex and the first hit (energy loss per mm in the
    dual-readout calorimeter and in the iron, from the material map below), calibrated on
    truth.

Units: mm, GeV, Tesla.
"""
import numpy as np
from scipy.optimize import minimize

# ---- detector: magnetic field of the simulation (IDEA_o1_v05 compact) ---------------------
FIELD_IN_CHAMBERS = 1.7          # T, |Bz| in the muon chambers, calorimeter and forward volume
FIELD_INSIDE_COIL = 2.0          # T, Bz inside the coil volume
COIL_RADIUS = 2300.0             # mm
FIELD_HALF_LENGTH = 2616.0       # mm, 1.2 x the coil half length

# ---- detector: material map (IDEA_o1_v05) ---------------------------------------------------
# dual-readout calorimeter = shell between the box (r < 2475, |z| < 2475) and the box
# (r < 4525, |z| < 4525), without the beam hole r < 252; iron everywhere outside it
CALORIMETER_INNER_BOX = 2475.0   # mm
CALORIMETER_OUTER_BOX = 4525.0   # mm
CALORIMETER_BEAM_HOLE = 252.0    # mm
# energy loss in GeV/mm, calibrated on truth: (p_true - p_measured) = 1.36 x the nominal loss
# (1.0 in the calorimeter, 1.45 in the iron); with it p_reco/p_true = 1.00 (median)
ENERGY_LOSS_PER_MM = {"calorimeter": 1.36e-3, "iron": 1.97e-3, "none": 0.0}

# ---- track model --------------------------------------------------------------------------
MINIMUM_XY_SPAN = 30.0           # mm, minimum xy extent of the hits for a circle fit
MAXIMUM_RADIUS = 2.0e6           # mm, a larger radius is treated as a straight line
MUON_MASS = 0.1056584            # GeV


def field_z_at(position):
    """Bz of the simulation field map at a point (x, y, z) in mm."""
    radius_xy = np.hypot(position[0], position[1])
    inside_coil_volume = radius_xy < COIL_RADIUS and abs(position[2]) < FIELD_HALF_LENGTH
    if inside_coil_volume:
        return FIELD_INSIDE_COIL
    return -FIELD_IN_CHAMBERS


def material_region(points):
    """For an array of points (..., 3): two boolean arrays, inside the calorimeter / inside the iron."""
    radius_xy = np.hypot(points[..., 0], points[..., 1])
    absolute_z = np.abs(points[..., 2])
    inside_outer_box = (radius_xy < CALORIMETER_OUTER_BOX) & (absolute_z < CALORIMETER_OUTER_BOX)
    inside_inner_box = (radius_xy < CALORIMETER_INNER_BOX) & (absolute_z < CALORIMETER_INNER_BOX)
    in_calorimeter = inside_outer_box & ~inside_inner_box & (radius_xy > CALORIMETER_BEAM_HOLE)
    in_iron = ~inside_outer_box
    return in_calorimeter, in_iron


def energy_loss_per_mm_at(position):
    """Energy loss [GeV/mm] at one point."""
    in_calorimeter, in_iron = material_region(np.asarray(position)[None, :])
    return (ENERGY_LOSS_PER_MM["calorimeter"] * in_calorimeter[0]
            + ENERGY_LOSS_PER_MM["iron"] * in_iron[0])


def energy_loss_on_segment(start, end, number_of_steps=60):
    """Ionisation loss [GeV] on the straight segment start -> end."""
    fractions = (np.arange(number_of_steps) + 0.5) / number_of_steps
    points = start[None, :] + fractions[:, None] * (end - start)[None, :]
    in_calorimeter, in_iron = material_region(points)
    step_length = np.linalg.norm(end - start) / number_of_steps
    return float(step_length * (ENERGY_LOSS_PER_MM["calorimeter"] * in_calorimeter.sum()
                                + ENERGY_LOSS_PER_MM["iron"] * in_iron.sum()))


class HelixTrack:
    """Helix (or straight line) fitted to the hits of one track."""

    def __init__(self, hits):
        self.is_valid = False
        hits = np.asarray(hits, float)
        # order the hits from the IP outward
        distance_from_origin = np.linalg.norm(hits, axis=1)
        hits = hits[np.argsort(distance_from_origin)]
        self.hits = hits
        self.number_of_hits = len(hits)
        if self.number_of_hits < 3:
            return
        self.innermost_hit = hits[0]
        hit_x = hits[:, 0]
        hit_y = hits[:, 1]
        hit_z = hits[:, 2]
        xy_span = np.hypot(hit_x - hit_x[0], hit_y - hit_y[0]).max()
        self.is_straight_line = True

        if xy_span >= MINIMUM_XY_SPAN:
            # algebraic circle fit: x^2 + y^2 = 2 x xc + 2 y yc + (R^2 - xc^2 - yc^2)
            design_matrix = np.c_[2 * hit_x, 2 * hit_y, np.ones(self.number_of_hits)]
            solution, *_ = np.linalg.lstsq(design_matrix, hit_x ** 2 + hit_y ** 2, rcond=None)
            centre_x = solution[0]
            centre_y = solution[1]
            radius = np.sqrt(max(solution[2] + centre_x * centre_x + centre_y * centre_y, 1e-9))

            # geometric refinement: Gauss-Newton on the distances of the hits to the circle
            for _ in range(8):
                offset_x = hit_x - centre_x
                offset_y = hit_y - centre_y
                distance_to_centre = np.hypot(offset_x, offset_y)
                jacobian = np.c_[-offset_x / distance_to_centre, -offset_y / distance_to_centre,
                                 -np.ones(self.number_of_hits)]
                residuals = distance_to_centre - radius
                try:
                    correction = np.linalg.lstsq(jacobian, -residuals, rcond=None)[0]
                except np.linalg.LinAlgError:
                    break
                centre_x = centre_x + correction[0]
                centre_y = centre_y + correction[1]
                radius = radius + correction[2]
                if np.abs(correction).max() < 1e-3:
                    break

            if np.isfinite(radius) and 0 < radius < MAXIMUM_RADIUS:
                # angle of every hit around the circle centre, without 2 pi jumps
                angle_around_centre = np.unwrap(np.arctan2(hit_y - centre_y, hit_x - centre_x))
                # rotation sense of the outward motion: +1 counter-clockwise, -1 clockwise
                if angle_around_centre[-1] > angle_around_centre[0]:
                    self.rotation_sense = 1.0
                else:
                    self.rotation_sense = -1.0
                # turning angle along the motion, 0 at the innermost hit
                turning_angle = self.rotation_sense * (angle_around_centre - angle_around_centre[0])
                if np.ptp(turning_angle) > 0:
                    z_per_radian = np.polyfit(turning_angle, hit_z, 1)[0]
                else:
                    z_per_radian = 0.0
                self.centre_x = centre_x
                self.centre_y = centre_y
                self.radius = radius
                self.z_per_radian = z_per_radian
                self.angle_at_innermost_hit = angle_around_centre[0]
                self.z_at_innermost_hit = float(np.mean(hit_z - z_per_radian * turning_angle))
                # 3D path length per radian of turning
                self.path_length_per_radian = np.sqrt(radius * radius + z_per_radian * z_per_radian)
                self.is_straight_line = False
                self.circle_chi2_xy = float(np.sum((np.hypot(hit_x - centre_x, hit_y - centre_y) - radius) ** 2))

        if self.is_straight_line:
            # principal axis of the hits, oriented outward
            hits_centre = hits.mean(axis=0)
            principal_axis = np.linalg.svd(hits - hits_centre)[2][0]
            if np.dot(hits[-1] - hits[0], principal_axis) < 0:
                principal_axis = -principal_axis
            self.line_centre = hits_centre
            self.line_direction = principal_axis / np.linalg.norm(principal_axis)
            self.radius = np.inf
            self.z_per_radian = 0.0
            self.rotation_sense = 0.0
        self.is_valid = True

        # measured momentum at the first hit
        if self.is_straight_line:
            self.measured_momentum = np.nan
        else:
            transverse_momentum = 0.3 * FIELD_IN_CHAMBERS * self.radius * 1e-3
            tan_lambda = self.z_per_radian / self.radius
            self.measured_momentum = transverse_momentum * np.sqrt(1.0 + tan_lambda ** 2)
        self.path_length_to_last_hit = float(np.linalg.norm(self.hits[-1] - self.innermost_hit))

    def point(self, path_length):
        """Point(s) of the trajectory at 3D path length(s) from the innermost hit."""
        path_length = np.asarray(path_length, float)
        if self.is_straight_line:
            # projection of the innermost hit on the line
            line_start = (self.line_centre
                          + np.dot(self.innermost_hit - self.line_centre, self.line_direction) * self.line_direction)
            return line_start + path_length[..., None] * self.line_direction
        turning_angle = path_length / self.path_length_per_radian
        angle = self.angle_at_innermost_hit + self.rotation_sense * turning_angle
        return np.stack([self.centre_x + self.radius * np.cos(angle),
                         self.centre_y + self.radius * np.sin(angle),
                         self.z_at_innermost_hit + self.z_per_radian * turning_angle], axis=-1)

    def tangent(self, path_length):
        """Unit vector of the outward motion at a path length."""
        if self.is_straight_line:
            return self.line_direction.copy()
        angle = self.angle_at_innermost_hit + self.rotation_sense * path_length / self.path_length_per_radian
        direction = np.array([-self.radius * np.sin(angle) * self.rotation_sense,
                              self.radius * np.cos(angle) * self.rotation_sense,
                              self.z_per_radian])
        return direction / np.linalg.norm(direction)


def propagated_trajectory(track, step_length=20.0, backward_factor=1.3, upstream_tolerance=100.0):
    """
    Polyline of the track, stepped INWARD from the innermost hit through the simulation field
    map while the ionisation loss of the material map is added back (the curvature follows the
    local field and momentum), plus the fitted helix up to upstream_tolerance mm outward.
    Returns: path lengths (<= upstream_tolerance), points, momenta, outward unit directions.
    """
    # outward part: the fitted helix
    outward_path_lengths = np.arange(upstream_tolerance, 0.0, -step_length)[::-1]
    if len(outward_path_lengths):
        outward_points = track.point(outward_path_lengths)
    else:
        outward_points = np.zeros((0, 3))
    outward_directions = np.array([track.tangent(length) for length in outward_path_lengths]).reshape(-1, 3)

    # inward part: stepping
    number_of_steps = int(backward_factor * np.linalg.norm(track.innermost_hit) / step_length)
    position = track.point(0.0).copy()
    direction = track.tangent(0.0)
    cos_lambda = np.hypot(direction[0], direction[1])
    tan_lambda = direction[2] / max(cos_lambda, 1e-12)
    azimuth = np.arctan2(direction[1], direction[0])
    if track.is_straight_line:
        momentum = np.nan
    else:
        momentum = track.measured_momentum
    points = [position.copy()]
    directions = [direction.copy()]
    momenta = [momentum]
    for _ in range(number_of_steps):
        if not track.is_straight_line:
            transverse_momentum = momentum * cos_lambda
            # signed xy curvature of the outward motion [1/mm]
            curvature = -track.rotation_sense * field_z_at(position) * 0.3e-3 / max(transverse_momentum, 1e-3)
            position = position - step_length * np.array([cos_lambda * np.cos(azimuth),
                                                          cos_lambda * np.sin(azimuth),
                                                          cos_lambda * tan_lambda])
            azimuth = azimuth - curvature * step_length * cos_lambda
            momentum = momentum + energy_loss_per_mm_at(position) * step_length
            direction = np.array([cos_lambda * np.cos(azimuth), cos_lambda * np.sin(azimuth), cos_lambda * tan_lambda])
        else:
            position = position - step_length * direction
        points.append(position.copy())
        directions.append(direction.copy())
        momenta.append(momentum)

    # order from the innermost point outward, then append the outward part
    points = np.array(points[::-1])
    directions = np.array(directions[::-1])
    momenta = np.array(momenta[::-1])
    path_lengths = np.r_[-step_length * np.arange(number_of_steps, -1, -1), outward_path_lengths]
    all_points = np.r_[points, outward_points]
    all_momenta = np.r_[momenta, np.full(len(outward_path_lengths), track.measured_momentum)]
    all_directions = np.r_[directions, outward_directions]
    return path_lengths, all_points, all_momenta, all_directions


def segment_closest_approach(first_start, first_end, second_start, second_end):
    """Closest approach of two segments: fraction along each segment and the distance."""
    first_vector = first_end - first_start
    second_vector = second_end - second_start
    between_starts = first_start - second_start
    first_length2 = first_vector @ first_vector
    dot_products = first_vector @ second_vector
    second_length2 = second_vector @ second_vector
    first_dot_offset = first_vector @ between_starts
    second_dot_offset = second_vector @ between_starts
    denominator = first_length2 * second_length2 - dot_products * dot_products
    if denominator > 1e-12 * first_length2 * second_length2:
        second_fraction = np.clip((first_length2 * second_dot_offset - dot_products * first_dot_offset) / denominator,
                                  0.0, 1.0)
    elif second_length2 > 0:
        second_fraction = np.clip(second_dot_offset / second_length2, 0.0, 1.0)
    else:
        second_fraction = 0.0
    if first_length2 > 0:
        first_fraction = np.clip((dot_products * second_fraction - first_dot_offset) / first_length2, 0.0, 1.0)
    else:
        first_fraction = 0.0
    if second_length2 > 0:
        second_fraction = np.clip((dot_products * first_fraction + second_dot_offset) / second_length2, 0.0, 1.0)
    else:
        second_fraction = 0.0
    distance = float(np.linalg.norm(first_start + first_fraction * first_vector
                                    - (second_start + second_fraction * second_vector)))
    return first_fraction, second_fraction, distance


def pair_kinematics(result, first_momentum_vector, second_momentum_vector):
    """Add the pair mass, the opening angle and the two momenta to a vertex result."""
    energy = sum(np.sqrt(vector @ vector + MUON_MASS ** 2)
                 for vector in (first_momentum_vector, second_momentum_vector))
    total_momentum = first_momentum_vector + second_momentum_vector
    result["mass"] = float(np.sqrt(max(energy * energy - total_momentum @ total_momentum, 0.0)))
    result["cos_opening"] = float(first_momentum_vector @ second_momentum_vector
                                  / (np.linalg.norm(first_momentum_vector) * np.linalg.norm(second_momentum_vector)))


def pair_vertex_stepped(first_track, second_track, upstream_tolerance=100.0, step_length=20.0):
    """
    Vertex of two tracks from their stepped trajectories (the method of the C++ finder):
    3D closest approach, vertex = midpoint; momenta (energy loss added back to the vertex)
    and directions at the vertex.
    """
    if not (first_track.is_valid and second_track.is_valid):
        return None
    (first_path, first_points, first_momenta, first_directions) = propagated_trajectory(
        first_track, step_length=step_length, upstream_tolerance=upstream_tolerance)
    (second_path, second_points, second_momenta, second_directions) = propagated_trajectory(
        second_track, step_length=step_length, upstream_tolerance=upstream_tolerance)

    # coarse search: the closest pair of trajectory points
    distances = np.linalg.norm(first_points[:, None, :] - second_points[None, :, :], axis=2)
    first_index, second_index = np.unravel_index(np.argmin(distances), distances.shape)

    # fine search: closest approach of the neighbouring segments
    best = (np.inf, 0, 0, 0.0, 0.0)
    for first_segment in (first_index - 1, first_index):
        for second_segment in (second_index - 1, second_index):
            first_ok = 0 <= first_segment < len(first_points) - 1
            second_ok = 0 <= second_segment < len(second_points) - 1
            if first_ok and second_ok:
                first_fraction, second_fraction, distance = segment_closest_approach(
                    first_points[first_segment], first_points[first_segment + 1],
                    second_points[second_segment], second_points[second_segment + 1])
                if distance < best[0]:
                    best = (distance, first_segment, second_segment, first_fraction, second_fraction)
    distance, first_segment, second_segment, first_fraction, second_fraction = best
    if not np.isfinite(distance):
        distance = float(distances[first_index, second_index])
        first_segment = min(first_index, len(first_points) - 2)
        second_segment = min(second_index, len(second_points) - 2)
        first_fraction = 0.0
        second_fraction = 0.0

    first_closest = first_points[first_segment] + first_fraction * (
        first_points[first_segment + 1] - first_points[first_segment])
    second_closest = second_points[second_segment] + second_fraction * (
        second_points[second_segment + 1] - second_points[second_segment])
    first_vertex_path = first_path[first_segment] + first_fraction * (
        first_path[first_segment + 1] - first_path[first_segment])
    second_vertex_path = second_path[second_segment] + second_fraction * (
        second_path[second_segment + 1] - second_path[second_segment])
    vertex = 0.5 * (first_closest + second_closest)

    # decay between the IP and the chambers: closer to the origin than every first hit, same side
    between_ip_and_hits = all(np.linalg.norm(vertex) < np.linalg.norm(track.innermost_hit)
                              and np.dot(vertex, track.innermost_hit) > 0
                              for track in (first_track, second_track))
    result = dict(vertex=vertex, distance_of_closest_approach=float(distance),
                  path_lengths=(float(first_vertex_path), float(second_vertex_path)),
                  upstream=bool(first_vertex_path < upstream_tolerance and second_vertex_path < upstream_tolerance),
                  between_ip_and_hits=bool(between_ip_and_hits),
                  radius_xy=float(np.hypot(vertex[0], vertex[1])), z=float(vertex[2]))

    first_direction = first_directions[first_segment]
    second_direction = second_directions[second_segment]
    if first_track.is_straight_line or second_track.is_straight_line:
        result["mass"] = np.nan
        result["cos_opening"] = float(first_direction @ second_direction)
        result["momenta"] = (np.nan, np.nan)
        return result
    first_momentum_vector = first_momenta[first_segment] * first_direction
    second_momentum_vector = second_momenta[second_segment] * second_direction
    pair_kinematics(result, first_momentum_vector, second_momentum_vector)
    result["momenta"] = (float(first_momenta[first_segment]), float(second_momenta[second_segment]))
    return result


def pair_vertex(first_track, second_track, upstream_tolerance=100.0, backward_factor=1.3, grid_points=90):
    """
    Vertex of two tracks from their fitted helices extrapolated inward (no field map, no
    stepping): grid search in the two path lengths, then a local minimisation of the distance.
    """
    if not (first_track.is_valid and second_track.is_valid):
        return None
    first_lowest = -backward_factor * np.linalg.norm(first_track.innermost_hit)
    second_lowest = -backward_factor * np.linalg.norm(second_track.innermost_hit)
    first_grid = np.linspace(first_lowest, upstream_tolerance, grid_points)
    second_grid = np.linspace(second_lowest, upstream_tolerance, grid_points)
    first_points = first_track.point(first_grid)
    second_points = second_track.point(second_grid)
    distances = np.linalg.norm(first_points[:, None, :] - second_points[None, :, :], axis=2)
    first_index, second_index = np.unravel_index(np.argmin(distances), distances.shape)

    def distance_between_tracks(path_lengths):
        return float(np.linalg.norm(first_track.point(path_lengths[0]) - second_track.point(path_lengths[1])))

    first_grid_step = first_grid[1] - first_grid[0]
    second_grid_step = second_grid[1] - second_grid[0]
    start = [first_grid[first_index], second_grid[second_index]]
    initial_simplex = [start,
                       [start[0] + first_grid_step, start[1]],
                       [start[0], start[1] + second_grid_step]]
    minimum = minimize(distance_between_tracks, x0=start, method="Nelder-Mead",
                       options=dict(xatol=0.5, fatol=0.05, maxiter=400, initial_simplex=initial_simplex))
    first_vertex_path, second_vertex_path = minimum.x
    first_closest = first_track.point(first_vertex_path)
    second_closest = second_track.point(second_vertex_path)
    vertex = 0.5 * (first_closest + second_closest)
    upstream = first_vertex_path < upstream_tolerance and second_vertex_path < upstream_tolerance
    between_ip_and_hits = all(np.linalg.norm(vertex) < np.linalg.norm(track.innermost_hit)
                              and np.dot(vertex, track.innermost_hit) > 0
                              for track in (first_track, second_track))
    result = dict(vertex=vertex, distance_of_closest_approach=float(np.linalg.norm(first_closest - second_closest)),
                  path_lengths=(float(first_vertex_path), float(second_vertex_path)),
                  upstream=bool(upstream), between_ip_and_hits=bool(between_ip_and_hits),
                  radius_xy=float(np.hypot(vertex[0], vertex[1])), z=float(vertex[2]))

    # kinematics at the vertex
    if first_track.is_straight_line or second_track.is_straight_line:
        result["mass"] = np.nan
        result["cos_opening"] = float(np.dot(first_track.tangent(first_vertex_path),
                                             second_track.tangent(second_vertex_path)))
        result["momenta"] = (first_track.measured_momentum, second_track.measured_momentum)
        return result
    momentum_vectors = []
    for track, vertex_path in ((first_track, first_vertex_path), (second_track, second_vertex_path)):
        momentum_at_vertex = track.measured_momentum + energy_loss_on_segment(vertex, track.innermost_hit)
        momentum_vectors.append(momentum_at_vertex * track.tangent(vertex_path))
    pair_kinematics(result, momentum_vectors[0], momentum_vectors[1])
    result["momenta"] = (float(np.linalg.norm(momentum_vectors[0])), float(np.linalg.norm(momentum_vectors[1])))
    return result
