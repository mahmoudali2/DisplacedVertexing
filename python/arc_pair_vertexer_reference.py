#!/usr/bin/env python3
"""
Python reference of the arc pair vertexer (ArcPairVertexFinder). Units: mm, GeV, Tesla.

Each track
  - The hits are ordered from the IP outward. A local frame is attached to the innermost hit:
      outward_direction_xy = outward direction in xy (from a first straight-line fit),
      normal_direction_xy  = z x outward_direction_xy,
      arc position s        = coordinate along outward_direction_xy,
      offset u              = coordinate along normal_direction_xy.
  - xy: circle through the hits, fitted as u(s) = a + b s + curvature s^2 / 2 (generalised
    least squares with the full hit covariance), then used as the exact circle with that
    tangent and curvature.
  - z: straight line z(s) = c + d s (generalised least squares).
  - Hit covariance = hit resolution (diagonal) + correlated multiple scattering from the
    material between the first hit and each hit (muon-system yoke).

Error of the extrapolated track at an arc position s < 0 (before the first hit)
  = fit covariance propagated with (1, s, s^2/2) in xy and (1, s) in z
  + multiple scattering in the material between that point and the first hit (calorimeter)
  + an optional field-mismatch term.

Vertex fit of a track pair
  Gauss-Newton on two residuals per track: the xy distance perpendicular to the track and the
  z difference at the same s. chi2 with ndf = 2 x tracks - 3.
"""
import numpy as np

HIGHLAND_CONSTANT = 13.6e-3   # GeV, multiple-scattering formula


class Config:
    # ---- detector: magnetic field --------------------------------------------------------
    field_for_momentum = 2.0            # T, converts the track omega into pT
    # the circle is fitted in the chambers' field and extrapolated through the calorimeter
    # region, whose field differs; calibrated on truth (pulls flat vs extrapolation): 0.5 T
    field_uncertainty = 0.5             # T
    # ---- detector: resolutions ----------------------------------------------------------
    momentum_resolution = 0.20          # sigma(p)/p of the standalone muon system
    hit_resolution = 1.0                # mm, per coordinate
    # ---- detector: material -------------------------------------------------------------
    use_material = True                 # multiple scattering in yoke and calorimeter
    calorimeter_radiation_lengths = 90.0
    calorimeter_inner_radius, calorimeter_outer_radius = 2400.0, 4500.0
    calorimeter_inner_half_z, calorimeter_outer_half_z = 2400.0, 4300.0
    yoke_radiation_lengths = 40.0       # muon-system yoke between the stations
    yoke_thickness = 2000.0             # mm
    # ---- physics ------------------------------------------------------------------------
    minimum_momentum = 2.0              # GeV, floor for the scattering estimate
    # ---- algorithm ----------------------------------------------------------------------
    fit_curvature = True
    minimum_hits = 3
    upstream_tolerance = 100.0          # mm, vertex may sit this far past the first hit
    integration_steps = 40
    # ---- vertex region: only vertices behind the solenoid are kept ------------------------
    solenoid_radius, solenoid_half_z = 2300.0, 2180.0


def calorimeter_radiation_lengths_per_mm(config):
    return config.calorimeter_radiation_lengths / (config.calorimeter_outer_radius - config.calorimeter_inner_radius)


def radiation_lengths_per_mm(position, config):
    """Radiation lengths per mm of path at a point."""
    radius_xy = np.hypot(position[0], position[1])
    absolute_z = abs(position[2])
    inside_outer_box = radius_xy < config.calorimeter_outer_radius and absolute_z < config.calorimeter_outer_half_z
    inside_inner_box = radius_xy < config.calorimeter_inner_radius and absolute_z < config.calorimeter_inner_half_z
    if inside_outer_box and not inside_inner_box:
        return calorimeter_radiation_lengths_per_mm(config)
    if not inside_outer_box:
        return config.yoke_radiation_lengths / config.yoke_thickness
    return 0.0


def logarithm_correction(total_radiation_lengths):
    """(1 + 0.038 ln x)^2 term of the Highland formula."""
    return (1.0 + 0.038 * np.log(max(total_radiation_lengths, 1e-3))) ** 2


class ArcTrack:
    """Local arc fit of one track."""

    def __init__(self, hits, omega, config):
        self.is_valid = False
        hits = hits[np.argsort(np.linalg.norm(hits, axis=1))]
        self.hits = hits
        if len(hits) < config.minimum_hits:
            return
        innermost_hit = hits[0]

        # first straight 3D line through the hits -> outward xy direction and tanLambda
        hits_centre = hits.mean(axis=0)
        line_direction = np.linalg.svd(hits - hits_centre)[2][0]
        if np.dot(hits[-1] - hits[0], line_direction) < 0:
            line_direction = -line_direction
        outward_xy = line_direction[:2] / np.linalg.norm(line_direction[:2])
        self.outward_direction_xy = outward_xy
        self.normal_direction_xy = np.array([-outward_xy[1], outward_xy[0]])
        offsets_from_first_hit = hits - innermost_hit
        arc_positions = offsets_from_first_hit[:, :2] @ self.outward_direction_xy
        normal_offsets = offsets_from_first_hit[:, :2] @ self.normal_direction_xy
        first_tan_lambda = line_direction[2] / np.linalg.norm(line_direction[:2])
        cos_lambda = 1.0 / np.sqrt(1.0 + first_tan_lambda ** 2)

        # momentum for the multiple-scattering estimate
        if omega:
            transverse_momentum = 0.3 * config.field_for_momentum * 1e-3 / abs(omega)
        else:
            transverse_momentum = 10.0
        momentum = max(transverse_momentum / cos_lambda, config.minimum_momentum)
        self.momentum = momentum
        self.transverse_momentum = momentum * cos_lambda
        self.cos_lambda = cos_lambda
        # expectation value of 1/p^2 with a Gaussian momentum resolution
        self.scattering_factor = (HIGHLAND_CONSTANT / momentum) ** 2 * (1.0 + 3.0 * config.momentum_resolution ** 2)

        # ---- hit covariance with correlated multiple scattering --------------------------
        path_lengths = arc_positions / cos_lambda      # 3D path length from the first hit
        number_of_hits = len(hits)
        scattering_covariance = np.zeros((number_of_hits, number_of_hits))
        # total radiation lengths along the track (11 samples), for the logarithm term
        samples = [radiation_lengths_per_mm(innermost_hit + (hits[-1] - innermost_hit) * fraction, config)
                   for fraction in np.linspace(0, 1, 11)]
        total_radiation_lengths = sum(samples) / 11 * max(path_lengths[-1], 1.0)
        log_term = logarithm_correction(total_radiation_lengths)
        for first in range(number_of_hits):
            for second in range(first, number_of_hits):
                shorter_path = min(path_lengths[first], path_lengths[second])
                if shorter_path <= 0:
                    continue
                step_edges = np.linspace(0.0, shorter_path, config.integration_steps + 1)
                step_centres = 0.5 * (step_edges[1:] + step_edges[:-1])
                step_size = step_edges[1] - step_edges[0]
                density = np.array([
                    radiation_lengths_per_mm(
                        innermost_hit + (hits[-1] - innermost_hit) * (centre / max(path_lengths[-1], 1e-9)), config)
                    for centre in step_centres])
                value = np.sum(self.scattering_factor * log_term * density
                               * (path_lengths[first] - step_centres) * (path_lengths[second] - step_centres)) * step_size
                if config.use_material:
                    scattering_covariance[first, second] = value
                    scattering_covariance[second, first] = value
                else:
                    scattering_covariance[first, second] = 0.0
                    scattering_covariance[second, first] = 0.0
        covariance_xy = scattering_covariance + np.eye(number_of_hits) * config.hit_resolution ** 2
        covariance_z = scattering_covariance / cos_lambda ** 2 + np.eye(number_of_hits) * config.hit_resolution ** 2

        # ---- generalised least-squares fits ------------------------------------------------
        if config.fit_curvature and number_of_hits >= 3:
            number_of_parameters = 3
        else:
            number_of_parameters = 2
        design_xy = np.stack([np.ones(number_of_hits), arc_positions, 0.5 * arc_positions ** 2][:number_of_parameters],
                             axis=1)
        weight_xy = np.linalg.inv(covariance_xy)
        parameter_covariance_xy = np.linalg.inv(design_xy.T @ weight_xy @ design_xy)
        self.parameters_xy = parameter_covariance_xy @ design_xy.T @ weight_xy @ normal_offsets
        self.parameter_covariance_xy = parameter_covariance_xy
        design_z = np.stack([np.ones(number_of_hits), arc_positions], axis=1)
        weight_z = np.linalg.inv(covariance_z)
        parameter_covariance_z = np.linalg.inv(design_z.T @ weight_z @ design_z)
        self.parameters_z = parameter_covariance_z @ design_z.T @ weight_z @ offsets_from_first_hit[:, 2]
        self.parameter_covariance_z = parameter_covariance_z
        self.innermost_hit = innermost_hit
        residual_xy = normal_offsets - design_xy @ self.parameters_xy
        self.chi2 = float(residual_xy @ weight_xy @ residual_xy)
        self.is_valid = True

    def offset_and_slope(self, arc_position):
        """Offset u and slope du/ds of the exact circle with tangent slope b and curvature at (0, a)."""
        offset_at_first_hit = self.parameters_xy[0]
        slope_at_first_hit = self.parameters_xy[1]
        if len(self.parameters_xy) > 2:
            curvature = self.parameters_xy[2]
        else:
            curvature = 0.0
        if abs(curvature) < 1e-9:
            return offset_at_first_hit + slope_at_first_hit * arc_position, slope_at_first_hit
        normalisation = np.sqrt(1.0 + slope_at_first_hit * slope_at_first_hit)
        signed_radius = 1.0 / curvature
        centre_s = -slope_at_first_hit / normalisation * signed_radius
        centre_u = offset_at_first_hit + signed_radius / normalisation
        under_root = signed_radius * signed_radius - (arc_position - centre_s) ** 2
        if under_root <= 0:
            # outside the circle's range in s: parabola instead
            return (offset_at_first_hit + slope_at_first_hit * arc_position + 0.5 * curvature * arc_position * arc_position,
                    slope_at_first_hit + curvature * arc_position)
        root = np.sqrt(under_root)
        offset = centre_u - np.sign(signed_radius) * root
        slope = np.sign(signed_radius) * (arc_position - centre_s) / root
        return offset, slope

    def z_at(self, arc_position):
        return self.innermost_hit[2] + self.parameters_z[0] + self.parameters_z[1] * arc_position

    def point(self, arc_position):
        offset, _ = self.offset_and_slope(arc_position)
        xy = self.innermost_hit[:2] + arc_position * self.outward_direction_xy + offset * self.normal_direction_xy
        return np.array([xy[0], xy[1], self.z_at(arc_position)])

    def local_coordinates(self, position):
        """Arc position s and normal offset u of a point in the track's local frame."""
        offset_xy = position[:2] - self.innermost_hit[:2]
        return float(offset_xy @ self.outward_direction_xy), float(offset_xy @ self.normal_direction_xy)

    def uncertainties(self, arc_position, config):
        """(sigma perpendicular to the track in xy, sigma in z) of the trajectory at arc position s."""
        jacobian_xy = np.array([1.0, arc_position, 0.5 * arc_position * arc_position])[: len(self.parameters_xy)]
        variance_xy = float(jacobian_xy @ self.parameter_covariance_xy @ jacobian_xy)
        jacobian_z = np.array([1.0, arc_position])
        variance_z = float(jacobian_z @ self.parameter_covariance_z @ jacobian_z)
        if arc_position < 0:
            # multiple scattering between the point and the first hit (calorimeter)
            path_to_first_hit = -arc_position / self.cos_lambda
            start = self.point(arc_position)
            step_edges = np.linspace(0.0, path_to_first_hit, config.integration_steps + 1)
            step_centres = 0.5 * (step_edges[1:] + step_edges[:-1])
            step_size = step_edges[1] - step_edges[0]
            density = np.array([
                radiation_lengths_per_mm(start + (self.innermost_hit - start) * (centre / path_to_first_hit), config)
                for centre in step_centres])
            total_radiation_lengths = float(np.sum(density) * step_size)
            if total_radiation_lengths > 0:
                log_term = logarithm_correction(total_radiation_lengths)
            else:
                log_term = 1.0
            # step_centres = distance of the scatterer from the point
            if config.use_material:
                scattering = float(np.sum(self.scattering_factor * log_term * density * step_centres ** 2) * step_size)
            else:
                scattering = 0.0
            variance_xy += scattering
            variance_z += scattering / self.cos_lambda ** 2
            # curvature error from the field mismatch: sagitta 0.5 * delta(curvature) * s^2
            variance_xy += (0.5 * 0.3e-3 * config.field_uncertainty / self.transverse_momentum
                            * arc_position * arc_position) ** 2
        return np.sqrt(variance_xy), np.sqrt(variance_z)


def normalised_residuals(tracks, vertex, config, fixed_uncertainties=None):
    """Two residuals per track divided by their uncertainties; also the uncertainties and arc positions."""
    residuals = []
    uncertainties = []
    arc_positions = []
    for index, track in enumerate(tracks):
        arc_position, normal_offset = track.local_coordinates(vertex)
        offset, slope = track.offset_and_slope(arc_position)
        z_of_track = track.z_at(arc_position)
        residual_xy = (normal_offset - offset) / np.sqrt(1.0 + slope * slope)
        residual_z = vertex[2] - z_of_track
        if fixed_uncertainties is None:
            sigma_xy, sigma_z = track.uncertainties(arc_position, config)
        else:
            sigma_xy, sigma_z = fixed_uncertainties[index]
        residuals += [residual_xy / sigma_xy, residual_z / sigma_z]
        uncertainties.append((sigma_xy, sigma_z))
        arc_positions.append(arc_position)
    return np.array(residuals), uncertainties, arc_positions


def tangent_lines_closest_approach(first_track, second_track):
    """Starting vertex: closest approach of the tangent lines at the first hits."""
    def tangent_line(track):
        direction = np.array([track.outward_direction_xy[0], track.outward_direction_xy[1], track.parameters_z[1]])
        direction = direction / np.linalg.norm(direction)
        return track.point(0.0), direction
    first_point, first_direction = tangent_line(first_track)
    second_point, second_direction = tangent_line(second_track)
    between_points = first_point - second_point
    cos_angle = first_direction @ second_direction
    denominator = 1 - cos_angle * cos_angle
    if denominator < 1e-8:
        return None
    first_dot = first_direction @ between_points
    second_dot = second_direction @ between_points
    first_distance = (cos_angle * second_dot - first_dot) / denominator
    second_distance = (second_dot - cos_angle * first_dot) / denominator
    return 0.5 * (first_point + first_distance * first_direction + second_point + second_distance * second_direction)


def numerical_jacobian(tracks, vertex, config, residuals, uncertainties, step=1.0):
    """Derivative of the residuals with respect to the three vertex coordinates."""
    jacobian = np.zeros((len(residuals), 3))
    for coordinate in range(3):
        shift = np.zeros(3)
        shift[coordinate] = step
        shifted_residuals = normalised_residuals(tracks, vertex + shift, config, uncertainties)[0]
        jacobian[:, coordinate] = (shifted_residuals - residuals) / step
    return jacobian


def fit_vertex(tracks, config, maximum_iterations=12):
    vertex = tangent_lines_closest_approach(tracks[0], tracks[1])
    if vertex is None:
        return None
    # Gauss-Newton iterations
    for _ in range(maximum_iterations):
        residuals, uncertainties, _ = normalised_residuals(tracks, vertex, config)
        jacobian = numerical_jacobian(tracks, vertex, config, residuals, uncertainties)
        try:
            correction = -np.linalg.solve(jacobian.T @ jacobian, jacobian.T @ residuals)
        except np.linalg.LinAlgError:
            return None
        vertex = vertex + correction
        if np.linalg.norm(correction) < 0.1:
            break
    # covariance and chi2 at the final vertex
    residuals, uncertainties, arc_positions = normalised_residuals(tracks, vertex, config)
    jacobian = numerical_jacobian(tracks, vertex, config, residuals, uncertainties)
    try:
        covariance = np.linalg.inv(jacobian.T @ jacobian)
    except np.linalg.LinAlgError:
        return None
    chi2 = float(residuals @ residuals)
    degrees_of_freedom = 2 * len(tracks) - 3
    closest_points = [track.point(arc_position) for track, arc_position in zip(tracks, arc_positions)]
    distance_of_closest_approach = float(np.linalg.norm(closest_points[0] - closest_points[1]))
    upstream = all(arc_position < config.upstream_tolerance for arc_position in arc_positions)
    radius_xy = np.hypot(vertex[0], vertex[1])
    behind_solenoid = radius_xy > config.solenoid_radius or abs(vertex[2]) > config.solenoid_half_z
    # the decay must lie between the IP and the chambers: closer to the origin than every
    # first hit and on the same side (rejects arcs extrapolated past the IP)
    between_ip_and_hits = all(np.linalg.norm(vertex) < np.linalg.norm(track.innermost_hit)
                              and np.dot(vertex, track.innermost_hit) > 0 for track in tracks)
    return dict(vertex=vertex, covariance=covariance, chi2=chi2, degrees_of_freedom=degrees_of_freedom,
                distance_of_closest_approach=distance_of_closest_approach, upstream=upstream,
                behind_solenoid=behind_solenoid, arc_positions=arc_positions, between_ip_and_hits=between_ip_and_hits)
