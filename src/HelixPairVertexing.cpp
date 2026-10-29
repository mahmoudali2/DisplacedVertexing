/*
 * Geometry-free 3D two-track vertexing, see HelixPairVertexing.h.
 * Author: Mahmoud Althakeel
 */
#include "HelixPairVertexing.h"

#include <TMatrixDSym.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace k4vertexing::helixpair {

namespace {
constexpr double kNotANumber = std::numeric_limits<double>::quiet_NaN();

/// Closest approach of the segments firstStart-firstEnd and secondStart-secondEnd:
/// the fraction along each segment and the distance.
double segmentClosestApproach(const TVector3& firstStart, const TVector3& firstEnd, const TVector3& secondStart,
                              const TVector3& secondEnd, double& firstFraction, double& secondFraction) {
  const TVector3 firstVector = firstEnd - firstStart;
  const TVector3 secondVector = secondEnd - secondStart;
  const TVector3 betweenStarts = firstStart - secondStart;
  const double firstLength2 = firstVector.Dot(firstVector);
  const double dotProduct = firstVector.Dot(secondVector);
  const double secondLength2 = secondVector.Dot(secondVector);
  const double firstDotOffset = firstVector.Dot(betweenStarts);
  const double secondDotOffset = secondVector.Dot(betweenStarts);
  const double denominator = firstLength2 * secondLength2 - dotProduct * dotProduct;
  auto clipToSegment = [](double fraction) { return std::clamp(fraction, 0.0, 1.0); };
  if (denominator > 1e-12 * firstLength2 * secondLength2)
    secondFraction = clipToSegment((firstLength2 * secondDotOffset - dotProduct * firstDotOffset) / denominator);
  else
    secondFraction = secondLength2 > 0 ? clipToSegment(secondDotOffset / secondLength2) : 0.0;
  firstFraction = firstLength2 > 0 ? clipToSegment((dotProduct * secondFraction - firstDotOffset) / firstLength2) : 0.0;
  secondFraction = secondLength2 > 0 ? clipToSegment((dotProduct * firstFraction + secondDotOffset) / secondLength2) : 0.0;
  return (firstStart + firstFraction * firstVector - (secondStart + secondFraction * secondVector)).Mag();
}

/// Trajectory point closest to a position (searched along the propagated polyline) and
/// the direction of the trajectory there.
void closestPointOnTrajectory(const HelixTrack::Trajectory& trajectory, const TVector3& position,
                              TVector3& closestPoint, TVector3& direction) {
  double bestDistance = std::numeric_limits<double>::infinity();
  for (std::size_t index = 0; index + 1 < trajectory.position.size(); ++index) {
    const TVector3 segment = trajectory.position[index + 1] - trajectory.position[index];
    const double segmentLength2 = segment.Mag2();
    double fraction = segmentLength2 > 0 ? (position - trajectory.position[index]).Dot(segment) / segmentLength2 : 0.0;
    fraction = std::clamp(fraction, 0.0, 1.0);
    const TVector3 candidate = trajectory.position[index] + fraction * segment;
    const double distance = (candidate - position).Mag();
    if (distance < bestDistance) {
      bestDistance = distance;
      closestPoint = candidate;
      direction = segment.Unit();
    }
  }
}

/// Covariance [mm^2] of the point of a track closest to a position, from the covariance
/// of its track state: numerical Jacobian (central differences, step = half a sigma of
/// each parameter) through the same build + propagation as the track itself.
/// Also returns the closest point and the track direction there.
bool trackPointCovariance(const HelixTrack& track, const TVector3& position, const Config& config,
                          double pointCovariance[3][3], TVector3& closestPoint, TVector3& direction) {
  const TrackStateParameters& state = track.state();
  if (!state.hasCovariance)
    return false;
  closestPointOnTrajectory(track.propagate(config), position, closestPoint, direction);
  // the rebuilt tracks start at the same point as the original one
  const std::vector<TVector3> startPoint{track.startPoint()};

  // derivative of the closest point with respect to each of the five parameters
  double derivative[5][3] = {};
  for (int parameter = 0; parameter < 5; ++parameter) {
    const double sigma = std::sqrt(std::max(state.covariance[parameter][parameter], 0.0));
    if (sigma == 0)
      continue;
    const double stepSize = 0.5 * sigma;
    TVector3 shiftedPoints[2];
    for (int side = 0; side < 2; ++side) {
      TrackStateParameters changedState = state;
      double* value[5] = {&changedState.d0, &changedState.phi, &changedState.omega, &changedState.z0,
                          &changedState.tanLambda};
      *value[parameter] += (side == 0 ? stepSize : -stepSize);
      HelixTrack rebuiltTrack;
      if (!rebuiltTrack.build(changedState, startPoint, config))
        return false;
      TVector3 unusedDirection;
      closestPointOnTrajectory(rebuiltTrack.propagate(config), position, shiftedPoints[side], unusedDirection);
    }
    const TVector3 slope = (shiftedPoints[0] - shiftedPoints[1]) * (1.0 / (2.0 * stepSize));
    derivative[parameter][0] = slope.X();
    derivative[parameter][1] = slope.Y();
    derivative[parameter][2] = slope.Z();
  }

  // point covariance = derivative^T x state covariance x derivative
  for (int row = 0; row < 3; ++row)
    for (int column = 0; column < 3; ++column) {
      double sum = 0;
      for (int first = 0; first < 5; ++first)
        for (int second = 0; second < 5; ++second)
          sum += derivative[first][row] * state.covariance[first][second] * derivative[second][column];
      pointCovariance[row][column] = sum;
    }
  return true;
}

/// Weight matrix of a track point: the inverse of its covariance in the plane
/// perpendicular to the track direction (the point is free along the track).
bool perpendicularWeight(const double pointCovariance[3][3], const TVector3& direction, double weight[3][3]) {
  // two unit vectors perpendicular to the track
  const TVector3 firstPerpendicular = direction.Orthogonal().Unit();
  const TVector3 secondPerpendicular = direction.Cross(firstPerpendicular).Unit();
  const TVector3 basis[2] = {firstPerpendicular, secondPerpendicular};
  // covariance in the perpendicular plane (2 x 2)
  double planeCovariance[2][2];
  for (int first = 0; first < 2; ++first)
    for (int second = 0; second < 2; ++second) {
      double sum = 0;
      for (int row = 0; row < 3; ++row)
        for (int column = 0; column < 3; ++column)
          sum += basis[first][row] * pointCovariance[row][column] * basis[second][column];
      planeCovariance[first][second] = sum;
    }
  const double determinant =
      planeCovariance[0][0] * planeCovariance[1][1] - planeCovariance[0][1] * planeCovariance[1][0];
  if (!(determinant > 0))
    return false;
  const double inverse[2][2] = {{planeCovariance[1][1] / determinant, -planeCovariance[0][1] / determinant},
                                {-planeCovariance[1][0] / determinant, planeCovariance[0][0] / determinant}};
  // back to 3 x 3
  for (int row = 0; row < 3; ++row)
    for (int column = 0; column < 3; ++column) {
      double sum = 0;
      for (int first = 0; first < 2; ++first)
        for (int second = 0; second < 2; ++second)
          sum += basis[first][row] * inverse[first][second] * basis[second][column];
      weight[row][column] = sum;
    }
  return true;
}

/// Linearised vertex fit of two tracks starting from the closest-approach midpoint
/// (three iterations). Each track contributes its closest point and weight; the vertex
/// minimises sum (vertex - point)^T weight (vertex - point).
/// Returns false if a covariance is missing or singular.
bool fitWithCovariance(const HelixTrack& firstTrack, const HelixTrack& secondTrack, const Config& config,
                       VertexResult& result) {
  TVector3 vertex = result.closestApproachMidpoint;
  for (int iteration = 0; iteration < 3; ++iteration) {
    TMatrixDSym sumOfWeights(3);
    double sumOfWeightedPoints[3] = {0, 0, 0};
    TVector3 closestPoints[2];
    double weights[2][3][3];
    int trackIndex = 0;
    for (const HelixTrack* track : {&firstTrack, &secondTrack}) {
      double pointCovariance[3][3];
      TVector3 direction;
      if (!trackPointCovariance(*track, vertex, config, pointCovariance, closestPoints[trackIndex], direction))
        return false;
      if (!perpendicularWeight(pointCovariance, direction, weights[trackIndex]))
        return false;
      for (int row = 0; row < 3; ++row)
        for (int column = 0; column < 3; ++column) {
          sumOfWeights(row, column) += weights[trackIndex][row][column];
          sumOfWeightedPoints[row] += weights[trackIndex][row][column] * closestPoints[trackIndex][column];
        }
      ++trackIndex;
    }
    double determinant = 0;
    TMatrixDSym vertexCovariance(sumOfWeights);
    vertexCovariance.Invert(&determinant);
    if (!(determinant > 0))
      return false;
    TVector3 fittedVertex;
    for (int row = 0; row < 3; ++row) {
      double sum = 0;
      for (int column = 0; column < 3; ++column)
        sum += vertexCovariance(row, column) * sumOfWeightedPoints[column];
      fittedVertex[row] = sum;
    }
    double chi2 = 0;
    for (int track = 0; track < 2; ++track) {
      const TVector3 residual = fittedVertex - closestPoints[track];
      for (int row = 0; row < 3; ++row)
        for (int column = 0; column < 3; ++column)
          chi2 += residual[row] * weights[track][row][column] * residual[column];
    }
    vertex = fittedVertex;
    result.chi2 = chi2;
    for (int row = 0; row < 3; ++row)
      for (int column = 0; column < 3; ++column)
        result.covariance[row][column] = vertexCovariance(row, column);
  }
  result.position = vertex;
  result.isFitted = true;
  return true;
}
} // namespace

double fieldZAt(const TVector3& position, const Config& config) {
  const bool insideFieldCylinder =
      position.Perp() < config.fieldRadius && std::fabs(position.Z()) < config.fieldHalfLength;
  return insideFieldCylinder ? config.innerField : config.outerField;
}

double energyLossPerMmAt(const TVector3& position, const Config& config) {
  const double radiusXY = position.Perp();
  const double absoluteZ = std::fabs(position.Z());
  const bool insideOuterBox = radiusXY < config.shellOuter && absoluteZ < config.shellOuter;
  const bool insideInnerBox = radiusXY < config.shellInner && absoluteZ < config.shellInner;
  // in the shell (outside the beam hole)
  if (insideOuterBox && !insideInnerBox && radiusXY > config.shellBeamHole)
    return config.energyLossInShell;
  // beyond the shell
  if (!insideOuterBox)
    return config.energyLossOutside;
  // inside the inner box or in the beam hole: no energy loss
  return 0.0;
}

bool HelixTrack::build(const TrackStateParameters& state, const std::vector<TVector3>& hits, const Config& config) {
  // the state describes the track; the hits (possibly none, e.g. a drift-chamber
  // track without silicon hits) only fix where the trajectory starts
  m_state = state;
  m_numberOfHits = static_cast<int>(hits.size());

  // point of closest approach to the reference point, EDM4hep convention
  const double sinPhi = std::sin(state.phi);
  const double cosPhi = std::cos(state.phi);
  const TVector3 closestApproach(state.referencePoint.X() - state.d0 * sinPhi,
                                 state.referencePoint.Y() + state.d0 * cosPhi, state.referencePoint.Z() + state.z0);

  // the trajectory starts at the innermost hit (closest to the origin)
  TVector3 startPoint = closestApproach;
  if (!hits.empty()) {
    startPoint = hits.front();
    for (const auto& hit : hits)
      if (hit.Mag() < startPoint.Mag())
        startPoint = hit;
  }
  m_startPoint = startPoint;
  m_fieldAtStart = fieldZAt(m_startPoint, config);

  // ---- straight line (no curvature or no field) ------------------------------
  if (state.omega == 0 || m_fieldAtStart == 0) {
    m_isStraightLine = true;
    m_lineDirection = TVector3(cosPhi, sinPhi, state.tanLambda).Unit();
    // the line starts at its point closest to the start point
    m_lineStart = closestApproach + (m_startPoint - closestApproach).Dot(m_lineDirection) * m_lineDirection;
    m_radius = std::numeric_limits<double>::infinity();
    m_zPerRadian = 0;
    m_rotationSense = 0;
    m_measuredMomentum = kNotANumber;
    return true;
  }

  // ---- helix ------------------------------------------------------------------
  // omega > 0: the circle centre lies to the right of the direction of motion,
  // i.e. the track turns clockwise (rotation sense -1)
  const double signedRadius = 1.0 / state.omega;
  m_radius = std::fabs(signedRadius);
  m_centreX = closestApproach.X() + signedRadius * sinPhi;
  m_centreY = closestApproach.Y() - signedRadius * cosPhi;
  m_rotationSense = state.omega > 0 ? -1.0 : 1.0;
  // z changes by tanLambda * radius per radian of turning
  m_zPerRadian = state.tanLambda * m_radius;
  m_pathLengthPerRadian = std::sqrt(m_radius * m_radius + m_zPerRadian * m_zPerRadian);

  // angle of the start point around the circle centre, and z there
  const double angleAtClosestApproach = std::atan2(closestApproach.Y() - m_centreY, closestApproach.X() - m_centreX);
  const double angleAtStart = std::atan2(m_startPoint.Y() - m_centreY, m_startPoint.X() - m_centreX);
  // turning angle from the closest approach to the start point, in the direction of motion
  double turningAngle = m_rotationSense * (angleAtStart - angleAtClosestApproach);
  while (turningAngle > M_PI)
    turningAngle -= 2 * M_PI;
  while (turningAngle < -M_PI)
    turningAngle += 2 * M_PI;
  m_angleAtStart = angleAtStart;
  m_zAtStart = closestApproach.Z() + m_zPerRadian * turningAngle;
  m_isStraightLine = false;

  // momentum from the curvature in the field at the start point
  const double transverseMomentum = 0.3 * std::fabs(m_fieldAtStart) * 1e-3 / std::fabs(state.omega);
  m_measuredMomentum = transverseMomentum * std::sqrt(1.0 + state.tanLambda * state.tanLambda);
  return true;
}

TVector3 HelixTrack::point(double pathLength) const {
  if (m_isStraightLine)
    return m_lineStart + pathLength * m_lineDirection;
  const double turningAngle = pathLength / m_pathLengthPerRadian;
  const double angle = m_angleAtStart + m_rotationSense * turningAngle;
  return TVector3(m_centreX + m_radius * std::cos(angle), m_centreY + m_radius * std::sin(angle),
                  m_zAtStart + m_zPerRadian * turningAngle);
}

TVector3 HelixTrack::tangent(double pathLength) const {
  if (m_isStraightLine)
    return m_lineDirection;
  const double angle = m_angleAtStart + m_rotationSense * pathLength / m_pathLengthPerRadian;
  return TVector3(-m_radius * std::sin(angle) * m_rotationSense, m_radius * std::cos(angle) * m_rotationSense,
                  m_zPerRadian)
      .Unit();
}

HelixTrack::Trajectory HelixTrack::propagate(const Config& config) const {
  Trajectory trajectory;
  // fitted helix up to the tolerance outward (path lengths tolerance, tolerance - step, ... > 0,
  // in increasing order)
  std::vector<double> outwardPathLengths;
  for (double pathLength = config.upstreamTolerance; pathLength > 0; pathLength -= config.stepLength)
    outwardPathLengths.push_back(pathLength);
  std::reverse(outwardPathLengths.begin(), outwardPathLengths.end());

  // inward stepping
  const int numberOfSteps = static_cast<int>(config.backwardFactor * m_startPoint.Mag() / config.stepLength);
  TVector3 position = point(0.0);
  TVector3 direction = tangent(0.0);
  const double cosLambda = std::hypot(direction.X(), direction.Y());
  const double tanLambda = direction.Z() / std::max(cosLambda, 1e-12);
  double azimuth = std::atan2(direction.Y(), direction.X());
  double momentum = m_measuredMomentum;
  const double signOfStartField = m_fieldAtStart < 0 ? -1.0 : 1.0;
  std::vector<TVector3> inwardPositions{position};
  std::vector<TVector3> inwardDirections{direction};
  std::vector<double> inwardMomenta{momentum};
  for (int step = 0; step < numberOfSteps; ++step) {
    if (!m_isStraightLine) {
      const double transverseMomentum = momentum * cosLambda;
      // signed xy curvature of the outward motion: rotation sense / radius where the hits
      // are, scales with B / pT
      const double curvature =
          m_rotationSense * signOfStartField * fieldZAt(position, config) * 0.3e-3 / std::max(transverseMomentum, 1e-3);
      position -= config.stepLength *
                  TVector3(cosLambda * std::cos(azimuth), cosLambda * std::sin(azimuth), cosLambda * tanLambda);
      azimuth -= curvature * config.stepLength * cosLambda;
      momentum += energyLossPerMmAt(position, config) * config.stepLength;
      direction = TVector3(cosLambda * std::cos(azimuth), cosLambda * std::sin(azimuth), cosLambda * tanLambda);
    } else {
      position -= config.stepLength * direction;
    }
    inwardPositions.push_back(position);
    inwardDirections.push_back(direction);
    inwardMomenta.push_back(momentum);
  }

  // innermost point first, then outward
  const std::size_t numberOfInwardPoints = inwardPositions.size();
  trajectory.pathLength.reserve(numberOfInwardPoints + outwardPathLengths.size());
  for (std::size_t index = 0; index < numberOfInwardPoints; ++index) {
    const std::size_t reversedIndex = numberOfInwardPoints - 1 - index;
    trajectory.pathLength.push_back(-config.stepLength * static_cast<double>(reversedIndex));
    trajectory.position.push_back(inwardPositions[reversedIndex]);
    trajectory.direction.push_back(inwardDirections[reversedIndex]);
    trajectory.momentum.push_back(inwardMomenta[reversedIndex]);
  }
  for (double pathLength : outwardPathLengths) {
    trajectory.pathLength.push_back(pathLength);
    trajectory.position.push_back(point(pathLength));
    trajectory.direction.push_back(tangent(pathLength));
    trajectory.momentum.push_back(m_measuredMomentum);
  }
  return trajectory;
}

std::optional<VertexResult> fitVertex(const HelixTrack& firstTrack, const HelixTrack& secondTrack,
                                      const Config& config) {
  const auto firstTrajectory = firstTrack.propagate(config);
  const auto secondTrajectory = secondTrack.propagate(config);
  const std::size_t firstSize = firstTrajectory.position.size();
  const std::size_t secondSize = secondTrajectory.position.size();
  if (firstSize < 2 || secondSize < 2)
    return std::nullopt;

  // coarse search: the closest pair of trajectory points
  std::size_t firstIndex = 0, secondIndex = 0;
  double bestDistance2 = std::numeric_limits<double>::infinity();
  for (std::size_t first = 0; first < firstSize; ++first)
    for (std::size_t second = 0; second < secondSize; ++second) {
      const double distance2 = (firstTrajectory.position[first] - secondTrajectory.position[second]).Mag2();
      if (distance2 < bestDistance2) {
        bestDistance2 = distance2;
        firstIndex = first;
        secondIndex = second;
      }
    }

  // fine search: closest approach of the neighbouring segments
  double distanceOfClosestApproach = std::numeric_limits<double>::infinity();
  double firstFraction = 0, secondFraction = 0;
  std::size_t firstSegment = 0, secondSegment = 0;
  for (long firstCandidate : {static_cast<long>(firstIndex) - 1, static_cast<long>(firstIndex)})
    for (long secondCandidate : {static_cast<long>(secondIndex) - 1, static_cast<long>(secondIndex)}) {
      const bool outsideFirst = firstCandidate < 0 || firstCandidate >= static_cast<long>(firstSize) - 1;
      const bool outsideSecond = secondCandidate < 0 || secondCandidate >= static_cast<long>(secondSize) - 1;
      if (outsideFirst || outsideSecond)
        continue;
      double candidateFirstFraction, candidateSecondFraction;
      const double distance = segmentClosestApproach(
          firstTrajectory.position[firstCandidate], firstTrajectory.position[firstCandidate + 1],
          secondTrajectory.position[secondCandidate], secondTrajectory.position[secondCandidate + 1],
          candidateFirstFraction, candidateSecondFraction);
      if (distance < distanceOfClosestApproach) {
        distanceOfClosestApproach = distance;
        firstSegment = firstCandidate;
        secondSegment = secondCandidate;
        firstFraction = candidateFirstFraction;
        secondFraction = candidateSecondFraction;
      }
    }
  if (!std::isfinite(distanceOfClosestApproach)) {
    distanceOfClosestApproach = std::sqrt(bestDistance2);
    firstSegment = std::min(firstIndex, firstSize - 2);
    secondSegment = std::min(secondIndex, secondSize - 2);
    firstFraction = secondFraction = 0;
  }

  const TVector3 firstClosest =
      firstTrajectory.position[firstSegment] +
      firstFraction * (firstTrajectory.position[firstSegment + 1] - firstTrajectory.position[firstSegment]);
  const TVector3 secondClosest =
      secondTrajectory.position[secondSegment] +
      secondFraction * (secondTrajectory.position[secondSegment + 1] - secondTrajectory.position[secondSegment]);
  VertexResult result;
  result.closestApproachMidpoint = 0.5 * (firstClosest + secondClosest);
  result.position = result.closestApproachMidpoint;
  result.pairDistanceOfClosestApproach = distanceOfClosestApproach;
  result.pathLengthFromStart[0] =
      firstTrajectory.pathLength[firstSegment] +
      firstFraction * (firstTrajectory.pathLength[firstSegment + 1] - firstTrajectory.pathLength[firstSegment]);
  result.pathLengthFromStart[1] =
      secondTrajectory.pathLength[secondSegment] +
      secondFraction * (secondTrajectory.pathLength[secondSegment + 1] - secondTrajectory.pathLength[secondSegment]);
  result.upstream = result.pathLengthFromStart[0] < config.upstreamTolerance &&
                    result.pathLengthFromStart[1] < config.upstreamTolerance;
  // closer to the origin than every start point, on the same side
  result.betweenIPAndHits = true;
  for (const HelixTrack* track : {&firstTrack, &secondTrack}) {
    const bool closerThanStart = result.closestApproachMidpoint.Mag() < track->startPoint().Mag();
    const bool sameSide = result.closestApproachMidpoint.Dot(track->startPoint()) > 0;
    if (!(closerThanStart && sameSide))
      result.betweenIPAndHits = false;
  }

  // pair kinematics at the vertex
  const TVector3& firstDirection = firstTrajectory.direction[firstSegment];
  const TVector3& secondDirection = secondTrajectory.direction[secondSegment];
  if (firstTrack.isStraightLine() || secondTrack.isStraightLine()) {
    result.mass = kNotANumber;
    result.cosOpening = firstDirection.Dot(secondDirection);
    result.momentumAtVertex[0] = result.momentumAtVertex[1] = kNotANumber;
  } else {
    const TVector3 firstMomentum = firstTrajectory.momentum[firstSegment] * firstDirection;
    const TVector3 secondMomentum = secondTrajectory.momentum[secondSegment] * secondDirection;
    const double mass2 = config.particleMass * config.particleMass;
    const double energy = std::sqrt(firstMomentum.Mag2() + mass2) + std::sqrt(secondMomentum.Mag2() + mass2);
    result.mass = std::sqrt(std::max(energy * energy - (firstMomentum + secondMomentum).Mag2(), 0.0));
    result.cosOpening = firstMomentum.Dot(secondMomentum) / (firstMomentum.Mag() * secondMomentum.Mag());
    result.momentumAtVertex[0] = firstTrajectory.momentum[firstSegment];
    result.momentumAtVertex[1] = secondTrajectory.momentum[secondSegment];
  }

  // vertex fit with the track-state covariances (keeps the closest-approach midpoint otherwise)
  if (firstTrack.state().hasCovariance && secondTrack.state().hasCovariance)
    fitWithCovariance(firstTrack, secondTrack, config, result);
  return result;
}

} // namespace k4vertexing::helixpair
