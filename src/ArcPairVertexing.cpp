/*
 * Two-track vertexing with local arc fits, see ArcPairVertexing.h.
 * Author: Mahmoud Althakeel
 */
#include "ArcPairVertexing.h"

#include <TMatrixDSymEigen.h>
#include <TVectorD.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace k4vertexing::arcpair {

namespace {
constexpr double kHighlandConstant = 13.6e-3; // GeV, multiple-scattering formula

/// (1 + 0.038 ln x)^2 term of the Highland formula, for x radiation lengths
double logarithmTerm(double radiationLengths) {
  return radiationLengths > 0 ? std::pow(1.0 + 0.038 * std::log(std::max(radiationLengths, 1e-3)), 2) : 1.0;
}

/// Generalised least squares: fitted parameters and their covariance.
bool generalisedLeastSquares(const TMatrixD& designMatrix, const TMatrixDSym& measurementCovariance,
                             const TVectorD& measurements, TVectorD& parameters, TMatrixDSym& parameterCovariance) {
  TMatrixDSym weight(measurementCovariance);
  double determinant = 0;
  weight.Invert(&determinant);
  if (determinant == 0)
    return false;
  TMatrixDSym normalMatrix(weight);
  normalMatrix.SimilarityT(designMatrix); // design^T x weight x design
  normalMatrix.Invert(&determinant);
  if (determinant == 0)
    return false;
  parameterCovariance.ResizeTo(normalMatrix);
  parameterCovariance = normalMatrix;
  TMatrixD designTransposed(TMatrixD::kTransposed, designMatrix);
  parameters.ResizeTo(designMatrix.GetNcols());
  parameters = TMatrixD(parameterCovariance, TMatrixD::kMult, TMatrixD(designTransposed, TMatrixD::kMult, weight)) *
               measurements;
  return true;
}
} // namespace

double radiationLengthsPerMm(const TVector3& position, const Config& config) {
  const double radiusXY = position.Perp();
  const double absoluteZ = std::fabs(position.Z());
  const bool insideOuterBox = radiusXY < config.shellOuterRadius && absoluteZ < config.shellOuterHalfLength;
  const bool insideInnerBox = radiusXY < config.shellInnerRadius && absoluteZ < config.shellInnerHalfLength;
  // in the shell: its radiation lengths spread over its radial depth
  if (insideOuterBox && !insideInnerBox)
    return config.shellRadiationLengths > 0
               ? config.shellRadiationLengths / (config.shellOuterRadius - config.shellInnerRadius)
               : 0.0;
  // beyond the shell: uniform density
  if (!insideOuterBox)
    return config.outsideRadiationLengths > 0 ? config.outsideRadiationLengths / config.outsideThickness : 0.0;
  // inside the inner box: no material
  return 0.0;
}

bool ArcTrack::build(std::vector<TVector3> hits, double omega, const Config& config) {
  if (hits.size() < config.minimumHits)
    return false;
  // order the hits from the IP outward
  std::sort(hits.begin(), hits.end(),
            [](const TVector3& left, const TVector3& right) { return left.Mag() < right.Mag(); });
  const int numberOfHits = static_cast<int>(hits.size());
  m_numberOfHits = numberOfHits;
  m_firstHit = hits.front();

  // first straight line: principal axis of the hits, oriented outward
  TVector3 hitsCentre(0, 0, 0);
  for (const auto& hit : hits)
    hitsCentre += hit;
  hitsCentre *= 1.0 / numberOfHits;
  TMatrixDSym scatterMatrix(3);
  for (const auto& hit : hits) {
    const TVector3 offset = hit - hitsCentre;
    for (int row = 0; row < 3; ++row)
      for (int column = 0; column < 3; ++column)
        scatterMatrix(row, column) += offset[row] * offset[column];
  }
  const TMatrixD eigenVectors = TMatrixDSymEigen(scatterMatrix).GetEigenVectors();
  TVector3 principalAxis(eigenVectors(0, 0), eigenVectors(1, 0), eigenVectors(2, 0));
  if ((hits.back() - hits.front()).Dot(principalAxis) < 0)
    principalAxis = -principalAxis;
  const double axisLengthXY = std::hypot(principalAxis.X(), principalAxis.Y());
  if (axisLengthXY < 1e-6)
    return false;
  m_outwardDirection[0] = principalAxis.X() / axisLengthXY;
  m_outwardDirection[1] = principalAxis.Y() / axisLengthXY;
  m_normalDirection[0] = -m_outwardDirection[1];
  m_normalDirection[1] = m_outwardDirection[0];
  const double tanLambda = principalAxis.Z() / axisLengthXY;
  m_cosLambda = 1.0 / std::sqrt(1.0 + tanLambda * tanLambda);

  // momentum for the scattering model; expectation of 1/p^2 for a Gaussian resolution
  const double transverseMomentum =
      omega != 0 ? 0.3 * config.fieldForMomentum * 1e-3 / std::fabs(omega) : config.defaultTransverseMomentum;
  m_momentum = std::max(transverseMomentum / m_cosLambda, config.minimumMomentum);
  m_transverseMomentum = m_momentum * m_cosLambda;
  m_scatteringFactor = std::pow(kHighlandConstant / m_momentum, 2) *
                       (1.0 + 3.0 * config.momentumResolution * config.momentumResolution);

  // local coordinates of the hits
  std::vector<double> arcPosition(numberOfHits), normalOffset(numberOfHits), zOffset(numberOfHits),
      pathLength(numberOfHits);
  for (int hit = 0; hit < numberOfHits; ++hit) {
    const TVector3 offset = hits[hit] - m_firstHit;
    arcPosition[hit] = offset.X() * m_outwardDirection[0] + offset.Y() * m_outwardDirection[1];
    normalOffset[hit] = offset.X() * m_normalDirection[0] + offset.Y() * m_normalDirection[1];
    zOffset[hit] = offset.Z();
    pathLength[hit] = arcPosition[hit] / m_cosLambda;
  }

  // hit covariance: resolution + correlated scattering between the hits
  const double lastPathLength = std::max(pathLength[numberOfHits - 1], 1e-9);
  const TVector3 chord = hits.back() - m_firstHit;
  double sumOfDensities = 0;
  for (int sample = 0; sample <= 10; ++sample)
    sumOfDensities += radiationLengthsPerMm(m_firstHit + chord * (sample / 10.0), config);
  const double logFactor = logarithmTerm(sumOfDensities / 11.0 * std::max(pathLength[numberOfHits - 1], 1.0));
  TMatrixDSym scatteringCovariance(numberOfHits);
  for (int first = 0; first < numberOfHits; ++first)
    for (int second = first; second < numberOfHits; ++second) {
      const double shorterPath = std::min(pathLength[first], pathLength[second]);
      if (shorterPath <= 0)
        continue;
      const double stepSize = shorterPath / config.integrationSteps;
      double sum = 0;
      for (int step = 0; step < config.integrationSteps; ++step) {
        const double scattererPath = (step + 0.5) * stepSize;
        const double density = radiationLengthsPerMm(m_firstHit + chord * (scattererPath / lastPathLength), config);
        sum += density * (pathLength[first] - scattererPath) * (pathLength[second] - scattererPath);
      }
      scatteringCovariance(first, second) = scatteringCovariance(second, first) =
          config.useMaterial ? m_scatteringFactor * logFactor * sum * stepSize : 0.0;
    }
  const double hitVariance = config.hitResolution * config.hitResolution;
  TMatrixDSym covarianceXY(scatteringCovariance), covarianceZ(numberOfHits);
  for (int row = 0; row < numberOfHits; ++row)
    for (int column = 0; column < numberOfHits; ++column)
      covarianceZ(row, column) = scatteringCovariance(row, column) / (m_cosLambda * m_cosLambda);
  for (int hit = 0; hit < numberOfHits; ++hit) {
    covarianceXY(hit, hit) += hitVariance;
    covarianceZ(hit, hit) += hitVariance;
  }

  // xy: u(s) = a + b s (+ curvature s^2 / 2);  z: dz(s) = c + d s
  const int numberOfParameters = (config.fitCurvature && numberOfHits >= 3) ? 3 : 2;
  TMatrixD designXY(numberOfHits, numberOfParameters), designZ(numberOfHits, 2);
  TVectorD measurementsXY(numberOfHits), measurementsZ(numberOfHits);
  for (int hit = 0; hit < numberOfHits; ++hit) {
    designXY(hit, 0) = 1;
    designXY(hit, 1) = arcPosition[hit];
    if (numberOfParameters == 3)
      designXY(hit, 2) = 0.5 * arcPosition[hit] * arcPosition[hit];
    designZ(hit, 0) = 1;
    designZ(hit, 1) = arcPosition[hit];
    measurementsXY[hit] = normalOffset[hit];
    measurementsZ[hit] = zOffset[hit];
  }
  TVectorD parametersXY, parametersZ;
  m_parameterCovarianceXY.ResizeTo(numberOfParameters, numberOfParameters);
  const bool fitXYWorked =
      generalisedLeastSquares(designXY, covarianceXY, measurementsXY, parametersXY, m_parameterCovarianceXY);
  if (!fitXYWorked ||
      !generalisedLeastSquares(designZ, covarianceZ, measurementsZ, parametersZ, m_parameterCovarianceZ))
    return false;
  m_parametersXY.assign(parametersXY.GetMatrixArray(), parametersXY.GetMatrixArray() + numberOfParameters);
  m_parametersZ[0] = parametersZ[0];
  m_parametersZ[1] = parametersZ[1];

  TVectorD residualsXY = measurementsXY - designXY * parametersXY;
  TMatrixDSym weightXY(covarianceXY);
  weightXY.Invert();
  m_chi2 = weightXY.Similarity(residualsXY);
  return true;
}

void ArcTrack::offsetAndSlope(double arcPosition, double& offset, double& slope) const {
  const double offsetAtFirstHit = m_parametersXY[0];
  const double slopeAtFirstHit = m_parametersXY[1];
  const double curvature = m_parametersXY.size() > 2 ? m_parametersXY[2] : 0.0;
  if (std::fabs(curvature) < 1e-9) {
    offset = offsetAtFirstHit + slopeAtFirstHit * arcPosition;
    slope = slopeAtFirstHit;
    return;
  }
  // exact circle with tangent slope b and the curvature at (0, a)
  const double normalisation = std::sqrt(1.0 + slopeAtFirstHit * slopeAtFirstHit);
  const double signedRadius = 1.0 / curvature;
  const double centreS = -slopeAtFirstHit / normalisation * signedRadius;
  const double centreU = offsetAtFirstHit + signedRadius / normalisation;
  const double underRoot = signedRadius * signedRadius - (arcPosition - centreS) * (arcPosition - centreS);
  if (underRoot <= 0) {
    // outside the circle's range in s: parabola instead
    offset = offsetAtFirstHit + slopeAtFirstHit * arcPosition + 0.5 * curvature * arcPosition * arcPosition;
    slope = slopeAtFirstHit + curvature * arcPosition;
    return;
  }
  const double root = std::sqrt(underRoot);
  const double radiusSign = signedRadius > 0 ? 1.0 : -1.0;
  offset = centreU - radiusSign * root;
  slope = radiusSign * (arcPosition - centreS) / root;
}

TVector3 ArcTrack::point(double arcPosition) const {
  double offset, slope;
  offsetAndSlope(arcPosition, offset, slope);
  return TVector3(m_firstHit.X() + arcPosition * m_outwardDirection[0] + offset * m_normalDirection[0],
                  m_firstHit.Y() + arcPosition * m_outwardDirection[1] + offset * m_normalDirection[1],
                  zAt(arcPosition));
}

void ArcTrack::localCoordinates(const TVector3& position, double& arcPosition, double& normalOffset) const {
  const double offsetX = position.X() - m_firstHit.X();
  const double offsetY = position.Y() - m_firstHit.Y();
  arcPosition = offsetX * m_outwardDirection[0] + offsetY * m_outwardDirection[1];
  normalOffset = offsetX * m_normalDirection[0] + offsetY * m_normalDirection[1];
}

void ArcTrack::uncertainties(double arcPosition, const Config& config, double& sigmaXY, double& sigmaZ) const {
  const int numberOfParameters = static_cast<int>(m_parametersXY.size());
  const double jacobian[3] = {1.0, arcPosition, 0.5 * arcPosition * arcPosition};
  double varianceXY = 0;
  for (int row = 0; row < numberOfParameters; ++row)
    for (int column = 0; column < numberOfParameters; ++column)
      varianceXY += jacobian[row] * m_parameterCovarianceXY(row, column) * jacobian[column];
  double varianceZ = m_parameterCovarianceZ(0, 0) + 2 * arcPosition * m_parameterCovarianceZ(0, 1) +
                     arcPosition * arcPosition * m_parameterCovarianceZ(1, 1);
  if (arcPosition < 0) {
    // scattering between the extrapolation point and the first hit
    const double pathToFirstHit = -arcPosition / m_cosLambda;
    const TVector3 start = point(arcPosition);
    const TVector3 towardsFirstHit = m_firstHit - start;
    const double stepSize = pathToFirstHit / config.integrationSteps;
    double totalRadiationLengths = 0, sum = 0;
    for (int step = 0; step < config.integrationSteps; ++step) {
      const double scattererDistance = (step + 0.5) * stepSize; // distance of the scatterer from the point
      const double density = radiationLengthsPerMm(start + towardsFirstHit * (scattererDistance / pathToFirstHit), config);
      totalRadiationLengths += density * stepSize;
      sum += density * scattererDistance * scattererDistance * stepSize;
    }
    const double scattering = config.useMaterial ? m_scatteringFactor * logarithmTerm(totalRadiationLengths) * sum : 0.0;
    varianceXY += scattering;
    varianceZ += scattering / (m_cosLambda * m_cosLambda);
    // curvature error from the field mismatch: sagitta 0.5 * delta(curvature) * s^2
    const double sagitta = 0.5 * 0.3e-3 * config.fieldUncertainty / m_transverseMomentum * arcPosition * arcPosition;
    varianceXY += sagitta * sagitta;
  }
  sigmaXY = std::sqrt(std::max(varianceXY, 0.0));
  sigmaZ = std::sqrt(std::max(varianceZ, 0.0));
}

namespace {
/// normalised residuals (two per track), the uncertainties and the arc position of the vertex
struct Residuals {
  std::vector<double> normalised, sigmaXY, sigmaZ, arcPosition;
};

void computeResiduals(const std::vector<const ArcTrack*>& tracks, const TVector3& vertex, const Config& config,
                      Residuals& residuals, bool updateUncertainties) {
  const std::size_t numberOfTracks = tracks.size();
  residuals.normalised.resize(2 * numberOfTracks);
  residuals.arcPosition.resize(numberOfTracks);
  if (updateUncertainties) {
    residuals.sigmaXY.resize(numberOfTracks);
    residuals.sigmaZ.resize(numberOfTracks);
  }
  for (std::size_t track = 0; track < numberOfTracks; ++track) {
    double arcPosition, normalOffset, offsetOfTrack, slope;
    tracks[track]->localCoordinates(vertex, arcPosition, normalOffset);
    tracks[track]->offsetAndSlope(arcPosition, offsetOfTrack, slope);
    if (updateUncertainties)
      tracks[track]->uncertainties(arcPosition, config, residuals.sigmaXY[track], residuals.sigmaZ[track]);
    residuals.normalised[2 * track] =
        (normalOffset - offsetOfTrack) / std::sqrt(1.0 + slope * slope) / residuals.sigmaXY[track];
    residuals.normalised[2 * track + 1] = (vertex.Z() - tracks[track]->zAt(arcPosition)) / residuals.sigmaZ[track];
    residuals.arcPosition[track] = arcPosition;
  }
}

/// Starting vertex: closest approach of the tangent lines at the first hits.
std::optional<TVector3> tangentLinesClosestApproach(const ArcTrack& firstTrack, const ArcTrack& secondTrack) {
  auto tangentLine = [](const ArcTrack& track, TVector3& linePoint, TVector3& lineDirection) {
    linePoint = track.point(0.0);
    lineDirection = TVector3(track.outwardDirectionX(), track.outwardDirectionY(), track.zSlope()).Unit();
  };
  TVector3 firstPoint, firstDirection, secondPoint, secondDirection;
  tangentLine(firstTrack, firstPoint, firstDirection);
  tangentLine(secondTrack, secondPoint, secondDirection);
  const TVector3 betweenPoints = firstPoint - secondPoint;
  const double cosAngle = firstDirection.Dot(secondDirection);
  const double denominator = 1.0 - cosAngle * cosAngle;
  if (denominator < 1e-8)
    return std::nullopt;
  const double firstDot = firstDirection.Dot(betweenPoints);
  const double secondDot = secondDirection.Dot(betweenPoints);
  const double firstDistance = (cosAngle * secondDot - firstDot) / denominator;
  const double secondDistance = (secondDot - cosAngle * firstDot) / denominator;
  return 0.5 * (firstPoint + firstDistance * firstDirection + secondPoint + secondDistance * secondDirection);
}

/// J^T J and J^T r from a numerical Jacobian at fixed uncertainties
void normalEquations(const std::vector<const ArcTrack*>& tracks, const TVector3& vertex, const Config& config,
                     const Residuals& atVertex, TMatrixDSym& jacobianProduct, TVectorD& gradient) {
  const double stepSize = config.jacobianStep; // mm
  const std::size_t numberOfResiduals = atVertex.normalised.size();
  std::vector<std::array<double, 3>> jacobian(numberOfResiduals);
  Residuals shifted = atVertex;
  for (int coordinate = 0; coordinate < 3; ++coordinate) {
    TVector3 shift(0, 0, 0);
    shift[coordinate] = stepSize;
    computeResiduals(tracks, vertex + shift, config, shifted, false);
    for (std::size_t residual = 0; residual < numberOfResiduals; ++residual)
      jacobian[residual][coordinate] = (shifted.normalised[residual] - atVertex.normalised[residual]) / stepSize;
  }
  jacobianProduct.ResizeTo(3, 3);
  jacobianProduct.Zero();
  gradient.ResizeTo(3);
  gradient.Zero();
  for (std::size_t residual = 0; residual < numberOfResiduals; ++residual)
    for (int row = 0; row < 3; ++row) {
      gradient[row] += jacobian[residual][row] * atVertex.normalised[residual];
      for (int column = 0; column < 3; ++column)
        jacobianProduct(row, column) += jacobian[residual][row] * jacobian[residual][column];
    }
}
} // namespace

std::optional<VertexResult> fitVertex(const std::vector<const ArcTrack*>& tracks, const Config& config) {
  if (tracks.size() < 2)
    return std::nullopt;
  auto startingVertex = tangentLinesClosestApproach(*tracks[0], *tracks[1]);
  if (!startingVertex)
    return std::nullopt;
  TVector3 vertex = *startingVertex;
  Residuals residuals;
  TMatrixDSym jacobianProduct(3);
  TVectorD gradient(3);

  // Gauss-Newton iterations
  for (int iteration = 0; iteration < config.maximumIterations; ++iteration) {
    computeResiduals(tracks, vertex, config, residuals, true);
    normalEquations(tracks, vertex, config, residuals, jacobianProduct, gradient);
    double determinant = 0;
    jacobianProduct.Invert(&determinant);
    if (determinant == 0 || !std::isfinite(determinant))
      return std::nullopt;
    const TVectorD correction = -1.0 * (jacobianProduct * gradient);
    vertex += TVector3(correction[0], correction[1], correction[2]);
    if (!std::isfinite(vertex.X()))
      return std::nullopt;
    if (std::sqrt(correction.Norm2Sqr()) < config.convergenceDistance)
      break;
  }

  // covariance and chi2 at the final vertex
  computeResiduals(tracks, vertex, config, residuals, true);
  normalEquations(tracks, vertex, config, residuals, jacobianProduct, gradient);
  double determinant = 0;
  jacobianProduct.Invert(&determinant);
  if (determinant == 0 || !std::isfinite(determinant))
    return std::nullopt;

  VertexResult result;
  result.position = vertex;
  result.covariance = jacobianProduct;
  for (double residual : residuals.normalised)
    result.chi2 += residual * residual;
  result.degreesOfFreedom = 2 * static_cast<int>(tracks.size()) - 3;
  result.arcPositionOfVertex = residuals.arcPosition;
  result.upstream = std::all_of(residuals.arcPosition.begin(), residuals.arcPosition.end(),
                                [&](double arcPosition) { return arcPosition < config.upstreamTolerance; });
  for (std::size_t first = 0; first < tracks.size(); ++first)
    for (std::size_t second = first + 1; second < tracks.size(); ++second) {
      const double distance = (tracks[first]->point(residuals.arcPosition[first]) -
                               tracks[second]->point(residuals.arcPosition[second]))
                                  .Mag();
      result.pairDistanceOfClosestApproach = std::max(result.pairDistanceOfClosestApproach, distance);
    }
  // the decay must lie between the IP and the first hits; an arc extrapolated far
  // enough can otherwise cross another track on the far side of the IP
  result.betweenIPAndHits = std::all_of(tracks.begin(), tracks.end(), [&](const ArcTrack* track) {
    return vertex.Mag() < track->firstHit().Mag() && vertex.Dot(track->firstHit()) > 0;
  });
  return result;
}

} // namespace k4vertexing::arcpair
