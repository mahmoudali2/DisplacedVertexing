/*
 * Gaudi component: two-track vertexing of long-lived-particle decays with local
 * arc fits.
 *
 * Each track is refitted locally from its own hits (circle in xy, line in z,
 * with multiple scattering between the hits) and extrapolated back only to the
 * decay, with the scattering in the material in front of the first hit added
 * to its error. Track pairs are vertexed with a chi2 fit. See
 * include/ArcPairVertexing.h.
 *
 * The detector (field used for the momentum, material, region) is given by
 * properties; the defaults describe no detector. Example settings are in
 * options/.
 *
 * Author: Mahmoud Althakeel
 *
 * Input : edm4hep::TrackCollection  (tracks with their hits attached)
 * Output: edm4hep::VertexCollection
 *   parameters = [ntracks(=2), index1, index2, pairDCA, arcPosition1, arcPosition2, upstream, inRegion,
 *                 chargeProduct, betweenIPAndHits]
 *   (indices into the input collection; s = distance of the vertex from each
 *    track's first hit along the track, negative = upstream)
 */
#include "Gaudi/Property.h"
#include "k4FWCore/Transformer.h"

#include <edm4hep/SenseWireHit.h>
#include <edm4hep/TrackCollection.h>
#include <edm4hep/TrackState.h>
#include <edm4hep/VertexCollection.h>

#include "ArcPairVertexing.h"
#include "VertexRegion.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

using namespace k4vertexing::arcpair;

struct ArcPairVertexFinder final
    : k4FWCore::Transformer<edm4hep::VertexCollection(const edm4hep::TrackCollection&)> {

  ArcPairVertexFinder(const std::string& name, ISvcLocator* serviceLocator)
      : Transformer(name, serviceLocator, {KeyValues("InputTracks", {"Tracks"})},
                    {KeyValues("OutputVertices", {"ArcPairVertices"})}) {}

  StatusCode initialize() override {
    const auto status = Transformer::initialize();
    if (!status.isSuccess())
      return status;

    // copy the properties into the kernel configuration
    m_config.fieldForMomentum = m_fieldForMomentum;
    m_config.defaultTransverseMomentum = m_defaultTransverseMomentum;
    m_config.momentumResolution = m_momentumResolution;
    m_config.minimumMomentum = m_minimumMomentum;
    m_config.hitResolution = m_hitResolution;
    m_config.fieldUncertainty = m_fieldUncertainty;
    m_config.useMaterial = m_useMaterial;
    m_config.shellRadiationLengths = m_shellRadiationLengths;
    m_config.shellInnerRadius = m_shellInnerRadius;
    m_config.shellOuterRadius = m_shellOuterRadius;
    m_config.shellInnerHalfLength = m_shellInnerHalfLength;
    m_config.shellOuterHalfLength = m_shellOuterHalfLength;
    m_config.outsideRadiationLengths = m_outsideRadiationLengths;
    m_config.outsideThickness = m_outsideThickness;
    m_config.fitCurvature = m_fitCurvature;
    m_config.minimumHits = m_minimumHits;
    m_config.upstreamTolerance = m_upstreamTolerance;
    m_config.integrationSteps = m_integrationSteps;
    m_config.maximumIterations = m_maximumIterations;
    m_config.jacobianStep = m_jacobianStep;
    m_config.convergenceDistance = m_convergenceDistance;

    // checks of the configuration
    if (m_fieldForMomentum == 0) {
      error() << "BFieldForMomentum must be set: the field [T] the input tracks' omega was computed with" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_shellRadiationLengths > 0 && m_shellOuterRadius <= m_shellInnerRadius) {
      error() << "ShellOuterR must be larger than ShellInnerR when ShellX0 > 0" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_outsideRadiationLengths > 0 && m_outsideThickness <= 0) {
      error() << "OutsideThickness must be positive when OutsideX0 > 0" << endmsg;
      return StatusCode::FAILURE;
    }
    if (!m_vertexRegion.configure(m_region, m_regionRadius, m_regionHalfLength)) {
      error() << "Region must be one of outside, barrel, endcap, inside, any (got '" << m_region.value() << "')"
              << endmsg;
      return StatusCode::FAILURE;
    }

    info() << "Arc pair vertexer: B for momentum " << m_fieldForMomentum.value() << " T, material "
           << (m_useMaterial ? "ON" : "OFF") << " (shell " << m_shellRadiationLengths.value() << " X0, outside "
           << m_outsideRadiationLengths.value() << " X0), field uncertainty " << m_fieldUncertainty.value()
           << " T, sigma(p)/p = " << m_momentumResolution.value() << ", hit resolution " << m_hitResolution.value()
           << " mm, Region = " << m_region.value() << ", MaxPairDCA = " << m_maxPairDistanceOfClosestApproach.value()
           << " mm, MaxChi2 = " << m_maximumChi2.value() << endmsg;
    return StatusCode::SUCCESS;
  }

  edm4hep::VertexCollection operator()(const edm4hep::TrackCollection& tracks) const override {
    edm4hep::VertexCollection vertices;

    // ---- local arc fit of every track from its hits -------------------------
    std::vector<ArcTrack> arcTracks;
    std::vector<int> indexInInputCollection;
    std::vector<int> chargeSign;
    int trackIndex = 0;
    for (const auto& track : tracks) {
      std::vector<TVector3> hitPositions;
      for (const auto& hit : track.getTrackerHits()) {
        // a drift-chamber hit stores the wire position, not a point on the track
        if (m_planarHitsOnly && hit.isA<edm4hep::SenseWireHit>())
          continue;
        hitPositions.emplace_back(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z);
      }
      // the curvature of the chosen state gives the momentum for the scattering model
      double omega = 0;
      for (const auto& trackState : track.getTrackStates())
        if (trackState.location == m_momentumStateLocation) {
          omega = trackState.omega;
          break;
        }
      ArcTrack arcTrack;
      if (arcTrack.build(hitPositions, omega, m_config)) {
        arcTracks.push_back(std::move(arcTrack));
        indexInInputCollection.push_back(trackIndex);
        // only the product of the two charges is used, so the sign convention
        // of the input state does not matter
        chargeSign.push_back(omega > 0 ? 1 : (omega < 0 ? -1 : 0));
      }
      ++trackIndex;
    }
    debug() << "Arc-fitted " << arcTracks.size() << " of " << tracks.size() << " tracks" << endmsg;

    // ---- vertex every pair, apply the selection ----------------------------
    struct VertexCandidate {
      VertexResult result;
      std::size_t firstTrack, secondTrack;
      bool inRegion;
    };
    std::vector<VertexCandidate> candidates;
    for (std::size_t first = 0; first < arcTracks.size(); ++first)
      for (std::size_t second = first + 1; second < arcTracks.size(); ++second) {
        auto result = fitVertex({&arcTracks[first], &arcTracks[second]}, m_config);
        if (!result)
          continue;
        const bool inRegion = m_vertexRegion.contains(result->position);
        if (!inRegion)
          continue;
        if (m_requireUpstream && !result->upstream)
          continue;
        if (m_requireBetweenIPAndHits && !result->betweenIPAndHits)
          continue;
        if (result->pairDistanceOfClosestApproach > m_maxPairDistanceOfClosestApproach || result->chi2 > m_maximumChi2)
          continue;
        // vertex must lie at least MinUpstreamDistance before BOTH first hits:
        // rejects one particle split into two nearby tracks, which "meet" at the hits
        const double largestArcPosition =
            *std::max_element(result->arcPositionOfVertex.begin(), result->arcPositionOfVertex.end());
        if (largestArcPosition > -m_minimumUpstreamDistance)
          continue;
        if (m_requireOppositeCharge && chargeSign[first] * chargeSign[second] >= 0)
          continue;
        candidates.push_back({std::move(*result), first, second, inRegion});
      }

    // ---- best chi2 first, each track used once ------------------------------
    std::sort(candidates.begin(), candidates.end(), [](const VertexCandidate& left, const VertexCandidate& right) {
      return left.result.chi2 < right.result.chi2;
    });
    std::vector<bool> trackIsUsed(arcTracks.size(), false);
    for (const auto& candidate : candidates) {
      const bool sharesTrack = trackIsUsed[candidate.firstTrack] || trackIsUsed[candidate.secondTrack];
      if (!m_keepOverlappingPairs && sharesTrack)
        continue;
      trackIsUsed[candidate.firstTrack] = true;
      trackIsUsed[candidate.secondTrack] = true;
      const VertexResult& result = candidate.result;
      auto vertex = vertices.create();
      vertex.setPosition({static_cast<float>(result.position.X()), static_cast<float>(result.position.Y()),
                          static_cast<float>(result.position.Z())});
      // lower triangle: xx, yx, yy, zx, zy, zz
      const TMatrixDSym& covariance = result.covariance;
      vertex.setCovMatrix({static_cast<float>(covariance(0, 0)), static_cast<float>(covariance(0, 1)),
                           static_cast<float>(covariance(1, 1)), static_cast<float>(covariance(0, 2)),
                           static_cast<float>(covariance(1, 2)), static_cast<float>(covariance(2, 2))});
      vertex.setChi2(static_cast<float>(result.chi2));
      vertex.setNdf(result.degreesOfFreedom);
      vertex.setPrimary(false);
      vertex.setAlgorithmType(m_algorithmType);
      const float parameters[] = {2.f,
                                  static_cast<float>(indexInInputCollection[candidate.firstTrack]),
                                  static_cast<float>(indexInInputCollection[candidate.secondTrack]),
                                  static_cast<float>(result.pairDistanceOfClosestApproach),
                                  static_cast<float>(result.arcPositionOfVertex[0]),
                                  static_cast<float>(result.arcPositionOfVertex[1]),
                                  result.upstream ? 1.f : 0.f,
                                  candidate.inRegion ? 1.f : 0.f,
                                  static_cast<float>(chargeSign[candidate.firstTrack] * chargeSign[candidate.secondTrack]),
                                  result.betweenIPAndHits ? 1.f : 0.f};
      for (float value : parameters)
        vertex.addToParameters(value);
      if (!m_allowMultipleVertices)
        break;
    }
    debug() << "Reconstructed " << vertices.size() << " arc pair vertices" << endmsg;
    return vertices;
  }

private:
  Config m_config;
  k4vertexing::VertexRegion m_vertexRegion;

  // --- momentum for the scattering model ---
  // ===== track input =====================================================================
  Gaudi::Property<int> m_momentumStateLocation{this, "MomentumStateLocation", edm4hep::TrackState::AtFirstHit,
                                               "Track state whose omega gives the momentum"};
  Gaudi::Property<unsigned> m_minimumHits{this, "MinHits", 3, "Minimum hits per track"};
  Gaudi::Property<bool> m_planarHitsOnly{this, "PlanarHitsOnly", true,
                                         "Use only space-point hits: skip drift-chamber wire hits (SenseWireHit)"};

  // ===== detector: magnetic field ========================================================
  Gaudi::Property<double> m_fieldForMomentum{
      this, "BFieldForMomentum", 0.0,
      "Field [T] the input tracking used to compute that state's omega (must be set)"};
  Gaudi::Property<double> m_fieldUncertainty{
      this, "FieldUncertainty", 0.0,
      "Field mismatch [T] between where the arc is fitted and where it is extrapolated "
      "(bending-plane error 0.5*0.3*dB/pT*s^2); 0 = off"};

  // ===== detector: resolutions ===========================================================
  Gaudi::Property<double> m_hitResolution{this, "HitResolution", 1.0, "Hit resolution per coordinate [mm]"};
  Gaudi::Property<double> m_momentumResolution{this, "MomentumResolution", 0.0, "sigma(p)/p of the input tracks"};

  // ===== detector: material (multiple scattering in the track errors) ====================
  Gaudi::Property<bool> m_useMaterial{this, "UseMaterialEffects", true, "Multiple scattering in the track errors"};
  Gaudi::Property<double> m_shellRadiationLengths{this, "ShellX0", 0.0,
                                    "Radiation lengths of the material shell (radial depth) in front of the hits"};
  Gaudi::Property<double> m_shellInnerRadius{this, "ShellInnerR", 0.0, "Inner box radius of the shell [mm]"};
  Gaudi::Property<double> m_shellOuterRadius{this, "ShellOuterR", 0.0, "Outer box radius of the shell [mm]"};
  Gaudi::Property<double> m_shellInnerHalfLength{this, "ShellInnerZ", 0.0, "Inner box half length of the shell [mm]"};
  Gaudi::Property<double> m_shellOuterHalfLength{this, "ShellOuterZ", 0.0, "Outer box half length of the shell [mm]"};
  Gaudi::Property<double> m_outsideRadiationLengths{this, "OutsideX0", 0.0,
                                      "Radiation lengths beyond the shell (spread over OutsideThickness)"};
  Gaudi::Property<double> m_outsideThickness{this, "OutsideThickness", 1.0,
                                             "Depth over which OutsideX0 is spread [mm]"};

  // ===== physics =========================================================================
  Gaudi::Property<double> m_defaultTransverseMomentum{this, "DefaultPt", 10.0, "pT [GeV] used when the state has omega = 0"};
  Gaudi::Property<double> m_minimumMomentum{this, "MinMomentum", 2.0, "Momentum floor for scattering [GeV]"};

  // ===== algorithm: track model and vertex fit ===========================================
  Gaudi::Property<bool> m_fitCurvature{this, "FitCurvature", true, "Circle in xy (true) or straight line (false)"};
  Gaudi::Property<double> m_upstreamTolerance{this, "UpstreamTolerance", 100.0, "Vertex may lie this far past a first hit [mm]"};
  Gaudi::Property<int> m_integrationSteps{this, "MaterialSteps", 40, "Steps of the material integrals"};
  Gaudi::Property<int> m_maximumIterations{this, "MaxIterations", 12, "Gauss-Newton iterations of the vertex fit"};
  Gaudi::Property<double> m_jacobianStep{this, "JacobianStep", 1.0,
                                         "Step of the numerical Jacobian [mm] (scale to the vertex resolution)"};
  Gaudi::Property<double> m_convergenceDistance{this, "ConvergenceDistance", 0.1,
                                                "Stop when the vertex moves less than this [mm]"};

  // ===== vertex selection ================================================================
  Gaudi::Property<double> m_maxPairDistanceOfClosestApproach{this, "MaxPairDCA", 100.0,
                                       "Max distance between the two trajectories at the vertex [mm]"};
  Gaudi::Property<double> m_maximumChi2{this, "MaxChi2", 1.0e9, "Max vertex chi2 (ndf = 1)"};
  Gaudi::Property<bool> m_requireOppositeCharge{this, "RequireOppositeCharge", true,
                                                "Keep only opposite-charge pairs"};
  Gaudi::Property<bool> m_requireUpstream{this, "RequireUpstream", true,
                                          "Vertex must lie before every track's first hit"};
  Gaudi::Property<double> m_minimumUpstreamDistance{this, "MinUpstreamDistance", 0.0,
                                                "Vertex at least this far before both first hits [mm]"};
  Gaudi::Property<bool> m_requireBetweenIPAndHits{this, "RequireBetweenIPAndHits", true,
                                         "Vertex closer to the origin than every first hit, on the same side"};
  Gaudi::Property<std::string> m_region{
      this, "Region", "any",
      "Vertex region w.r.t. the cylinder (RegionR, RegionHalfZ): outside, barrel, endcap, inside, any"};
  Gaudi::Property<double> m_regionRadius{this, "RegionR", 0.0, "Radius of the region cylinder [mm]"};
  Gaudi::Property<double> m_regionHalfLength{this, "RegionHalfZ", 0.0, "Half length of the region cylinder [mm]"};

  // ===== output ==========================================================================
  Gaudi::Property<bool> m_allowMultipleVertices{this, "AllowMultipleVertices", true,
                                                "Write every non-overlapping pair, not only the best"};
  Gaudi::Property<bool> m_keepOverlappingPairs{this, "KeepOverlappingPairs", false,
                                          "Write pairs sharing a track too (for cut studies)"};
  Gaudi::Property<int> m_algorithmType{this, "AlgorithmType", 3, "edm4hep::Vertex algorithmType tag"};
};

DECLARE_COMPONENT(ArcPairVertexFinder)
