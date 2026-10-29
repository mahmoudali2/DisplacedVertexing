/*
 * Gaudi component: geometry-free 3D two-track vertexing of long-lived-particle
 * decays.
 *
 * Each track is described by a helix taken from its fitted track state and
 * propagated inward in 3D through a two-zone field map, with the energy loss of a
 * material map added back; pairs are vertexed at the 3D point of closest approach.
 * Nothing assumes a barrel geometry, so tracks that cross endcap discs are
 * vertexed like barrel tracks. See include/HelixPairVertexing.h.
 *
 * The detector (field map, material, region) is given by properties; the
 * defaults describe no detector. Example settings are in options/.
 * Cuts on the pair mass and opening angle (MinMass/MaxMass/MinCosOpening) are
 * off by default.
 *
 * Each track is described by its track state at StateLocation (d0, phi, omega,
 * z0, tanLambda); its innermost space-point hit only fixes where the trajectory
 * starts (drift-chamber wire hits are not space points and are skipped). With
 * StateReference = referencePoint the parameters refer to the state's
 * referencePoint (EDM4hep standard); "origin" is for producers that give the
 * parameters w.r.t. (0,0,0) whatever the referencePoint.
 *
 * Properties are grouped by use: track input, detector (field, material),
 * physics, algorithm, vertex selection, output.
 *
 * Author: Mahmoud Althakeel
 *
 * Input : edm4hep::TrackCollection  (tracks with a track state)
 * Output: edm4hep::VertexCollection, best (smallest DCA) pairs first, each track used once
 *   parameters = [ntracks(=2), index1, index2, pairDCA, pathLength1, pathLength2, upstream, inRegion,
 *                 chargeProduct, betweenIPAndHits, mass, cosOpening, momentum1, momentum2,
 *                 closestApproachX, closestApproachY, closestApproachZ]
 *   (indices into the input collection; pathLength = path length of the vertex from each
 *    track's start point, negative = upstream; chargeProduct = product of the xy
 *    rotation senses, -1 = opposite charge, 0 = a straight track; momenta at the vertex)
 *   When both track states carry a covariance the position is the fitted vertex with
 *   its covariance, chi2 and ndf = 1; otherwise it is the closest-approach midpoint,
 *   covMatrix is zero, chi2 = -1 and ndf = 0.
 */
#include "Gaudi/Property.h"
#include "k4FWCore/Transformer.h"

#include <edm4hep/SenseWireHit.h>
#include <edm4hep/TrackCollection.h>
#include <edm4hep/TrackState.h>
#include <edm4hep/VertexCollection.h>

#include "HelixPairVertexing.h"
#include "VertexRegion.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace k4vertexing::helixpair;

struct HelixPairVertexFinder final
    : k4FWCore::Transformer<edm4hep::VertexCollection(const edm4hep::TrackCollection&)> {

  HelixPairVertexFinder(const std::string& name, ISvcLocator* serviceLocator)
      : Transformer(name, serviceLocator, {KeyValues("InputTracks", {"Tracks"})},
                    {KeyValues("OutputVertices", {"HelixPairVertices"})}) {}

  StatusCode initialize() override {
    const auto status = Transformer::initialize();
    if (!status.isSuccess())
      return status;

    // copy the properties into the kernel configuration
    m_config.innerField = m_innerField;
    m_config.outerField = m_outerField;
    m_config.fieldRadius = m_fieldRadius;
    m_config.fieldHalfLength = m_fieldHalfLength;
    m_config.shellInner = m_shellInner;
    m_config.shellOuter = m_shellOuter;
    m_config.shellBeamHole = m_shellBeamHole;
    m_config.energyLossInShell = m_energyLossInShell;
    m_config.energyLossOutside = m_energyLossOutside;
    m_config.particleMass = m_particleMass;
    m_config.stepLength = m_stepLength;
    m_config.backwardFactor = m_backwardFactor;
    m_config.upstreamTolerance = m_upstreamTolerance;

    // checks of the configuration
    if (m_innerField == 0 && m_outerField == 0) {
      error() << "Set the field map: InnerField (and OuterField, FieldR, FieldHalfZ for a two-zone map) [T]"
              << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_energyLossInShell > 0 && m_shellOuter <= m_shellInner) {
      error() << "ShellOuter must be larger than ShellInner when dEdxShell > 0" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_stepLength <= 0) {
      error() << "Step must be positive" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_stateReference.value() != "referencePoint" && m_stateReference.value() != "origin") {
      error() << "StateReference must be referencePoint or origin" << endmsg;
      return StatusCode::FAILURE;
    }
    if (!m_vertexRegion.configure(m_region, m_regionRadius, m_regionHalfLength)) {
      error() << "Region must be one of outside, barrel, endcap, inside, any (got '" << m_region.value() << "')"
              << endmsg;
      return StatusCode::FAILURE;
    }

    info() << "Helix pair vertexer: field " << m_innerField.value() << " T (r < " << m_fieldRadius.value()
           << ", |z| < " << m_fieldHalfLength.value() << " mm) / " << m_outerField.value() << " T, dE/dx shell "
           << m_energyLossInShell.value() << " outside " << m_energyLossOutside.value() << " GeV/mm, step "
           << m_stepLength.value() << " mm, Region = " << m_region.value()
           << ", MaxPairDCA = " << m_maxPairDistanceOfClosestApproach.value()
           << " mm, StateLocation = " << m_stateLocation.value() << endmsg;
    return StatusCode::SUCCESS;
  }

  edm4hep::VertexCollection operator()(const edm4hep::TrackCollection& tracks) const override {
    edm4hep::VertexCollection vertices;

    // ---- helix of every track from its track state ----------------------------
    std::vector<HelixTrack> helixTracks;
    std::vector<int> indexInInputCollection;
    int trackIndex = 0;
    for (const auto& track : tracks) {
      // space-point hits: the innermost one is where the trajectory starts
      std::vector<TVector3> hitPositions;
      for (const auto& hit : track.getTrackerHits()) {
        // a drift-chamber hit stores the wire position, not a point on the track
        if (hit.isA<edm4hep::SenseWireHit>())
          continue;
        hitPositions.emplace_back(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z);
      }
      // find the requested state; tracks without it are skipped
      HelixTrack helixTrack;
      bool isBuilt = false;
      for (const auto& trackState : track.getTrackStates()) {
        if (trackState.location != m_stateLocation)
          continue;
        TrackStateParameters parameters;
        parameters.d0 = trackState.D0;
        parameters.phi = trackState.phi;
        parameters.omega = trackState.omega;
        parameters.z0 = trackState.Z0;
        parameters.tanLambda = trackState.tanLambda;
        // covariance of (d0, phi, omega, z0, tanLambda)
        const edm4hep::TrackParams parameterOrder[5] = {edm4hep::TrackParams::d0, edm4hep::TrackParams::phi,
                                                        edm4hep::TrackParams::omega, edm4hep::TrackParams::z0,
                                                        edm4hep::TrackParams::tanLambda};
        double sumOfVariances = 0;
        for (int row = 0; row < 5; ++row)
          for (int column = 0; column < 5; ++column) {
            parameters.covariance[row][column] = trackState.getCovMatrix(parameterOrder[row], parameterOrder[column]);
            if (row == column)
              sumOfVariances += parameters.covariance[row][column];
          }
        parameters.hasCovariance = m_useCovariance && sumOfVariances > 0;
        if (m_stateReversed) {
          // the state describes the track moving inward: the same helix in the
          // EDM4hep (outward) convention has phi + pi, -omega and -tanLambda
          parameters.phi = trackState.phi + M_PI;
          parameters.omega = -trackState.omega;
          parameters.tanLambda = -trackState.tanLambda;
          // omega and tanLambda change sign: so do their correlations with the others
          const double signChange[5] = {1, 1, -1, 1, -1};
          for (int row = 0; row < 5; ++row)
            for (int column = 0; column < 5; ++column)
              parameters.covariance[row][column] *= signChange[row] * signChange[column];
        }
        if (m_stateReference.value() == "referencePoint")
          parameters.referencePoint = TVector3(trackState.referencePoint[0], trackState.referencePoint[1],
                                               trackState.referencePoint[2]);
        else
          parameters.referencePoint = TVector3(0, 0, 0);
        isBuilt = helixTrack.build(parameters, hitPositions, m_config);
        break;
      }
      if (isBuilt) {
        helixTracks.push_back(std::move(helixTrack));
        indexInInputCollection.push_back(trackIndex);
      }
      ++trackIndex;
    }
    debug() << "Helix built for " << helixTracks.size() << " of " << tracks.size() << " tracks" << endmsg;

    // ---- vertex every pair, apply the selection ----------------------------
    struct VertexCandidate {
      VertexResult result;
      std::size_t firstTrack, secondTrack;
      int chargeProduct;
      bool inRegion;
    };
    std::vector<VertexCandidate> candidates;
    for (std::size_t first = 0; first < helixTracks.size(); ++first)
      for (std::size_t second = first + 1; second < helixTracks.size(); ++second) {
        const int chargeProduct = helixTracks[first].rotationSense() * helixTracks[second].rotationSense();
        // a straight track has no charge information: accepted by the charge cut
        if (m_requireOppositeCharge && chargeProduct > 0)
          continue;
        auto result = fitVertex(helixTracks[first], helixTracks[second], m_config);
        if (!result)
          continue;
        const bool inRegion = m_vertexRegion.contains(result->position);
        if (m_requireUpstream && !result->upstream)
          continue;
        if (m_requireBetweenIPAndHits && !result->betweenIPAndHits)
          continue;
        if (!inRegion || result->pairDistanceOfClosestApproach > m_maxPairDistanceOfClosestApproach)
          continue;
        const bool massOutsideWindow = result->mass < m_minimumMass || result->mass > m_maximumMass;
        if (std::isfinite(result->mass) && massOutsideWindow)
          continue;
        if (result->cosOpening < m_minimumCosOpening)
          continue;
        if (result->isFitted && result->chi2 > m_maximumChi2)
          continue;
        candidates.push_back({std::move(*result), first, second, chargeProduct, inRegion});
      }

    // ---- smallest DCA first, each track used once ---------------------------
    std::sort(candidates.begin(), candidates.end(), [](const VertexCandidate& left, const VertexCandidate& right) {
      return left.result.pairDistanceOfClosestApproach < right.result.pairDistanceOfClosestApproach;
    });
    std::vector<bool> trackIsUsed(helixTracks.size(), false);
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
      if (result.isFitted) {
        // lower triangle: xx, yx, yy, zx, zy, zz
        vertex.setCovMatrix({static_cast<float>(result.covariance[0][0]), static_cast<float>(result.covariance[0][1]),
                             static_cast<float>(result.covariance[1][1]), static_cast<float>(result.covariance[0][2]),
                             static_cast<float>(result.covariance[1][2]), static_cast<float>(result.covariance[2][2])});
        vertex.setChi2(static_cast<float>(result.chi2));
        vertex.setNdf(1);
      } else {
        vertex.setChi2(-1.f);
        vertex.setNdf(0);
      }
      vertex.setPrimary(false);
      vertex.setAlgorithmType(m_algorithmType);
      const float parameters[] = {2.f,
                                  static_cast<float>(indexInInputCollection[candidate.firstTrack]),
                                  static_cast<float>(indexInInputCollection[candidate.secondTrack]),
                                  static_cast<float>(result.pairDistanceOfClosestApproach),
                                  static_cast<float>(result.pathLengthFromStart[0]),
                                  static_cast<float>(result.pathLengthFromStart[1]),
                                  result.upstream ? 1.f : 0.f,
                                  candidate.inRegion ? 1.f : 0.f,
                                  static_cast<float>(candidate.chargeProduct),
                                  result.betweenIPAndHits ? 1.f : 0.f,
                                  static_cast<float>(result.mass),
                                  static_cast<float>(result.cosOpening),
                                  static_cast<float>(result.momentumAtVertex[0]),
                                  static_cast<float>(result.momentumAtVertex[1]),
                                  static_cast<float>(result.closestApproachMidpoint.X()),
                                  static_cast<float>(result.closestApproachMidpoint.Y()),
                                  static_cast<float>(result.closestApproachMidpoint.Z())};
      for (float value : parameters)
        vertex.addToParameters(value);
      if (!m_allowMultipleVertices)
        break;
    }
    debug() << "Reconstructed " << vertices.size() << " helix pair vertices" << endmsg;
    return vertices;
  }

private:
  Config m_config;
  k4vertexing::VertexRegion m_vertexRegion;

  // ===== track input ===========================================================
  Gaudi::Property<int> m_stateLocation{this, "StateLocation", edm4hep::TrackState::AtFirstHit,
                                       "Track state that describes each track (1 AtIP, 2 AtFirstHit, 3 AtLastHit)"};
  Gaudi::Property<std::string> m_stateReference{
      this, "StateReference", "referencePoint",
      "Point the state parameters refer to: referencePoint (EDM4hep) or origin"};
  Gaudi::Property<bool> m_stateReversed{
      this, "StateReversed", false,
      "The state describes the track moving inward: phi, omega and tanLambda are turned to the outward "
      "direction"};
  Gaudi::Property<bool> m_useCovariance{this, "UseCovariance", true,
                                        "Fit the vertex with the state covariances (chi2, vertex covariance)"};

  // ===== detector: magnetic field ==============================================
  Gaudi::Property<double> m_innerField{this, "InnerField", 0.0,
                                       "Bz inside the cylinder (FieldR, FieldHalfZ) [T] (must be set)"};
  Gaudi::Property<double> m_outerField{this, "OuterField", 0.0, "Bz outside the cylinder [T]"};
  Gaudi::Property<double> m_fieldRadius{this, "FieldR", 1.0e12,
                                        "Radius of the inner-field cylinder [mm] (default: InnerField everywhere)"};
  Gaudi::Property<double> m_fieldHalfLength{this, "FieldHalfZ", 1.0e12,
                                            "Half length of the inner-field cylinder [mm]"};

  // ===== detector: material (energy loss added back while propagating) ===========
  Gaudi::Property<double> m_shellInner{this, "ShellInner", 0.0, "Inner box (r and |z|) of the material shell [mm]"};
  Gaudi::Property<double> m_shellOuter{this, "ShellOuter", 0.0, "Outer box (r and |z|) of the material shell [mm]"};
  Gaudi::Property<double> m_shellBeamHole{this, "ShellBeamHole", 0.0, "Beam-hole radius of the shell [mm]"};
  Gaudi::Property<double> m_energyLossInShell{this, "dEdxShell", 0.0,
                                              "Energy loss in the shell [GeV/mm]; 0 = none"};
  Gaudi::Property<double> m_energyLossOutside{this, "dEdxOutside", 0.0,
                                              "Energy loss beyond the shell [GeV/mm]; 0 = none"};

  // ===== physics =================================================================
  Gaudi::Property<double> m_particleMass{this, "ParticleMass", 0.1056584,
                                         "Mass hypothesis of both tracks for the pair mass [GeV]"};

  // ===== algorithm: propagation and vertex search ===================================
  Gaudi::Property<double> m_stepLength{this, "Step", 20.0, "Propagation step [mm]"};
  Gaudi::Property<double> m_backwardFactor{this, "BackFactor", 1.3,
                                           "Propagate inward up to this factor x distance of the first hit"};
  Gaudi::Property<double> m_upstreamTolerance{this, "UpstreamTolerance", 100.0,
                                              "Vertex may lie this far past a first hit [mm]"};

  // ===== vertex selection ==========================================================
  Gaudi::Property<double> m_maxPairDistanceOfClosestApproach{
      this, "MaxPairDCA", 300.0, "Max distance between the two trajectories at the vertex [mm]"};
  Gaudi::Property<bool> m_requireOppositeCharge{this, "RequireOppositeCharge", true,
                                                "Keep only pairs with opposite rotation sense"};
  Gaudi::Property<bool> m_requireUpstream{this, "RequireUpstream", true,
                                          "Vertex must lie before every track's first hit"};
  Gaudi::Property<bool> m_requireBetweenIPAndHits{this, "RequireBetweenIPAndHits", true,
                                                  "Vertex closer to the origin than every first hit, on the same side"};
  Gaudi::Property<std::string> m_region{
      this, "Region", "any",
      "Vertex region w.r.t. the cylinder (RegionR, RegionHalfZ): outside, barrel, endcap, inside, any"};
  Gaudi::Property<double> m_regionRadius{this, "RegionR", 0.0, "Radius of the region cylinder [mm]"};
  Gaudi::Property<double> m_regionHalfLength{this, "RegionHalfZ", 0.0, "Half length of the region cylinder [mm]"};
  Gaudi::Property<double> m_minimumMass{this, "MinMass", 0.0, "Min pair mass at the vertex [GeV]"};
  Gaudi::Property<double> m_maximumMass{this, "MaxMass", 1.0e9, "Max pair mass at the vertex [GeV]"};
  Gaudi::Property<double> m_minimumCosOpening{this, "MinCosOpening", -1.0, "Min cosine of the opening angle"};
  Gaudi::Property<double> m_maximumChi2{this, "MaxChi2", 1.0e9, "Max vertex chi2 (ndf = 1), fitted vertices only"};

  // ===== output ====================================================================
  Gaudi::Property<bool> m_allowMultipleVertices{this, "AllowMultipleVertices", true,
                                                "Write every non-overlapping pair, not only the best"};
  Gaudi::Property<bool> m_keepOverlappingPairs{this, "KeepOverlappingPairs", false,
                                               "Write pairs sharing a track too (for cut studies)"};
  Gaudi::Property<int> m_algorithmType{this, "AlgorithmType", 4, "edm4hep::Vertex algorithmType tag"};
};

DECLARE_COMPONENT(HelixPairVertexFinder)
