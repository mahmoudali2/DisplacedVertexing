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
 *   parameters = [ntracks(=2), index1, index2, pairDCA, s1, s2, upstream, inRegion,
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
    m_cfg.bFieldForMomentum = m_bFieldForMomentum;
    m_cfg.defaultPt = m_defaultPt;
    m_cfg.momentumResolution = m_momentumResolution;
    m_cfg.minMomentum = m_minMomentum;
    m_cfg.hitResolution = m_hitResolution;
    m_cfg.fieldUncertainty = m_fieldUncertainty;
    m_cfg.useMaterial = m_useMaterial;
    m_cfg.shellX0 = m_shellX0;
    m_cfg.shellInnerR = m_shellInnerR;
    m_cfg.shellOuterR = m_shellOuterR;
    m_cfg.shellInnerZ = m_shellInnerZ;
    m_cfg.shellOuterZ = m_shellOuterZ;
    m_cfg.outsideX0 = m_outsideX0;
    m_cfg.outsideThickness = m_outsideThickness;
    m_cfg.fitCurvature = m_fitCurvature;
    m_cfg.minHits = m_minHits;
    m_cfg.upstreamTolerance = m_upstreamTolerance;
    m_cfg.nSteps = m_materialSteps;
    m_cfg.maxIterations = m_maxIterations;
    m_cfg.jacobianStep = m_jacobianStep;
    m_cfg.convergenceDistance = m_convergenceDistance;

    // checks of the configuration
    if (m_bFieldForMomentum == 0) {
      error() << "BFieldForMomentum must be set: the field [T] the input tracks' omega was computed with" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_shellX0 > 0 && m_shellOuterR <= m_shellInnerR) {
      error() << "ShellOuterR must be larger than ShellInnerR when ShellX0 > 0" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_outsideX0 > 0 && m_outsideThickness <= 0) {
      error() << "OutsideThickness must be positive when OutsideX0 > 0" << endmsg;
      return StatusCode::FAILURE;
    }
    if (!m_regionCut.configure(m_region, m_regionR, m_regionHalfZ)) {
      error() << "Region must be one of outside, barrel, endcap, inside, any (got '" << m_region.value() << "')"
              << endmsg;
      return StatusCode::FAILURE;
    }

    info() << "Arc pair vertexer: B for momentum " << m_bFieldForMomentum.value() << " T, material "
           << (m_useMaterial ? "ON" : "OFF") << " (shell " << m_shellX0.value() << " X0, outside "
           << m_outsideX0.value() << " X0), field uncertainty " << m_fieldUncertainty.value()
           << " T, sigma(p)/p = " << m_momentumResolution.value() << ", hit resolution " << m_hitResolution.value()
           << " mm, Region = " << m_region.value() << ", MaxPairDCA = " << m_maxPairDCA.value()
           << " mm, MaxChi2 = " << m_maxChi2.value() << endmsg;
    return StatusCode::SUCCESS;
  }

  edm4hep::VertexCollection operator()(const edm4hep::TrackCollection& tracks) const override {
    edm4hep::VertexCollection vertices;

    // ---- local arc fit of every track from its hits -------------------------
    std::vector<ArcTrack> arcs;
    std::vector<int> sourceIndex;
    std::vector<int> charge;
    int index = 0;
    for (const auto& track : tracks) {
      std::vector<TVector3> hits;
      for (const auto& hit : track.getTrackerHits()) {
        // a drift-chamber hit stores the wire position, not a point on the track
        if (m_planarHitsOnly && hit.isA<edm4hep::SenseWireHit>())
          continue;
        hits.emplace_back(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z);
      }
      double omega = 0;
      for (const auto& state : track.getTrackStates())
        if (state.location == m_momentumStateLocation) {
          omega = state.omega;
          break;
        }
      ArcTrack arc;
      if (arc.build(hits, omega, m_cfg)) {
        arcs.push_back(std::move(arc));
        sourceIndex.push_back(index);
        // only the product of the two charges is used, so the sign convention
        // of the input state does not matter
        charge.push_back(omega > 0 ? 1 : (omega < 0 ? -1 : 0));
      }
      ++index;
    }
    debug() << "Arc-fitted " << arcs.size() << " of " << tracks.size() << " tracks" << endmsg;

    // ---- vertex every pair, apply the selection ----------------------------
    struct Candidate {
      VertexResult fit;
      std::size_t a, b;
      bool region;
    };
    std::vector<Candidate> candidates;
    for (std::size_t i = 0; i < arcs.size(); ++i)
      for (std::size_t j = i + 1; j < arcs.size(); ++j) {
        auto fit = fitVertex({&arcs[i], &arcs[j]}, m_cfg);
        if (!fit)
          continue;
        const bool region = m_regionCut.contains(fit->position);
        if (!region)
          continue;
        if (m_requireUpstream && !fit->upstream)
          continue;
        if (m_requireBetween && !fit->betweenIPAndHits)
          continue;
        if (fit->pairDCA > m_maxPairDCA || fit->chi2 > m_maxChi2)
          continue;
        // vertex must lie at least MinUpstreamDistance before BOTH first hits:
        // rejects one particle split into two nearby tracks, which "meet" at the hits
        if (*std::max_element(fit->sFromFirst.begin(), fit->sFromFirst.end()) > -m_minUpstreamDistance)
          continue;
        if (m_requireOppositeCharge && charge[i] * charge[j] >= 0)
          continue;
        candidates.push_back({std::move(*fit), i, j, region});
      }

    // ---- best chi2 first, each track used once ------------------------------
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate& x, const Candidate& y) { return x.fit.chi2 < y.fit.chi2; });
    std::vector<bool> used(arcs.size(), false);
    for (const auto& cand : candidates) {
      if (!m_keepOverlapping && (used[cand.a] || used[cand.b]))
        continue;
      used[cand.a] = used[cand.b] = true;
      const auto& f = cand.fit;
      auto vertex = vertices.create();
      vertex.setPosition({static_cast<float>(f.position.X()), static_cast<float>(f.position.Y()),
                          static_cast<float>(f.position.Z())});
      const TMatrixDSym& c = f.covariance;
      vertex.setCovMatrix({static_cast<float>(c(0, 0)), static_cast<float>(c(0, 1)), static_cast<float>(c(1, 1)),
                           static_cast<float>(c(0, 2)), static_cast<float>(c(1, 2)), static_cast<float>(c(2, 2))});
      vertex.setChi2(static_cast<float>(f.chi2));
      vertex.setNdf(f.ndf);
      vertex.setPrimary(false);
      vertex.setAlgorithmType(m_algorithmType);
      for (float value : {2.f, static_cast<float>(sourceIndex[cand.a]), static_cast<float>(sourceIndex[cand.b]),
                          static_cast<float>(f.pairDCA), static_cast<float>(f.sFromFirst[0]),
                          static_cast<float>(f.sFromFirst[1]), f.upstream ? 1.f : 0.f, cand.region ? 1.f : 0.f,
                          static_cast<float>(charge[cand.a] * charge[cand.b]), f.betweenIPAndHits ? 1.f : 0.f})
        vertex.addToParameters(value);
      if (!m_allowMultipleVertices)
        break;
    }
    debug() << "Reconstructed " << vertices.size() << " arc pair vertices" << endmsg;
    return vertices;
  }

private:
  Config m_cfg;
  k4vertexing::VertexRegion m_regionCut;

  // --- momentum for the scattering model ---
  // ===== track input =====================================================================
  Gaudi::Property<int> m_momentumStateLocation{this, "MomentumStateLocation", edm4hep::TrackState::AtFirstHit,
                                               "Track state whose omega gives the momentum"};
  Gaudi::Property<unsigned> m_minHits{this, "MinHits", 3, "Minimum hits per track"};
  Gaudi::Property<bool> m_planarHitsOnly{this, "PlanarHitsOnly", true,
                                         "Use only space-point hits: skip drift-chamber wire hits (SenseWireHit)"};

  // ===== detector: magnetic field ========================================================
  Gaudi::Property<double> m_bFieldForMomentum{
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
  Gaudi::Property<double> m_shellX0{this, "ShellX0", 0.0,
                                    "Radiation lengths of the material shell (radial depth) in front of the hits"};
  Gaudi::Property<double> m_shellInnerR{this, "ShellInnerR", 0.0, "Inner box radius of the shell [mm]"};
  Gaudi::Property<double> m_shellOuterR{this, "ShellOuterR", 0.0, "Outer box radius of the shell [mm]"};
  Gaudi::Property<double> m_shellInnerZ{this, "ShellInnerZ", 0.0, "Inner box half length of the shell [mm]"};
  Gaudi::Property<double> m_shellOuterZ{this, "ShellOuterZ", 0.0, "Outer box half length of the shell [mm]"};
  Gaudi::Property<double> m_outsideX0{this, "OutsideX0", 0.0,
                                      "Radiation lengths beyond the shell (spread over OutsideThickness)"};
  Gaudi::Property<double> m_outsideThickness{this, "OutsideThickness", 1.0,
                                             "Depth over which OutsideX0 is spread [mm]"};

  // ===== physics =========================================================================
  Gaudi::Property<double> m_defaultPt{this, "DefaultPt", 10.0, "pT [GeV] used when the state has omega = 0"};
  Gaudi::Property<double> m_minMomentum{this, "MinMomentum", 2.0, "Momentum floor for scattering [GeV]"};

  // ===== algorithm: track model and vertex fit ===========================================
  Gaudi::Property<bool> m_fitCurvature{this, "FitCurvature", true, "Circle in xy (true) or straight line (false)"};
  Gaudi::Property<double> m_upstreamTolerance{this, "UpstreamTolerance", 100.0, "Vertex may lie this far past a first hit [mm]"};
  Gaudi::Property<int> m_materialSteps{this, "MaterialSteps", 40, "Steps of the material integrals"};
  Gaudi::Property<int> m_maxIterations{this, "MaxIterations", 12, "Gauss-Newton iterations of the vertex fit"};
  Gaudi::Property<double> m_jacobianStep{this, "JacobianStep", 1.0,
                                         "Step of the numerical Jacobian [mm] (scale to the vertex resolution)"};
  Gaudi::Property<double> m_convergenceDistance{this, "ConvergenceDistance", 0.1,
                                                "Stop when the vertex moves less than this [mm]"};

  // ===== vertex selection ================================================================
  Gaudi::Property<double> m_maxPairDCA{this, "MaxPairDCA", 100.0,
                                       "Max distance between the two trajectories at the vertex [mm]"};
  Gaudi::Property<double> m_maxChi2{this, "MaxChi2", 1.0e9, "Max vertex chi2 (ndf = 1)"};
  Gaudi::Property<bool> m_requireOppositeCharge{this, "RequireOppositeCharge", true,
                                                "Keep only opposite-charge pairs"};
  Gaudi::Property<bool> m_requireUpstream{this, "RequireUpstream", true,
                                          "Vertex must lie before every track's first hit"};
  Gaudi::Property<double> m_minUpstreamDistance{this, "MinUpstreamDistance", 0.0,
                                                "Vertex at least this far before both first hits [mm]"};
  Gaudi::Property<bool> m_requireBetween{this, "RequireBetweenIPAndHits", true,
                                         "Vertex closer to the origin than every first hit, on the same side"};
  Gaudi::Property<std::string> m_region{
      this, "Region", "any",
      "Vertex region w.r.t. the cylinder (RegionR, RegionHalfZ): outside, barrel, endcap, inside, any"};
  Gaudi::Property<double> m_regionR{this, "RegionR", 0.0, "Radius of the region cylinder [mm]"};
  Gaudi::Property<double> m_regionHalfZ{this, "RegionHalfZ", 0.0, "Half length of the region cylinder [mm]"};

  // ===== output ==========================================================================
  Gaudi::Property<bool> m_allowMultipleVertices{this, "AllowMultipleVertices", true,
                                                "Write every non-overlapping pair, not only the best"};
  Gaudi::Property<bool> m_keepOverlapping{this, "KeepOverlappingPairs", false,
                                          "Write pairs sharing a track too (for cut studies)"};
  Gaudi::Property<int> m_algorithmType{this, "AlgorithmType", 3, "edm4hep::Vertex algorithmType tag"};
};

DECLARE_COMPONENT(ArcPairVertexFinder)
