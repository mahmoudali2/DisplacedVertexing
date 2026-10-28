/*
 * Gaudi component: geometry-free 3D two-track vertexing of long-lived-particle
 * decays.
 *
 * Companion of ArcPairVertexFinder (same input, same output layout plus the pair
 * kinematics). Each track is a helix fitted to its own hits and propagated
 * inward in 3D through a two-zone field map with the energy loss of a material
 * map added back; pairs are vertexed at the 3D point of closest approach.
 * Nothing assumes a barrel geometry, so tracks that cross endcap discs are
 * vertexed like barrel tracks. See include/HelixPairVertexing.h.
 *
 * The detector (field map, material, region) is given by properties; the
 * defaults describe no detector. Example settings are in options/.
 * Cuts on the pair mass and opening angle (MinMass/MaxMass/MinCosOpening) are
 * off by default.
 *
 * TrackInput selects what describes each track:
 *   hits  : helix fitted to the track's own hits (default)
 *   state : helix from the track state at StateLocation (d0, phi, omega, z0,
 *           tanLambda); the hits only fix where the trajectory starts. With
 *           StateReference = referencePoint the parameters refer to the state's
 *           referencePoint (EDM4hep standard); "origin" is for producers that
 *           give the parameters w.r.t. (0,0,0) whatever the referencePoint.
 *
 * Author: Mahmoud Althakeel
 *
 * Input : edm4hep::TrackCollection  (tracks with their hits attached)
 * Output: edm4hep::VertexCollection, best (smallest DCA) pairs first, each track used once
 *   parameters = [ntracks(=2), index1, index2, pairDCA, t1, t2, upstream, inRegion,
 *                 chargeProduct, betweenIPAndHits, mass, cosOpening, p1, p2]
 *   (indices into the input collection; t = path length of the vertex from each
 *    track's first hit, negative = upstream; chargeProduct = product of the xy
 *    rotation senses, -1 = opposite charge, 0 = a straight track; p at the vertex)
 *   No covariance is computed: covMatrix is zero, chi2 = pairDCA [mm], ndf = 0.
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
    m_cfg.innerField = m_innerField;
    m_cfg.outerField = m_outerField;
    m_cfg.fieldR = m_fieldR;
    m_cfg.fieldHalfZ = m_fieldHalfZ;
    m_cfg.shellInner = m_shellInner;
    m_cfg.shellOuter = m_shellOuter;
    m_cfg.shellHole = m_shellHole;
    m_cfg.dEdxShell = m_dEdxShell;
    m_cfg.dEdxOutside = m_dEdxOutside;
    m_cfg.particleMass = m_particleMass;
    m_cfg.minHits = m_minHits;
    m_cfg.maxHitsInFit = m_maxHitsInFit;
    m_cfg.minXYSpan = m_minXYSpan;
    m_cfg.maxRadius = m_maxRadius;
    m_cfg.step = m_step;
    m_cfg.backFactor = m_backFactor;
    m_cfg.upstreamTolerance = m_upstreamTolerance;

    // checks of the configuration
    if (m_innerField == 0 && m_outerField == 0) {
      error() << "Set the field map: InnerField (and OuterField, FieldR, FieldHalfZ for a two-zone map) [T]"
              << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_dEdxShell > 0 && m_shellOuter <= m_shellInner) {
      error() << "ShellOuter must be larger than ShellInner when dEdxShell > 0" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_step <= 0) {
      error() << "Step must be positive" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_trackInput.value() != "hits" && m_trackInput.value() != "state") {
      error() << "TrackInput must be hits or state (got '" << m_trackInput.value() << "')" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_stateReference.value() != "referencePoint" && m_stateReference.value() != "origin") {
      error() << "StateReference must be referencePoint or origin" << endmsg;
      return StatusCode::FAILURE;
    }
    if (!m_regionCut.configure(m_region, m_regionR, m_regionHalfZ)) {
      error() << "Region must be one of outside, barrel, endcap, inside, any (got '" << m_region.value() << "')"
              << endmsg;
      return StatusCode::FAILURE;
    }

    info() << "Helix pair vertexer: field " << m_innerField.value() << " T (r < " << m_fieldR.value()
           << ", |z| < " << m_fieldHalfZ.value() << " mm) / " << m_outerField.value() << " T, dE/dx shell "
           << m_dEdxShell.value() << " outside " << m_dEdxOutside.value() << " GeV/mm, step " << m_step.value()
           << " mm, Region = " << m_region.value() << ", MaxPairDCA = " << m_maxPairDCA.value()
           << " mm, TrackInput = " << m_trackInput.value() << endmsg;
    return StatusCode::SUCCESS;
  }

  edm4hep::VertexCollection operator()(const edm4hep::TrackCollection& tracks) const override {
    edm4hep::VertexCollection vertices;

    // ---- helix fit of every track from its hits -----------------------------
    std::vector<HelixTrack> helices;
    std::vector<int> sourceIndex;
    int index = 0;
    for (const auto& track : tracks) {
      std::vector<TVector3> hits;
      for (const auto& hit : track.getTrackerHits()) {
        // a drift-chamber hit stores the wire position, not a point on the track
        if (m_planarHitsOnly && hit.isA<edm4hep::SenseWireHit>())
          continue;
        hits.emplace_back(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z);
      }
      HelixTrack helix;
      bool ok = false;
      if (m_trackInput.value() == "hits") {
        ok = helix.build(hits, m_cfg);
      } else {
        // find the requested state; tracks without it are skipped
        for (const auto& state : track.getTrackStates()) {
          if (state.location != m_stateLocation)
            continue;
          StateParams params;
          params.d0 = state.D0;
          params.phi = state.phi;
          params.omega = state.omega;
          params.z0 = state.Z0;
          params.tanLambda = state.tanLambda;
          if (m_stateReversed) {
            // the state describes the track moving inward: the same helix in the
            // EDM4hep (outward) convention has phi + pi, -omega and -tanLambda
            params.phi = state.phi + M_PI;
            params.omega = -state.omega;
            params.tanLambda = -state.tanLambda;
          }
          if (m_stateReference.value() == "referencePoint")
            params.reference = TVector3(state.referencePoint[0], state.referencePoint[1], state.referencePoint[2]);
          else
            params.reference = TVector3(0, 0, 0);
          ok = helix.buildFromState(params, hits, m_cfg);
          break;
        }
      }
      if (ok) {
        helices.push_back(std::move(helix));
        sourceIndex.push_back(index);
      }
      ++index;
    }
    debug() << "Helix-fitted " << helices.size() << " of " << tracks.size() << " tracks" << endmsg;

    // ---- vertex every pair, apply the selection ----------------------------
    struct Candidate {
      VertexResult fit;
      std::size_t a, b;
      int chargeProduct;
      bool region;
    };
    std::vector<Candidate> candidates;
    for (std::size_t i = 0; i < helices.size(); ++i)
      for (std::size_t j = i + 1; j < helices.size(); ++j) {
        const int q = helices[i].rotation() * helices[j].rotation();
        // a straight track has no charge information: accepted by the charge cut
        if (m_requireOppositeCharge && q > 0)
          continue;
        auto fit = fitVertex(helices[i], helices[j], m_cfg);
        if (!fit)
          continue;
        const bool region = m_regionCut.contains(fit->position);
        if (m_requireUpstream && !fit->upstream)
          continue;
        if (m_requireBetween && !fit->betweenIPAndHits)
          continue;
        if (!region || fit->pairDCA > m_maxPairDCA)
          continue;
        if (std::isfinite(fit->mass) && (fit->mass < m_minMass || fit->mass > m_maxMass))
          continue;
        if (fit->cosOpening < m_minCosOpening)
          continue;
        candidates.push_back({std::move(*fit), i, j, q, region});
      }

    // ---- smallest DCA first, each track used once ---------------------------
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate& x, const Candidate& y) { return x.fit.pairDCA < y.fit.pairDCA; });
    std::vector<bool> used(helices.size(), false);
    for (const auto& cand : candidates) {
      if (!m_keepOverlapping && (used[cand.a] || used[cand.b]))
        continue;
      used[cand.a] = used[cand.b] = true;
      const auto& f = cand.fit;
      auto vertex = vertices.create();
      vertex.setPosition({static_cast<float>(f.position.X()), static_cast<float>(f.position.Y()),
                          static_cast<float>(f.position.Z())});
      vertex.setChi2(static_cast<float>(f.pairDCA));
      vertex.setNdf(0);
      vertex.setPrimary(false);
      vertex.setAlgorithmType(m_algorithmType);
      for (float value : {2.f, static_cast<float>(sourceIndex[cand.a]), static_cast<float>(sourceIndex[cand.b]),
                          static_cast<float>(f.pairDCA), static_cast<float>(f.t[0]), static_cast<float>(f.t[1]),
                          f.upstream ? 1.f : 0.f, cand.region ? 1.f : 0.f, static_cast<float>(cand.chargeProduct),
                          f.betweenIPAndHits ? 1.f : 0.f, static_cast<float>(f.mass), static_cast<float>(f.cosOpening),
                          static_cast<float>(f.p[0]), static_cast<float>(f.p[1])})
        vertex.addToParameters(value);
      if (!m_allowMultipleVertices)
        break;
    }
    debug() << "Reconstructed " << vertices.size() << " helix pair vertices" << endmsg;
    return vertices;
  }

private:
  Config m_cfg;
  k4vertexing::VertexRegion m_regionCut;

  // --- field map ---
  Gaudi::Property<double> m_innerField{this, "InnerField", 0.0,
                                       "Bz inside the cylinder (FieldR, FieldHalfZ) [T] (must be set)"};
  Gaudi::Property<double> m_outerField{this, "OuterField", 0.0, "Bz outside the cylinder [T]"};
  Gaudi::Property<double> m_fieldR{this, "FieldR", 1.0e12,
                                   "Radius of the inner-field cylinder [mm] (default: InnerField everywhere)"};
  Gaudi::Property<double> m_fieldHalfZ{this, "FieldHalfZ", 1.0e12, "Half length of the inner-field cylinder [mm]"};
  // --- material map (energy loss) ---
  Gaudi::Property<double> m_shellInner{this, "ShellInner", 0.0, "Inner box (r and |z|) of the material shell [mm]"};
  Gaudi::Property<double> m_shellOuter{this, "ShellOuter", 0.0, "Outer box (r and |z|) of the material shell [mm]"};
  Gaudi::Property<double> m_shellHole{this, "ShellBeamHole", 0.0, "Beam-hole radius of the shell [mm]"};
  Gaudi::Property<double> m_dEdxShell{this, "dEdxShell", 0.0, "Energy loss in the shell [GeV/mm]; 0 = none"};
  Gaudi::Property<double> m_dEdxOutside{this, "dEdxOutside", 0.0,
                                        "Energy loss beyond the shell [GeV/mm]; 0 = none"};
  // --- track description ---
  Gaudi::Property<std::string> m_trackInput{this, "TrackInput", "hits",
                                            "hits: helix fitted to the hits; state: helix from a track state"};
  Gaudi::Property<int> m_stateLocation{this, "StateLocation", edm4hep::TrackState::AtFirstHit,
                                       "Track state used when TrackInput = state (1 AtIP, 2 AtFirstHit, 3 AtLastHit)"};
  Gaudi::Property<std::string> m_stateReference{
      this, "StateReference", "referencePoint",
      "Point the state parameters refer to: referencePoint (EDM4hep) or origin"};
  Gaudi::Property<bool> m_stateReversed{
      this, "StateReversed", false,
      "The state describes the track moving inward: phi, omega and tanLambda are turned to the outward "
      "direction"};
  // --- decay products ---
  Gaudi::Property<double> m_particleMass{this, "ParticleMass", 0.1056584,
                                         "Mass hypothesis of both tracks for the pair mass [GeV]"};
  Gaudi::Property<bool> m_planarHitsOnly{this, "PlanarHitsOnly", true,
                                         "Use only space-point hits: skip drift-chamber wire hits (SenseWireHit)"};
  // --- track model and propagation ---
  Gaudi::Property<unsigned> m_minHits{this, "MinHits", 3, "Minimum hits per track"};
  Gaudi::Property<unsigned> m_maxHitsInFit{this, "MaxHitsInFit", 0,
                                           "Fit only the N innermost hits (0 = all)"};
  Gaudi::Property<double> m_minXYSpan{this, "MinXYSpan", 30.0,
                                      "Minimum xy extent of the hits for a circle fit, else a 3D line [mm]"};
  Gaudi::Property<double> m_maxRadius{this, "MaxRadius", 2.0e6, "Larger fitted radius = straight line [mm]"};
  Gaudi::Property<double> m_step{this, "Step", 20.0, "Propagation step [mm]"};
  Gaudi::Property<double> m_backFactor{this, "BackFactor", 1.3,
                                       "Propagate inward up to this factor x distance of the first hit"};
  Gaudi::Property<double> m_upstreamTolerance{this, "UpstreamTolerance", 100.0,
                                              "Vertex may lie this far past a first hit [mm]"};
  // --- selection ---
  Gaudi::Property<std::string> m_region{
      this, "Region", "any",
      "Vertex region w.r.t. the cylinder (RegionR, RegionHalfZ): outside, barrel, endcap, inside, any"};
  Gaudi::Property<double> m_regionR{this, "RegionR", 0.0, "Radius of the region cylinder [mm]"};
  Gaudi::Property<double> m_regionHalfZ{this, "RegionHalfZ", 0.0, "Half length of the region cylinder [mm]"};
  Gaudi::Property<bool> m_requireUpstream{this, "RequireUpstream", true,
                                          "Vertex must lie before every track's first hit"};
  Gaudi::Property<bool> m_requireBetween{this, "RequireBetweenIPAndHits", true,
                                         "Vertex closer to the origin than every first hit, on the same side"};
  Gaudi::Property<double> m_maxPairDCA{this, "MaxPairDCA", 300.0,
                                       "Max distance between the two trajectories at the vertex [mm]"};
  Gaudi::Property<bool> m_requireOppositeCharge{this, "RequireOppositeCharge", true,
                                                "Keep only pairs with opposite rotation sense"};
  Gaudi::Property<double> m_minMass{this, "MinMass", 0.0, "Min pair mass at the vertex [GeV]"};
  Gaudi::Property<double> m_maxMass{this, "MaxMass", 1.0e9, "Max pair mass at the vertex [GeV]"};
  Gaudi::Property<double> m_minCosOpening{this, "MinCosOpening", -1.0, "Min cosine of the opening angle"};
  Gaudi::Property<bool> m_keepOverlapping{this, "KeepOverlappingPairs", false,
                                          "Write pairs sharing a track too (for cut studies)"};
  Gaudi::Property<bool> m_allowMultipleVertices{this, "AllowMultipleVertices", true,
                                                "Write every non-overlapping pair, not only the best"};
  Gaudi::Property<int> m_algorithmType{this, "AlgorithmType", 4, "edm4hep::Vertex algorithmType tag"};
};

DECLARE_COMPONENT(HelixPairVertexFinder)
