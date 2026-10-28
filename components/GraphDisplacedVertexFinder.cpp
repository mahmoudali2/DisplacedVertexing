/*
 * Gaudi component: graph-based displaced-vertex finding on EDM4hep tracks.
 *
 * Key4hep port of the Delphes module by Jonathan Kriewald:
 *   J. Kriewald, "No Track left behind: Graph-based Vertexing for long-lived
 *   Particle Reconstruction", arXiv:2510.00856
 *   https://github.com/jkriewald/delphes-LLP
 * The algorithm is his; this file is only the EDM4hep/Gaudi interface.
 *
 * Port author: Mahmoud Althakeel
 *
 * Input : edm4hep::TrackCollection  (fitted tracks, e.g. GenfitTrackFitter output)
 * Output: edm4hep::VertexCollection (displaced vertices) and the primary vertex
 */
#include "Gaudi/Property.h"
#include "k4FWCore/Transformer.h"

#include <edm4hep/TrackCollection.h>
#include <edm4hep/VertexCollection.h>

#include "EDM4hepHelixConverter.h"
#include "GraphTrack.h"
#include "GraphVertexKernel.h"

#include <cmath>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

using namespace k4vertexing;

struct GraphDisplacedVertexFinder final
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::VertexCollection, edm4hep::VertexCollection>(
          const edm4hep::TrackCollection&)> {

  GraphDisplacedVertexFinder(const std::string& name, ISvcLocator* serviceLocator)
      : MultiTransformer(name, serviceLocator, {KeyValues("InputTracks", {"FittedTracks"})},
                         {KeyValues("OutputDisplacedVertices", {"DisplacedVertices"}),
                          KeyValues("OutputPrimaryVertex", {"GraphPrimaryVertex"})}) {}

  StatusCode initialize() override {
    const auto status = MultiTransformer::initialize();
    if (!status.isSuccess())
      return status;
    if (m_minTracks < 2) {
      error() << "MinTracks must be at least 2" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_magneticFieldZ == 0.) {
      error() << "MagneticFieldZ must be non-zero" << endmsg;
      return StatusCode::FAILURE;
    }
    info() << "Graph displaced-vertex finder configured with MinTracks = " << m_minTracks.value()
           << " (2 = two-track vertices such as LLP -> mu mu are accepted)" << endmsg;
    return StatusCode::SUCCESS;
  }

  /**
   * Bring a StandaloneMuonTracking state into the EDM4hep perigee convention.
   *
   * The outer-segment state built by StandaloneMuonTracking has:
   *  - d0 and phi at the perigee with respect to the origin (correct point), but
   *    d0 and omega with the opposite sign to EDM4hep: its circle only passes
   *    through the first muon hit once both signs are flipped (measured: median
   *    0.8 mm with the flip, ~1 m without);
   *  - z0 and tanLambda from a straight-line fit of z versus arc length with
   *    s = 0 at the FIRST muon hit, i.e. z0 is the z at that hit (which is also
   *    the state's referencePoint), not the z at the perigee.
   * The vertexer needs all five parameters at the same point, so:
   *    d0 -> -d0,  omega -> -omega,  z0 -> z0 - tanLambda * s*
   * where s* is the transverse arc length from the perigee to the reference
   * point along the (sign-corrected) circle.
   *
   * Covariance: the sign flips are a diagonal Jacobian diag(-1,1,-1,1,1); the z0
   * shift adds dz0'/dtanLambda = -s* (dependence of s* on the circle parameters
   * is neglected).
   */
  static void correctMuonSystemState(edm4hep::TrackState& state) {
    const double d0 = -state.D0;
    const double omega = -state.omega;
    const double phi = state.phi;
    const double tanLambda = state.tanLambda;

    // circle centre of the EDM4hep helix about the origin
    const double radius = 1.0 / omega; // signed
    const double centreX = (radius - d0) * std::sin(phi);
    const double centreY = -(radius - d0) * std::cos(phi);
    // position angle about the centre at the perigee and at the reference point;
    // along the helix that angle advances as (-omega * s)
    const double perigeeX = -d0 * std::sin(phi), perigeeY = d0 * std::cos(phi);
    const double anglePerigee = std::atan2(perigeeY - centreY, perigeeX - centreX);
    const double angleRef = std::atan2(state.referencePoint.y - centreY, state.referencePoint.x - centreX);
    double dAngle = anglePerigee - angleRef;
    while (dAngle > M_PI)
      dAngle -= 2.0 * M_PI;
    while (dAngle < -M_PI)
      dAngle += 2.0 * M_PI;
    const double sStar = dAngle / omega;

    // covariance: unpack, apply J, repack
    TMatrixDSym cov = unpackCovariance(state);
    TMatrixD jacobian(kNPerigee, kNPerigee);
    jacobian.UnitMatrix();
    jacobian(kD0, kD0) = -1.0;
    jacobian(kOmega, kOmega) = -1.0;
    jacobian(kZ0, kTanLambda) = -sStar;
    cov.Similarity(jacobian); // J * C * J^T

    state.D0 = static_cast<float>(d0);
    state.omega = static_cast<float>(omega);
    state.Z0 = static_cast<float>(state.Z0 - tanLambda * sStar);
    int packed = 0;
    for (int row = 0; row < kNPerigee; ++row)
      for (int column = 0; column <= row; ++column, ++packed)
        state.covMatrix[packed] = static_cast<float>(cov(row, column));
  }

  /// Build the kernel's track view from one EDM4hep track.
  bool makeTrack(const edm4hep::Track& track, int index, GraphTrack& out) const {
    if (m_skipFailedFits && (track.getNdf() <= 0 || track.getChi2() < 0.))
      return false;

    const edm4hep::TrackState* state = nullptr;
    for (const auto& candidate : track.getTrackStates())
      if (candidate.location == m_trackStateLocation) {
        state = &candidate;
        break;
      }
    if (!state || state->omega == 0.f)
      return false;

    // Muon-system tracks need two corrections before they describe one
    // consistent helix in the EDM4hep convention (see correctMuonSystemState).
    edm4hep::TrackState corrected = *state;
    if (m_muonSystemConvention)
      correctMuonSystemState(corrected);
    else
      moveStateToOrigin(corrected); // the Delphes perigee is about the origin

    const TVectorD parameters = toDelphesPerigee(corrected);
    const TMatrixDSym covariance = toDelphesPerigeeCovariance(corrected);
    if (!hasUsableCovariance(covariance))
      return false;

    out.D0 = parameters(0);
    out.Phi = parameters(1);
    out.C = parameters(2);
    out.DZ = parameters(3);
    out.CtgTheta = parameters(4);
    out.TrackCovariance = covariance;

    const double pt = transverseMomentum(state->omega, m_magneticFieldZ);
    const double pz = pt * state->tanLambda;
    const double mass = m_particleMass;
    out.PT = pt;
    out.Mass = mass;
    out.Charge = (corrected.omega > 0.f) ? +1 : -1;
    out.Momentum.SetXYZM(pt * std::cos(corrected.phi), pt * std::sin(corrected.phi), pz, mass);
    // No per-hit timing in this input; ErrorT <= 0 tells the kernel to skip timing.
    out.Position.SetXYZT(state->referencePoint.x, state->referencePoint.y, state->referencePoint.z, 0.0);
    out.ErrorT = -1.0;
    out.sourceIndex = index;
    return true;
  }

  std::tuple<edm4hep::VertexCollection, edm4hep::VertexCollection>
  operator()(const edm4hep::TrackCollection& tracks) const override {
    edm4hep::VertexCollection displacedVertices;
    edm4hep::VertexCollection primaryVertices;

    // ---- build the kernel's track views -----------------------------------
    std::vector<std::unique_ptr<GraphTrack>> storage;
    std::vector<const Candidate*> view;
    storage.reserve(tracks.size());
    view.reserve(tracks.size());
    int index = 0;
    for (const auto& track : tracks) {
      auto entry = std::make_unique<GraphTrack>();
      if (makeTrack(track, index, *entry)) {
        view.push_back(entry.get());
        storage.push_back(std::move(entry));
      }
      ++index;
    }
    debug() << "Using " << view.size() << " of " << tracks.size() << " tracks" << endmsg;
    if (view.size() < static_cast<std::size_t>(m_minTracks))
      return std::make_tuple(std::move(displacedVertices), std::move(primaryVertices));

    // ---- configure and run -------------------------------------------------
    GraphVertexKernel kernel;
    configure(kernel);
    kernel.Init();
    kernel.SetTracks(view);
    kernel.Process();

    // ---- primary vertex ----------------------------------------------------
    if (kernel.HasPrimaryVertex()) {
      const TVector3& position = kernel.GetPrimaryVertex();
      const TMatrixDSym& covariance = kernel.GetPrimaryVertexCov();
      auto vertex = primaryVertices.create();
      vertex.setPosition({static_cast<float>(position.X()), static_cast<float>(position.Y()),
                          static_cast<float>(position.Z())});
      vertex.setCovMatrix({static_cast<float>(covariance(0, 0)), static_cast<float>(covariance(0, 1)),
                           static_cast<float>(covariance(1, 1)), static_cast<float>(covariance(0, 2)),
                           static_cast<float>(covariance(1, 2)), static_cast<float>(covariance(2, 2))});
      vertex.setPrimary(true);
      vertex.setAlgorithmType(m_algorithmType);
    }

    // ---- displaced vertices -------------------------------------------------
    for (const auto& cluster : kernel.GetClusters()) {
      if (cluster.ndf <= 0 || cluster.tracks.size() < static_cast<std::size_t>(m_minTracks))
        continue;
      auto vertex = displacedVertices.create();
      vertex.setPosition({static_cast<float>(cluster.fittedPos(0)), static_cast<float>(cluster.fittedPos(1)),
                          static_cast<float>(cluster.fittedPos(2))});
      const TMatrixDSym& covariance = cluster.fittedCov;
      vertex.setCovMatrix({static_cast<float>(covariance(0, 0)), static_cast<float>(covariance(0, 1)),
                           static_cast<float>(covariance(1, 1)), static_cast<float>(covariance(0, 2)),
                           static_cast<float>(covariance(1, 2)), static_cast<float>(covariance(2, 2))});
      vertex.setChi2(static_cast<float>(cluster.chi2));
      vertex.setNdf(cluster.ndf);
      vertex.setPrimary(false);
      vertex.setAlgorithmType(m_algorithmType);
      // edm4hep::Vertex relates to ReconstructedParticle, not to Track, so the
      // association back to the input tracks is recorded as the track indices
      // in the vertex parameters. parameters[0] is the number of tracks,
      // followed by that many indices into the input TrackCollection.
      vertex.addToParameters(static_cast<float>(cluster.tracks.size()));
      for (const auto* clusterTrack : cluster.tracks) {
        if (clusterTrack && clusterTrack->sourceIndex >= 0 &&
            clusterTrack->sourceIndex < static_cast<int>(tracks.size()))
          vertex.addToParameters(static_cast<float>(clusterTrack->sourceIndex));
      }
    }
    debug() << "Reconstructed " << displacedVertices.size() << " displaced vertices" << endmsg;
    return std::make_tuple(std::move(displacedVertices), std::move(primaryVertices));
  }

private:
  /// Copy the Gaudi properties onto the kernel's public configuration.
  void configure(GraphVertexKernel& kernel) const {
    kernel.fMinsDisp = m_minTrackIPSig;
    kernel.fMinD0 = m_minD0;
    kernel.fMinZ0 = m_minZ0;
    kernel.fMinTrackPT = m_minTrackPT;

    kernel.fBeamSpotSigmaX = m_beamSpotSigmaX;
    kernel.fBeamSpotSigmaY = m_beamSpotSigmaY;
    kernel.fBeamSpotSigmaZ = m_beamSpotSigmaZ;
    kernel.fPVChi2NDFMax = m_pvChi2NDFMax;
    kernel.fUsePVCut = m_usePVCut;
    kernel.fPVCutChi2 = m_pvCutChi2;

    kernel.fSeedSelector = m_seedSelector;
    kernel.fMinSeedPT = m_minSeedPT;
    kernel.fMinSeedIPSig = m_minSeedIPSig;
    kernel.fMinTracks = m_minTracks;

    kernel.fkNN = m_kNN;
    kernel.fMinSupport = m_minSupport;
    kernel.fUsePairGuard = m_usePairGuard;
    kernel.fPruneBridges = m_pruneBridges;

    kernel.fchi2PairCut = m_chi2PairCut;
    kernel.fchi2TripCut = m_chi2TripletCut;
    kernel.fweightCut = m_weightCut;
    kernel.fMTChi2NDFMax = m_chi2NDFMax;
    kernel.fBridgeCut = m_bridgeCut;

    kernel.fUseTiming = m_useTiming;
    kernel.fTimeGate = m_timeGate;
    kernel.fTimeVariance = m_timeVarianceCap;
    kernel.fc0 = m_chi2_0;
    kernel.fbeta = m_beta;
    kernel.fFitOpts.useSelfSeeding = m_useSelfSeeding;

    kernel.fChi2AssociationMax = m_chi2AssociationMax;
    kernel.fCorrectMomenta = m_correctMomenta;
    kernel.fVerbose = m_verbose;
    kernel.fTimingR = m_timingR;
    kernel.fTimingZ = m_timingZ;
  }

  // ===== track input =====================================================================
  Gaudi::Property<int> m_trackStateLocation{this, "TrackStateLocation", edm4hep::TrackState::AtIP,
                                            "TrackState to read (1 = AtIP, 3 = AtLastHit)"};
  Gaudi::Property<bool> m_muonSystemConvention{
      this, "MuonSystemTrackConvention", false,
      "Only for files from older StandaloneMuonTracking versions, whose AtLastHit state has inverted d0/omega "
      "signs and z0 at the first hit (see correctMuonSystemState). Off for CombinatorialTrackFinder output"};
  Gaudi::Property<bool> m_skipFailedFits{this, "SkipFailedFitTracks", true,
                                         "Ignore tracks whose upstream fit failed (chi2 < 0 or ndf <= 0)"};

  // ===== detector ========================================================================
  Gaudi::Property<double> m_magneticFieldZ{this, "MagneticFieldZ", 2.0, "Solenoid field [T]"};
  Gaudi::Property<double> m_beamSpotSigmaX{this, "BeamSpotSigmaX", 0.1, "Beam-spot sigma x [mm]"};
  Gaudi::Property<double> m_beamSpotSigmaY{this, "BeamSpotSigmaY", 0.1, "Beam-spot sigma y [mm]"};
  Gaudi::Property<double> m_beamSpotSigmaZ{this, "BeamSpotSigmaZ", 1.0, "Beam-spot sigma z [mm]"};
  Gaudi::Property<double> m_timingR{this, "TimingR", 2.25, "Timing layer radius [m]"};
  Gaudi::Property<double> m_timingZ{this, "TimingZ", 2.5, "Timing layer half-length [m]"};

  // ===== physics =========================================================================
  Gaudi::Property<double> m_particleMass{this, "ParticleMass", 0.1056583755, "Mass hypothesis [GeV]"};

  // ===== track preselection ==============================================================
  Gaudi::Property<double> m_minTrackIPSig{this, "MinTrackIPSig", 3.0, "Minimum impact-parameter significance"};
  Gaudi::Property<double> m_minD0{this, "MinD0", 0.0, "Minimum |d0| [mm]"};
  Gaudi::Property<double> m_minZ0{this, "MinZ0", 0.0, "Minimum |z0| [mm]"};
  Gaudi::Property<double> m_minTrackPT{this, "MinTrackPT", 0.1, "Minimum track pT [GeV]"};

  // ===== algorithm: primary vertex =======================================================
  Gaudi::Property<double> m_pvChi2NDFMax{this, "PVChi2NDFMax", 9.0, "Max chi2/ndf for the primary-vertex fit"};
  Gaudi::Property<bool> m_usePVCut{this, "UsePVCut", true, "Split prompt/displaced tracks using the PV"};
  Gaudi::Property<double> m_pvCutChi2{this, "PVCutChi2", 9.0, "chi2 for associating a track to the PV"};

  // ===== algorithm: clustering ===========================================================
  Gaudi::Property<int> m_seedSelector{this, "SeedSelector", 0, "0: min pT, 1: IP significance, 2: high-pT lepton"};
  Gaudi::Property<double> m_minSeedPT{this, "MinSeedPT", 1.0, "Minimum seed-track pT [GeV]"};
  Gaudi::Property<double> m_minSeedIPSig{this, "MinSeedIPSig", 5.0, "Minimum seed IP significance"};
  Gaudi::Property<int> m_minTracks{this, "MinTracks", 2,
                                   "Minimum tracks per vertex. 2 keeps two-track LLP decays such as "
                                   "N -> mu mu nu; the original paper's analysis used 3"};
  Gaudi::Property<int> m_kNN{this, "kNN", -1, "kNN graph pruning (-1 = off)"};
  Gaudi::Property<int> m_minSupport{this, "MinSupport", 0,
                                    "Minimum triplet support for pruning (0 = off; must stay 0 for two-track "
                                    "vertices, which have no triplet support)"};
  Gaudi::Property<bool> m_usePairGuard{this, "UsePairGuard", true, "First-hit guard on pairs"};
  Gaudi::Property<bool> m_pruneBridges{this, "PruneBridges", true, "Prune weakly supported graph bridges"};

  // ===== algorithm: vertex fit ===========================================================
  Gaudi::Property<double> m_chi2PairCut{this, "Chi2PairCut", 9.21, "chi2 cut on two-track fits"};
  Gaudi::Property<double> m_chi2TripletCut{this, "Chi2TripletCut", 9.21, "chi2 cut on three-track fits"};
  Gaudi::Property<double> m_weightCut{this, "WeightCut", 0.1, "Minimum robust track weight"};
  Gaudi::Property<double> m_chi2NDFMax{this, "Chi2NDFMax", 9.21, "Max chi2/ndf of an accepted vertex"};
  Gaudi::Property<double> m_bridgeCut{this, "BridgeCut", 3.0, "chi2 cut for graph bridges"};
  Gaudi::Property<double> m_chi2_0{this, "Chi2Sigmoid", 9.21, "Sigmoid midpoint of the robust weight"};
  Gaudi::Property<double> m_beta{this, "Beta", 2.0, "Sigmoid width of the robust weight"};
  Gaudi::Property<bool> m_useSelfSeeding{this, "UseSelfSeeding", true, "Euclidean self-seeding (needed for pairs)"};
  Gaudi::Property<double> m_chi2AssociationMax{this, "Chi2AssociationMax", 4.0, "Max chi2 to attach leftover tracks"};
  Gaudi::Property<bool> m_correctMomenta{this, "CorrectMomenta", true, "Refit track momenta at the vertex"};

  // ===== algorithm: timing (off by default: fitted tracks carry no per-hit time) =========
  Gaudi::Property<bool> m_useTiming{this, "UseTiming", false, "Use track timing in the fit"};
  Gaudi::Property<double> m_timeGate{this, "TimeGate", 9.21, "chi2 gate on the time residual"};
  Gaudi::Property<double> m_timeVarianceCap{this, "TimeVarianceCap", 4.0, "Cap on the time variance"};

  // ===== output ==========================================================================
  Gaudi::Property<int> m_algorithmType{this, "AlgorithmType", 2, "Value written to Vertex::algorithmType"};
  Gaudi::Property<bool> m_verbose{this, "VerboseKernel", false, "Verbose printout from the ported kernel"
};
};

DECLARE_COMPONENT(GraphDisplacedVertexFinder)
