"""
Run ONLY the graph-based displaced vertex finder on an existing collection of
fitted tracks.

Use this when tracking has already been run: it reads `FittedTracks` from the
input file and writes the reconstructed displaced vertices, which makes it fast
to iterate on the finder's configuration without re-running digitisation,
track finding and track fitting each time.

Algorithm: J. Kriewald, arXiv:2510.00856 (github.com/jkriewald/delphes-LLP)
Port author: Mahmoud Althakeel

Example
-------
    k4run runGraphVertexing.py \
        --inputFile  reco_with_fitted_tracks.root \
        --outputFile graph_vertices.root

Useful switches
---------------
    --openChi2       disable every chi2 gate (diagnostic: see the full chi2
                     distribution instead of only what survives the cuts)
    --verboseKernel  print the finder's internal step-by-step decisions:
                     prompt/displaced split, and for each track pair whether it
                     was accepted or which guard rejected it
"""
from Gaudi.Configuration import INFO
from Configurables import EventDataSvc, UniqueIDGenSvc
from k4FWCore import IOSvc, ApplicationMgr
from k4FWCore.parseArgs import parser

parser.add_argument("--inputFile", required=True, help="file containing FittedTracks")
parser.add_argument("--outputFile", required=True)
parser.add_argument("--inputTracks", default="FittedTracks", help="input track collection")
parser.add_argument("--openChi2", action="store_true", help="accept all chi2 (diagnostic)")
parser.add_argument("--verboseKernel", action="store_true", help="print the finder's internal steps")
parser.add_argument("--trackState", type=int, default=1,
                    help="TrackState location: 1 = AtIP, 2 = AtFirstHit, 3 = AtLastHit (states at a hit are moved to the origin)")
parser.add_argument("--minTracks", type=int, default=2, help="minimum tracks per vertex")
parser.add_argument("--muonConvention", action="store_true",
                    help="input from an older StandaloneMuonTracking version: fix d0/omega sign and move z0 to the perigee")
parser.add_argument("--noPromptSplit", action="store_true",
                    help="send every track to the graph (no IP-significance / PV-association veto). "
                         "Needed for muon-system tracks: their cm-scale d0/z0 errors make boosted "
                         "LLP daughters look PV-compatible")
args = parser.parse_args()

svc = IOSvc("IOSvc")
svc.Input = args.inputFile
svc.Output = args.outputFile

id_service = UniqueIDGenSvc("UniqueIDGenSvc")

try:
    from Configurables import GraphDisplacedVertexFinder
except ImportError:
    # used when the Configurables database is not available
    from DisplacedVertexing.DisplacedVertexingConf import GraphDisplacedVertexFinder

# With --openChi2 every chi2 threshold is pushed far beyond any realistic value,
# so nothing is rejected on chi2 and the whole distribution can be inspected.
CHI2 = 1.0e12 if args.openChi2 else 9.21

finder = GraphDisplacedVertexFinder(
    "GraphDisplacedVertexFinder",
    InputTracks=[args.inputTracks],
    OutputDisplacedVertices=["DisplacedVertices"],
    OutputPrimaryVertex=["GraphPrimaryVertex"],

    # ---- track input -----------------------------------------------------------
    TrackStateLocation=args.trackState,
    MuonSystemTrackConvention=args.muonConvention,
    SkipFailedFitTracks=True,

    # ---- detector --------------------------------------------------------------------
    MagneticFieldZ=2.0,                # T
    BeamSpotSigmaX=0.1,                # mm
    BeamSpotSigmaY=0.1,                # mm
    BeamSpotSigmaZ=1.0,                # mm
    TimingR=2.25,                      # m
    TimingZ=2.5,                       # m

    # ---- physics ---------------------------------------------------------------------------
    ParticleMass=0.1056583755,         # GeV

    # ---- track preselection -----------------------------------------------------------------
    MinTrackIPSig=3.0,
    MinD0=0.0,                         # mm
    MinZ0=0.0,                         # mm
    MinTrackPT=0.1,                    # GeV

    # ---- algorithm: primary vertex -------------------------------------------------------------
    PVChi2NDFMax=CHI2,
    UsePVCut=True,
    PVCutChi2=CHI2,

    # ---- algorithm: clustering (settings for two-track decays) -----------------------------------
    SeedSelector=0,
    MinSeedPT=1.0,                     # GeV
    MinSeedIPSig=5.0,
    MinTracks=args.minTracks,
    kNN=-1,                            # no kNN pruning
    MinSupport=0,                      # a pair has no triplet support; must stay 0 for MinTracks=2
    UsePairGuard=True,
    PruneBridges=True,

    # ---- algorithm: vertex fit --------------------------------------------------------------------
    Chi2PairCut=CHI2,
    Chi2TripletCut=CHI2,
    WeightCut=0.1,
    Chi2NDFMax=CHI2,
    BridgeCut=(1.0e12 if args.openChi2 else 3.0),
    Chi2Sigmoid=9.21,
    Beta=2.0,
    UseSelfSeeding=True,               # the seeding path that works for two tracks
    Chi2AssociationMax=(1.0e12 if args.openChi2 else 4.0),
    CorrectMomenta=True,

    # ---- algorithm: timing (off: fitted tracks carry no per-hit time) -------------------------------
    UseTiming=False,
    TimeGate=9.21,
    TimeVarianceCap=4.0,

    # ---- output ----------------------------------------------------------------------------
    AlgorithmType=2,
    VerboseKernel=args.verboseKernel,
    OutputLevel=INFO,
)
if args.noPromptSplit:
    finder.MinTrackIPSig = -1.0
    finder.MinD0 = -1.0
    finder.MinZ0 = -1.0
    finder.UsePVCut = False

ApplicationMgr(
    TopAlg=[finder],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc"), id_service],
    OutputLevel=INFO,
)
