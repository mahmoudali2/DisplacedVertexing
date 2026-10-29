"""
Run ArcPairVertexFinder: two-track vertices from local arc fits to the tracks' hits.

Every property is listed below with its default value, grouped by use, so the
steering shows the full configuration. The defaults suit a silicon tracker in a
uniform solenoid field; change the detector block for another detector.

Example
-------
    k4run runArcPairVertexing.py \
        --inputFile  tracks.root \
        --outputFile vertices.root \
        --inputTracks SiTracks_Refitted --field 2.0
"""

from Gaudi.Configuration import INFO
from Configurables import EventDataSvc, UniqueIDGenSvc
from k4FWCore import IOSvc, ApplicationMgr
from k4FWCore.parseArgs import parser

parser.add_argument("--inputFile", required=True)
parser.add_argument("--outputFile", required=True)
parser.add_argument("--inputTracks", default="Tracks", help="input track collection")
parser.add_argument("--outputVertices", default="ArcPairVertices", help="output vertex collection")
parser.add_argument("--field", type=float, default=2.0,
                    help="field [T] the tracking used to compute the tracks' omega")
parser.add_argument("--stateLocation", type=int, default=2,
                    help="track state whose omega gives the momentum (1 AtIP, 2 AtFirstHit, 3 AtLastHit)")
args = parser.parse_args()

inputOutputService = IOSvc("IOSvc")
inputOutputService.Input = args.inputFile
inputOutputService.Output = args.outputFile

try:
    from Configurables import ArcPairVertexFinder
except ImportError:
    from DisplacedVertexing.DisplacedVertexingConf import ArcPairVertexFinder

arcPairVertexFinder = ArcPairVertexFinder(
    "ArcPairVertexFinder",
    InputTracks=[args.inputTracks],
    OutputVertices=[args.outputVertices],

    # ---- track input -----------------------------------------------------------
    MomentumStateLocation=args.stateLocation,
    MinHits=3,
    PlanarHitsOnly=True,               # skip drift-chamber wire hits

    # ---- detector: magnetic field ----------------------------------------------------
    BFieldForMomentum=args.field,      # T
    FieldUncertainty=0.0,              # T, 0 = uniform field

    # ---- detector: resolutions ----------------------------------------------------------
    HitResolution=0.01,                # mm per coordinate
    MomentumResolution=0.01,           # sigma(p)/p of the input tracks

    # ---- detector: material (multiple scattering in the track errors) ------------------
    # off: inside a tracker there is little material between the decay and the hits
    UseMaterialEffects=False,
    ShellX0=0.0,
    ShellInnerR=0.0,                   # mm
    ShellOuterR=0.0,                   # mm
    ShellInnerZ=0.0,                   # mm
    ShellOuterZ=0.0,                   # mm
    OutsideX0=0.0,
    OutsideThickness=1.0,              # mm

    # ---- physics ---------------------------------------------------------------------------
    DefaultPt=10.0,                    # GeV, used when a state has omega = 0
    MinMomentum=2.0,                   # GeV, floor for the scattering estimate

    # ---- algorithm: track model and vertex fit -----------------------------------------------
    FitCurvature=True,
    UpstreamTolerance=300.0,           # mm, vertex may lie this far past a first hit
    MaterialSteps=40,
    MaxIterations=12,
    JacobianStep=0.01,                 # mm, scale to the vertex resolution
    ConvergenceDistance=0.001,         # mm

    # ---- vertex selection ------------------------------------------------------------------
    MaxPairDCA=5.0,                    # mm
    MaxChi2=1.0e9,
    RequireOppositeCharge=True,
    RequireUpstream=False,
    MinUpstreamDistance=0.0,           # mm
    RequireBetweenIPAndHits=False,
    Region="any",                      # outside, barrel, endcap, inside, any
    RegionR=0.0,                       # mm
    RegionHalfZ=0.0,                   # mm

    # ---- output ----------------------------------------------------------------------------
    AllowMultipleVertices=True,
    KeepOverlappingPairs=False,
    AlgorithmType=3,
    OutputLevel=INFO,
)

ApplicationMgr(
    TopAlg=[arcPairVertexFinder],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc"), UniqueIDGenSvc("UniqueIDGenSvc")],
    OutputLevel=INFO,
)
