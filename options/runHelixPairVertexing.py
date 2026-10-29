"""
Run HelixPairVertexFinder: two-track vertices from the tracks' fitted track states.

Every property is listed below with its default value, grouped by use, so the
steering shows the full configuration. The defaults suit a silicon tracker in a
uniform solenoid field; change the detector block for another detector.

Example
-------
    k4run runHelixPairVertexing.py \
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
parser.add_argument("--outputVertices", default="HelixPairVertices", help="output vertex collection")
parser.add_argument("--field", type=float, default=2.0, help="solenoid field [T]")
parser.add_argument("--stateLocation", type=int, default=2,
                    help="track state used (1 AtIP, 2 AtFirstHit, 3 AtLastHit)")
args = parser.parse_args()

inputOutputService = IOSvc("IOSvc")
inputOutputService.Input = args.inputFile
inputOutputService.Output = args.outputFile

try:
    from Configurables import HelixPairVertexFinder
except ImportError:
    from DisplacedVertexing.DisplacedVertexingConf import HelixPairVertexFinder

helixPairVertexFinder = HelixPairVertexFinder(
    "HelixPairVertexFinder",
    InputTracks=[args.inputTracks],
    OutputVertices=[args.outputVertices],

    # ---- track input -----------------------------------------------------------
    StateLocation=args.stateLocation,
    StateReference="referencePoint",   # EDM4hep convention
    StateReversed=False,
    UseCovariance=True,                # vertex fit with the state covariances (chi2, covariance)

    # ---- detector: magnetic field ----------------------------------------------------
    # one field everywhere: leave FieldR and FieldHalfZ at their (very large) defaults
    InnerField=args.field,             # T
    OuterField=0.0,                    # T, used outside (FieldR, FieldHalfZ)
    FieldR=1.0e12,                     # mm
    FieldHalfZ=1.0e12,                 # mm

    # ---- detector: material (energy loss added back while propagating) --------------
    # off: inside a tracker the fitted state already includes the material
    ShellInner=0.0,                    # mm
    ShellOuter=0.0,                    # mm
    ShellBeamHole=0.0,                 # mm
    dEdxShell=0.0,                     # GeV/mm
    dEdxOutside=0.0,                   # GeV/mm

    # ---- physics ---------------------------------------------------------------------------
    ParticleMass=0.1056584,            # GeV, mass hypothesis of both tracks (pair mass only)

    # ---- algorithm: propagation and vertex search ------------------------------------------
    Step=5.0,                          # mm
    BackFactor=1.3,                    # propagate inward up to 1.3 x distance of the first hit
    UpstreamTolerance=300.0,           # mm, vertex may lie this far past a first hit

    # ---- vertex selection ------------------------------------------------------------------
    MaxPairDCA=5.0,                    # mm
    RequireOppositeCharge=True,
    RequireUpstream=False,
    RequireBetweenIPAndHits=False,
    Region="any",                      # outside, barrel, endcap, inside, any
    RegionR=0.0,                       # mm
    RegionHalfZ=0.0,                   # mm
    MinMass=0.0,                       # GeV
    MaxMass=1.0e9,                     # GeV
    MinCosOpening=-1.0,
    MaxChi2=1.0e9,                     # fitted vertices only (ndf = 1)

    # ---- output ----------------------------------------------------------------------------
    AllowMultipleVertices=True,
    KeepOverlappingPairs=False,
    AlgorithmType=4,
    OutputLevel=INFO,
)

ApplicationMgr(
    TopAlg=[helixPairVertexFinder],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc"), UniqueIDGenSvc("UniqueIDGenSvc")],
    OutputLevel=INFO,
)
