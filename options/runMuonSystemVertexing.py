"""
Run the two pair vertexers on IDEA muon-system tracks (decays outside the solenoid):
  HelixPairVertexFinder : 3D stepped helix, barrel and endcap
  ArcPairVertexFinder   : local arc fit with a chi2 vertex fit

All the IDEA numbers (field map, material, region, resolutions) are set here;
the C++ defaults describe no detector.

Author: Mahmoud Althakeel

Example
-------
    k4run runMuonSystemVertexing.py \
        --inputFile  hnl_tracks.root \
        --outputFile muon_vertices.root

Output collections:
  MuonSystemHelixVertices  (helix, algorithmType 4, unless --noHelix)
  MuonSystemVertices       (arc,   algorithmType 3, unless --noArc)

    --region  outside (default) | barrel | endcap | inside | any   (vertex region, both finders)
    --loose   no DCA / charge / upstream / between / region cut, all pairs kept
              (validation and cut scans; the quantities are in the vertex parameters)
"""

from Gaudi.Configuration import INFO, DEBUG
from Configurables import EventDataSvc, UniqueIDGenSvc
from k4FWCore import IOSvc, ApplicationMgr
from k4FWCore.parseArgs import parser

parser.add_argument("--inputFile", required=True)
parser.add_argument("--outputFile", required=True)
parser.add_argument("--inputTracks", default="StandaloneMuonTracks")
parser.add_argument("--region", default="outside", choices=["outside", "barrel", "endcap", "inside", "any"])
parser.add_argument("--loose", action="store_true", help="no selection, keep all pairs (validation, cut scans)")
parser.add_argument("--noArc", action="store_true", help="run only the helix vertexer")
parser.add_argument("--noHelix", action="store_true", help="run only the arc vertexer")
parser.add_argument("--noMaterial", action="store_true",
                    help="arc: no multiple scattering in the track errors (hit resolution + field term only)")
parser.add_argument("--fieldUncertainty", type=float, default=0.5,
                    help="arc: field mismatch [T] in the bending-plane error; 0 = off")
parser.add_argument("--noCurvature", action="store_true", help="arc: straight line in xy instead of a circle")
parser.add_argument("--debug", action="store_true", help="print the per-event counts of both finders")
args = parser.parse_args()

svc = IOSvc("IOSvc")
svc.Input = args.inputFile
svc.Output = args.outputFile

try:
    from Configurables import HelixPairVertexFinder, ArcPairVertexFinder
except ImportError:
    from DisplacedVertexing.DisplacedVertexingConf import HelixPairVertexFinder, ArcPairVertexFinder

# ---------------------------------------------------------------------------
# IDEA geometry used by both finders
# ---------------------------------------------------------------------------
SOLENOID_R = 2300.0       # mm, vertex region "outside" = beyond this cylinder
SOLENOID_HALF_Z = 2180.0  # mm

# ---------------------------------------------------------------------------
# 3D helix vertexer
# ---------------------------------------------------------------------------
helix = HelixPairVertexFinder(
    "HelixPairVertexFinder",
    InputTracks=[args.inputTracks],
    OutputVertices=["MuonSystemHelixVertices"],
    # field map of the IDEA_o1_v05 simulation: +2 T inside the coil volume,
    # -1.7 T return field everywhere else (where the muon chambers are)
    InnerField=2.0,
    OuterField=-1.7,
    FieldR=2300.0,
    FieldHalfZ=2616.0,        # 1.2 x coil half length
    # energy loss: dual-readout calorimeter shell + iron yoke beyond it
    # (dE/dx calibrated on truth, 1.36 x nominal in the calorimeter)
    ShellInner=2475.0,
    ShellOuter=4525.0,
    ShellBeamHole=252.0,
    dEdxShell=1.36e-3,        # GeV/mm
    dEdxOutside=1.97e-3,      # GeV/mm
    ParticleMass=0.1056584,   # muons
    # muon-system scale of the track model and propagation
    MinHits=3,
    MinXYSpan=30.0,
    Step=20.0,
    BackFactor=1.3,
    UpstreamTolerance=100.0,
    # selection
    Region=args.region,
    RegionR=SOLENOID_R,
    RegionHalfZ=SOLENOID_HALF_Z,
    MaxPairDCA=300.0,
    OutputLevel=INFO,
)

# ---------------------------------------------------------------------------
# arc vertexer
# ---------------------------------------------------------------------------
arc = ArcPairVertexFinder(
    "ArcPairVertexFinder",
    InputTracks=[args.inputTracks],
    OutputVertices=["MuonSystemVertices"],
    # StandaloneMuonTracks compute omega with 2 T (although the chambers sit in -1.7 T)
    BFieldForMomentum=2.0,
    MomentumResolution=0.20,  # sigma(p)/p of the standalone muon system
    HitResolution=1.0,        # mm
    FieldUncertainty=args.fieldUncertainty,
    # multiple scattering: dual-readout calorimeter shell + muon-system yoke
    UseMaterialEffects=not args.noMaterial,
    ShellX0=90.0,
    ShellInnerR=2400.0,
    ShellOuterR=4500.0,
    ShellInnerZ=2400.0,
    ShellOuterZ=4300.0,
    OutsideX0=40.0,
    OutsideThickness=2000.0,
    FitCurvature=not args.noCurvature,
    MinHits=3,
    # muon-system scale of the numerics
    JacobianStep=1.0,
    ConvergenceDistance=0.1,
    UpstreamTolerance=100.0,
    # selection (HNL mN = 20 GeV working point)
    Region=args.region,
    RegionR=SOLENOID_R,
    RegionHalfZ=SOLENOID_HALF_Z,
    MaxPairDCA=100.0,
    OutputLevel=INFO,
)

if args.loose:
    for f in (helix, arc):
        f.MaxPairDCA = 1.0e9
        f.RequireOppositeCharge = False
        f.AllowMultipleVertices = True
        f.KeepOverlappingPairs = True
        f.RequireBetweenIPAndHits = False
        f.RequireUpstream = False
        f.Region = "any"
    arc.MaxChi2 = 1.0e12
    arc.MinUpstreamDistance = -1.0e9

if args.debug:
    helix.OutputLevel = DEBUG
    arc.OutputLevel = DEBUG

algorithms = []
if not args.noHelix:
    algorithms.append(helix)
if not args.noArc:
    algorithms.append(arc)

ApplicationMgr(
    TopAlg=algorithms,
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc"), UniqueIDGenSvc("UniqueIDGenSvc")],
    OutputLevel=INFO,
)
