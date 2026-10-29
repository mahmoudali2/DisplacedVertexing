"""
Run the two pair vertexers on IDEA muon-system tracks (decays outside the solenoid):
  HelixPairVertexFinder : 3D stepped helix from the track state, barrel and endcap
  ArcPairVertexFinder   : local arc fit to the hits with a chi2 vertex fit

The muon-system tracks must come from the current CombinatorialTrackFinder, whose
stored track state follows the EDM4hep convention.

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
parser.add_argument("--useCovariance", action="store_true",
                    help="helix: fit the vertex with the track-state covariances (tracks fitted with GenFit)")
parser.add_argument("--stateLocation", type=int, default=2,
                    help="helix: track state used (1 AtIP, 2 AtFirstHit, 3 AtLastHit)")
args = parser.parse_args()

inputOutputService = IOSvc("IOSvc")
inputOutputService.Input = args.inputFile
inputOutputService.Output = args.outputFile

try:
    from Configurables import HelixPairVertexFinder, ArcPairVertexFinder
except ImportError:
    from DisplacedVertexing.DisplacedVertexingConf import HelixPairVertexFinder, ArcPairVertexFinder

# ---------------------------------------------------------------------------
# IDEA geometry used by both finders
# ---------------------------------------------------------------------------
SOLENOID_RADIUS = 2300.0       # mm, vertex region "outside" = beyond this cylinder
SOLENOID_HALF_Z = 2180.0  # mm

# ---------------------------------------------------------------------------
# 3D helix vertexer
# ---------------------------------------------------------------------------
helixPairVertexFinder = HelixPairVertexFinder(
    "HelixPairVertexFinder",
    InputTracks=[args.inputTracks],
    OutputVertices=["MuonSystemHelixVertices"],

    # ---- track input -----------------------------------------------------------
    StateLocation=args.stateLocation,
    # the CTF circle-fit state carries fixed covariance floors, not the covariance of the
    # fit: weighting the tracks with it makes the vertex worse, so the closest-approach
    # midpoint is used (--useCovariance for tracks fitted with GenFit)
    UseCovariance=args.useCovariance,

    # ---- detector: magnetic field ----------------------------------------------------
    # field map of the IDEA_o1_v05 simulation: +2 T inside the coil volume,
    # -1.7 T return field everywhere else (where the muon chambers are)
    InnerField=2.0,                    # T
    OuterField=-1.7,                   # T
    FieldR=2300.0,                     # mm
    FieldHalfZ=2616.0,                 # mm, 1.2 x coil half length

    # ---- detector: material (energy loss added back while propagating) --------------
    # dual-readout calorimeter shell + iron yoke beyond it
    # (dE/dx calibrated on truth, 1.36 x nominal in the calorimeter)
    ShellInner=2475.0,                 # mm
    ShellOuter=4525.0,                 # mm
    ShellBeamHole=252.0,               # mm
    dEdxShell=1.36e-3,                 # GeV/mm
    dEdxOutside=1.97e-3,               # GeV/mm

    # ---- physics ---------------------------------------------------------------------------
    ParticleMass=0.1056584,            # GeV, muons

    # ---- algorithm: propagation and vertex search ------------------------------------------
    Step=20.0,                         # mm
    BackFactor=1.3,
    UpstreamTolerance=100.0,           # mm

    # ---- vertex selection ------------------------------------------------------------------
    MaxPairDCA=300.0,                  # mm
    Region=args.region,
    RegionR=SOLENOID_RADIUS,
    RegionHalfZ=SOLENOID_HALF_Z,
    OutputLevel=INFO,
)

# ---------------------------------------------------------------------------
# arc vertexer
# ---------------------------------------------------------------------------
arcPairVertexFinder = ArcPairVertexFinder(
    "ArcPairVertexFinder",
    InputTracks=[args.inputTracks],
    OutputVertices=["MuonSystemVertices"],

    # ---- track input -----------------------------------------------------------
    MinHits=3,

    # ---- detector: magnetic field ----------------------------------------------------
    BFieldForMomentum=2.0,             # T, field used by the tracking to compute omega
    FieldUncertainty=args.fieldUncertainty,

    # ---- detector: resolutions ----------------------------------------------------------
    HitResolution=1.0,                 # mm
    MomentumResolution=0.20,           # sigma(p)/p of the standalone muon system

    # ---- detector: material (multiple scattering in the track errors) ------------------
    # dual-readout calorimeter shell + muon-system yoke
    UseMaterialEffects=not args.noMaterial,
    ShellX0=90.0,
    ShellInnerR=2400.0,                # mm
    ShellOuterR=4500.0,                # mm
    ShellInnerZ=2400.0,                # mm
    ShellOuterZ=4300.0,                # mm
    OutsideX0=40.0,
    OutsideThickness=2000.0,           # mm

    # ---- algorithm: track model and vertex fit -----------------------------------------------
    FitCurvature=not args.noCurvature,
    UpstreamTolerance=100.0,           # mm
    JacobianStep=1.0,                  # mm
    ConvergenceDistance=0.1,           # mm

    # ---- vertex selection ------------------------------------------------------------------
    MaxPairDCA=100.0,                  # mm
    Region=args.region,
    RegionR=SOLENOID_RADIUS,
    RegionHalfZ=SOLENOID_HALF_Z,
    OutputLevel=INFO,
)

if args.loose:
    for vertexFinder in (helixPairVertexFinder, arcPairVertexFinder):
        vertexFinder.MaxPairDCA = 1.0e9
        vertexFinder.RequireOppositeCharge = False
        vertexFinder.AllowMultipleVertices = True
        vertexFinder.KeepOverlappingPairs = True
        vertexFinder.RequireBetweenIPAndHits = False
        vertexFinder.RequireUpstream = False
        vertexFinder.Region = "any"
    arcPairVertexFinder.MaxChi2 = 1.0e12
    arcPairVertexFinder.MinUpstreamDistance = -1.0e9

if args.debug:
    helixPairVertexFinder.OutputLevel = DEBUG
    arcPairVertexFinder.OutputLevel = DEBUG

algorithms = []
if not args.noHelix:
    algorithms.append(helixPairVertexFinder)
if not args.noArc:
    algorithms.append(arcPairVertexFinder)

ApplicationMgr(
    TopAlg=algorithms,
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc"), UniqueIDGenSvc("UniqueIDGenSvc")],
    OutputLevel=INFO,
)
