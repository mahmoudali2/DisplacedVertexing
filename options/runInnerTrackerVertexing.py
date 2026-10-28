"""
Run the 3D helix pair vertexer on inner-tracker (silicon) tracks.

Settings: uniform 2 T, no energy loss, helix from the 8 innermost hits,
opposite charge, DCA < 2 mm, cos(opening) > -0.95, vertex r_xy > 10 mm.

Author: Mahmoud Althakeel

Example
-------
    k4run runInnerTrackerVertexing.py \
        --inputFile  it_tracks.root \
        --outputFile it_vertices.root

Output collection: InnerTrackerHelixVertices (algorithmType 4).
"""

from Gaudi.Configuration import INFO
from Configurables import EventDataSvc, UniqueIDGenSvc
from k4FWCore import IOSvc, ApplicationMgr
from k4FWCore.parseArgs import parser

parser.add_argument("--inputFile", required=True)
parser.add_argument("--outputFile", required=True)
parser.add_argument("--inputTracks", default="GlobalTrackerTracks")
parser.add_argument("--field", type=float, default=2.0, help="solenoid field [T]")
parser.add_argument("--loose", action="store_true", help="no selection, keep all pairs (validation, cut scans)")
args = parser.parse_args()

svc = IOSvc("IOSvc")
svc.Input = args.inputFile
svc.Output = args.outputFile

try:
    from Configurables import HelixPairVertexFinder
except ImportError:
    from DisplacedVertexing.DisplacedVertexingConf import HelixPairVertexFinder

helix = HelixPairVertexFinder(
    "HelixPairVertexFinder",
    InputTracks=[args.inputTracks],
    OutputVertices=["InnerTrackerHelixVertices"],
    # uniform solenoid field (FieldR / FieldHalfZ left at "everywhere")
    InnerField=args.field,
    # no energy loss in the tracker (dEdxShell = dEdxOutside = 0 by default)
    ParticleMass=0.1056584,   # muons
    # silicon scale of the track model and propagation
    MinHits=3,
    MaxHitsInFit=8,           # the innermost hits describe the track near the decay best
    MinXYSpan=5.0,            # mm
    Step=5.0,                 # mm (chord error < 0.01 mm for pT > 1 GeV)
    BackFactor=1.3,
    UpstreamTolerance=300.0,  # mm, as the s range [-2200, 300] mm of the Python analysis
    # selection
    RequireUpstream=False,
    RequireBetweenIPAndHits=False,
    RequireOppositeCharge=True,
    MaxPairDCA=2.0,           # mm
    MinCosOpening=-0.95,
    Region="barrel",          # vertex r_xy > RegionR
    RegionR=10.0,             # mm
    RegionHalfZ=1.0e12,
    OutputLevel=INFO,
)

if args.loose:
    helix.MaxPairDCA = 1.0e9
    helix.RequireOppositeCharge = False
    helix.AllowMultipleVertices = True
    helix.KeepOverlappingPairs = True
    helix.MinCosOpening = -1.0
    helix.Region = "any"

ApplicationMgr(
    TopAlg=[helix],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc"), UniqueIDGenSvc("UniqueIDGenSvc")],
    OutputLevel=INFO,
)
