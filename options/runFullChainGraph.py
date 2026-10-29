"""
Full reconstruction for the displaced 2-track LLP scan, ending in the ported
graph-based displaced vertex ggtfTrackFinder.

  digi (VTX + Si wrapper + DCH)  ->  GGTFTrackFinder  ->  GenfitTrackFitter (DAF)
                                 ->  GraphDisplacedVertexFinder

The genfitTrackFitter uses the Deterministic Annealing Filter (DAF), which down-weights
outlier hits during the fit - the appropriate choice for displaced tracks where
hit-to-track association is less certain than for prompt tracks.

Vertexing input is the FITTED track collection.
"""
import os
import math
from Gaudi.Configuration import INFO
from Configurables import EventDataSvc, GeoSvc, UniqueIDGenSvc
from k4FWCore import IOSvc, ApplicationMgr
from k4FWCore.parseArgs import parser

parser.add_argument("--inputFile", required=True)
parser.add_argument("--outputFile", required=True)
parser.add_argument("--modelPath", required=True, help="ONNX model for GGTF")
parser.add_argument("--fitter", default="DAF", help="DAF (default) or KALMAN_REF")
parser.add_argument("--compact", default="FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml")
args = parser.parse_args()

inputOutputService = IOSvc("IOSvc")
inputOutputService.Input = args.inputFile
inputOutputService.Output = args.outputFile

uniqueIdService = UniqueIDGenSvc("UniqueIDGenSvc")

geometryService = GeoSvc("GeoSvc")
geometryService.detectors = [os.path.join(os.environ["K4GEO"], args.compact)]
geometryService.OutputLevel = INFO
geometryService.EnableGeant4Geo = False

from Configurables import DDPlanarDigi, DCHdigi_v02

# ---------------- vertex detector ----------------
innerVertexResolution_x = 0.003
innerVertexResolution_y = 0.003
innerVertexResolution_t = 1000
outerVertexResolution_x = 0.050 / math.sqrt(12)
outerVertexResolution_y = 0.150 / math.sqrt(12)
outerVertexResolution_t = 1000

vertexBarrelDigitizer = DDPlanarDigi("VTXBdigitizer")
vertexBarrelDigitizer.SubDetectorName = "Vertex"
vertexBarrelDigitizer.IsStrip = False
vertexBarrelDigitizer.ResolutionU = [innerVertexResolution_x] * 3 + [outerVertexResolution_x] * 2
vertexBarrelDigitizer.ResolutionV = [innerVertexResolution_y] * 3 + [outerVertexResolution_y] * 2
vertexBarrelDigitizer.ResolutionT = [innerVertexResolution_t] * 3 + [outerVertexResolution_t] * 2
vertexBarrelDigitizer.SimTrackHitCollectionName = ["VertexBarrelCollection"]
vertexBarrelDigitizer.SimTrkHitRelCollection = ["VTXBSimDigiLinks"]
vertexBarrelDigitizer.TrackerHitCollectionName = ["VTXBDigis"]
vertexBarrelDigitizer.ForceHitsOntoSurface = True
vertexBarrelDigitizer.CellIDBits = 32

vertexEndcapDigitizer = DDPlanarDigi("VTXDdigitizer")
vertexEndcapDigitizer.SubDetectorName = "Vertex"
vertexEndcapDigitizer.IsStrip = False
vertexEndcapDigitizer.ResolutionU = [outerVertexResolution_x] * 3
vertexEndcapDigitizer.ResolutionV = [outerVertexResolution_y] * 3
vertexEndcapDigitizer.ResolutionT = [outerVertexResolution_t] * 3
vertexEndcapDigitizer.SimTrackHitCollectionName = ["VertexEndcapCollection"]
vertexEndcapDigitizer.SimTrkHitRelCollection = ["VTXDSimDigiLinks"]
vertexEndcapDigitizer.TrackerHitCollectionName = ["VTXDDigis"]
vertexEndcapDigitizer.ForceHitsOntoSurface = True
vertexEndcapDigitizer.CellIDBits = 32

# ---------------- silicon wrapper ----------------
siWrapperResolution_x = 0.050 / math.sqrt(12)
siWrapperResolution_y = 1.0 / math.sqrt(12)
siWrapperResolution_t = 0.040

siliconWrapperBarrelDigitizer = DDPlanarDigi("SiWrBdigitizer")
siliconWrapperBarrelDigitizer.SubDetectorName = "SiWrB"
siliconWrapperBarrelDigitizer.IsStrip = False
siliconWrapperBarrelDigitizer.ResolutionU = [siWrapperResolution_x] * 2
siliconWrapperBarrelDigitizer.ResolutionV = [siWrapperResolution_y] * 2
siliconWrapperBarrelDigitizer.ResolutionT = [siWrapperResolution_t] * 2
siliconWrapperBarrelDigitizer.SimTrackHitCollectionName = ["SiWrBCollection"]
siliconWrapperBarrelDigitizer.SimTrkHitRelCollection = ["SiWrBSimDigiLinks"]
siliconWrapperBarrelDigitizer.TrackerHitCollectionName = ["SiWrBDigis"]
siliconWrapperBarrelDigitizer.ForceHitsOntoSurface = True
siliconWrapperBarrelDigitizer.CellIDBits = 32

siliconWrapperEndcapDigitizer = DDPlanarDigi("SiWrDdigitizer")
siliconWrapperEndcapDigitizer.SubDetectorName = "SiWrD"
siliconWrapperEndcapDigitizer.IsStrip = False
siliconWrapperEndcapDigitizer.ResolutionU = [siWrapperResolution_x] * 2
siliconWrapperEndcapDigitizer.ResolutionV = [siWrapperResolution_y] * 2
siliconWrapperEndcapDigitizer.ResolutionT = [siWrapperResolution_t] * 2
siliconWrapperEndcapDigitizer.SimTrackHitCollectionName = ["SiWrDCollection"]
siliconWrapperEndcapDigitizer.SimTrkHitRelCollection = ["SiWrDSimDigiLinks"]
siliconWrapperEndcapDigitizer.TrackerHitCollectionName = ["SiWrDDigis"]
siliconWrapperEndcapDigitizer.ForceHitsOntoSurface = True
siliconWrapperEndcapDigitizer.CellIDBits = 32

# ---------------- drift chamber ----------------
driftChamberDigitizer = DCHdigi_v02(
    "DCHdigi2",
    InputSimHitCollection=["DCHCollection"],
    OutputDigihitCollection=["DCH_DigiCollection"],
    OutputLinkCollection=["DCH_DigiSimAssociationCollection"],
    DCH_name="DCH_v2",
    zResolution_mm=30.0,
    xyResolution_mm=0.1,
    Deadtime_ns=400.0,
    GasType=0,
    ReadoutWindowStartTime_ns=1.0,
    ReadoutWindowDuration_ns=450.0,
    DriftVelocity_um_per_ns=-1.0,
    SignalVelocity_mm_per_ns=200.0,
    OutputLevel=INFO,
)

# ---------------- track finding ----------------
from Configurables import GGTFTrackFinder

ggtfTrackFinder = GGTFTrackFinder(
    "GGTFTrackFinder",
    InputPlanarHitCollections=["VTXBDigis", "VTXDDigis", "SiWrBDigis", "SiWrDDigis"],
    InputWireHitCollections=["DCH_DigiCollection"],
    OutputTracksGGTF=["GGTFTracks"],
    ModelPath=args.modelPath,
    OutputLevel=INFO,
)

# ---------------- track fitting (DAF) ----------------
from Configurables import GenfitTrackFitter

genfitTrackFitter = GenfitTrackFitter(
    "GenfitTrackFitter",
    InputTracks=["GGTFTracks"],
    OutputFittedTracks=["FittedTracks"],
    OutputFittedTracksWithFilteredHits=["FittedTracksWithFilteredHits"],
    OutputFittedHits=["FittedHits"],
    OutputLevel=INFO,
)
genfitTrackFitter.FitterType = args.genfitTrackFitter        # "DAF": deterministic annealing filter
genfitTrackFitter.BetaInit = 100.0                # annealing schedule (DAF only)
genfitTrackFitter.BetaFinal = 0.05
genfitTrackFitter.BetaSteps = 15

# ---------------- graph-based displaced vertexing ----------------
from Configurables import GraphDisplacedVertexFinder

graphVertexFinder = GraphDisplacedVertexFinder(
    "GraphDisplacedVertexFinder",
    InputTracks=["FittedTracks"],
    OutputDisplacedVertices=["DisplacedVertices"],
    OutputPrimaryVertex=["GraphPrimaryVertex"],
    MagneticFieldZ=2.0,
    # two-track LLP decays (X -> mu mu)
    MinTracks=2,
    MinSupport=0,        # no triplet support exists for a pair
    kNN=-1,
    UseSelfSeeding=True, # EuclideanSeed is the seeding path that works for N=2
    SkipFailedFitTracks=True,
    OutputLevel=INFO,
)

ApplicationMgr(
    TopAlg=[vertexBarrelDigitizer, vertexEndcapDigitizer, siliconWrapperBarrelDigitizer, siliconWrapperEndcapDigitizer, driftChamberDigitizer, ggtfTrackFinder, genfitTrackFitter, graphVertexFinder],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc"), uniqueIdService],
    OutputLevel=INFO,
)
