"""
Full reconstruction for the displaced 2-track LLP scan, ending in the ported
graph-based displaced vertex finder.

  digi (VTX + Si wrapper + DCH)  ->  GGTFTrackFinder  ->  GenfitTrackFitter (DAF)
                                 ->  GraphDisplacedVertexFinder

The fitter uses the Deterministic Annealing Filter (DAF), which down-weights
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

svc = IOSvc("IOSvc")
svc.Input = args.inputFile
svc.Output = args.outputFile

id_service = UniqueIDGenSvc("UniqueIDGenSvc")

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [os.path.join(os.environ["K4GEO"], args.compact)]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

from Configurables import DDPlanarDigi, DCHdigi_v02

# ---------------- vertex detector ----------------
innerVertexResolution_x = 0.003
innerVertexResolution_y = 0.003
innerVertexResolution_t = 1000
outerVertexResolution_x = 0.050 / math.sqrt(12)
outerVertexResolution_y = 0.150 / math.sqrt(12)
outerVertexResolution_t = 1000

vtxb = DDPlanarDigi("VTXBdigitizer")
vtxb.SubDetectorName = "Vertex"
vtxb.IsStrip = False
vtxb.ResolutionU = [innerVertexResolution_x] * 3 + [outerVertexResolution_x] * 2
vtxb.ResolutionV = [innerVertexResolution_y] * 3 + [outerVertexResolution_y] * 2
vtxb.ResolutionT = [innerVertexResolution_t] * 3 + [outerVertexResolution_t] * 2
vtxb.SimTrackHitCollectionName = ["VertexBarrelCollection"]
vtxb.SimTrkHitRelCollection = ["VTXBSimDigiLinks"]
vtxb.TrackerHitCollectionName = ["VTXBDigis"]
vtxb.ForceHitsOntoSurface = True
vtxb.CellIDBits = 32

vtxd = DDPlanarDigi("VTXDdigitizer")
vtxd.SubDetectorName = "Vertex"
vtxd.IsStrip = False
vtxd.ResolutionU = [outerVertexResolution_x] * 3
vtxd.ResolutionV = [outerVertexResolution_y] * 3
vtxd.ResolutionT = [outerVertexResolution_t] * 3
vtxd.SimTrackHitCollectionName = ["VertexEndcapCollection"]
vtxd.SimTrkHitRelCollection = ["VTXDSimDigiLinks"]
vtxd.TrackerHitCollectionName = ["VTXDDigis"]
vtxd.ForceHitsOntoSurface = True
vtxd.CellIDBits = 32

# ---------------- silicon wrapper ----------------
siWrapperResolution_x = 0.050 / math.sqrt(12)
siWrapperResolution_y = 1.0 / math.sqrt(12)
siWrapperResolution_t = 0.040

siwrb = DDPlanarDigi("SiWrBdigitizer")
siwrb.SubDetectorName = "SiWrB"
siwrb.IsStrip = False
siwrb.ResolutionU = [siWrapperResolution_x] * 2
siwrb.ResolutionV = [siWrapperResolution_y] * 2
siwrb.ResolutionT = [siWrapperResolution_t] * 2
siwrb.SimTrackHitCollectionName = ["SiWrBCollection"]
siwrb.SimTrkHitRelCollection = ["SiWrBSimDigiLinks"]
siwrb.TrackerHitCollectionName = ["SiWrBDigis"]
siwrb.ForceHitsOntoSurface = True
siwrb.CellIDBits = 32

siwrd = DDPlanarDigi("SiWrDdigitizer")
siwrd.SubDetectorName = "SiWrD"
siwrd.IsStrip = False
siwrd.ResolutionU = [siWrapperResolution_x] * 2
siwrd.ResolutionV = [siWrapperResolution_y] * 2
siwrd.ResolutionT = [siWrapperResolution_t] * 2
siwrd.SimTrackHitCollectionName = ["SiWrDCollection"]
siwrd.SimTrkHitRelCollection = ["SiWrDSimDigiLinks"]
siwrd.TrackerHitCollectionName = ["SiWrDDigis"]
siwrd.ForceHitsOntoSurface = True
siwrd.CellIDBits = 32

# ---------------- drift chamber ----------------
dch = DCHdigi_v02(
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

finder = GGTFTrackFinder(
    "GGTFTrackFinder",
    InputPlanarHitCollections=["VTXBDigis", "VTXDDigis", "SiWrBDigis", "SiWrDDigis"],
    InputWireHitCollections=["DCH_DigiCollection"],
    OutputTracksGGTF=["GGTFTracks"],
    ModelPath=args.modelPath,
    OutputLevel=INFO,
)

# ---------------- track fitting (DAF) ----------------
from Configurables import GenfitTrackFitter

fitter = GenfitTrackFitter(
    "GenfitTrackFitter",
    InputTracks=["GGTFTracks"],
    OutputFittedTracks=["FittedTracks"],
    OutputFittedTracksWithFilteredHits=["FittedTracksWithFilteredHits"],
    OutputFittedHits=["FittedHits"],
    OutputLevel=INFO,
)
fitter.FitterType = args.fitter        # "DAF": deterministic annealing filter
fitter.BetaInit = 100.0                # annealing schedule (DAF only)
fitter.BetaFinal = 0.05
fitter.BetaSteps = 15

# ---------------- graph-based displaced vertexing ----------------
from Configurables import GraphDisplacedVertexFinder

vertexer = GraphDisplacedVertexFinder(
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
    TopAlg=[vtxb, vtxd, siwrb, siwrd, dch, finder, fitter, vertexer],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc"), id_service],
    OutputLevel=INFO,
)
