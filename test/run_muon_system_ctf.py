# Muon-system tracking (CombinatorialTrackFinder) of the HNL low-coupling sample, used to
# test the helix vertexer on track states. Input/output: environment variables MS_IN, MS_OUT.
#

import os
from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventDataSvc
from Configurables import CombinatorialTrackFinder
from Configurables import GeoSvc

# Geometry service: gives access to the detector geometry
geometryService = GeoSvc("GeoSvc")
geometryService.detectors = ["/eos/user/c/cnaborpa/k4geo/FCCee/IDEA/compact/IDEA_o1_v05/IDEA_o1_v05.xml"]
geometryService.OutputLevel = INFO
geometryService.EnableGeant4Geo = False

# Track finding: CombinatorialTrackFinder on the muon-system hits
combinatorialTrackFinder = CombinatorialTrackFinder("CombinatorialTrackFinder")
combinatorialTrackFinder.DetectorName = "Muon-System"           # Name of the detector to process
combinatorialTrackFinder.ParticleType = "muon"                  # Particle type for material effects
combinatorialTrackFinder.EncodingStringParameterName = "MuonSystemReadoutID"
combinatorialTrackFinder.OutputLevel = INFO

# GenFit configuration
combinatorialTrackFinder.UseGenFit = False                      # Enable GenFit track fitting
combinatorialTrackFinder.MaxFitIterations = 100                 # Maximum iterations for the Kalman fit
# Inner propagation (disabled: solenoid boundary extrapolation not yet calibrated)
combinatorialTrackFinder.DoInnerPropagation = False
combinatorialTrackFinder.InnerPropTargetRadius = 100.0          # target radius in cm if enabled

# ---- N-hit combinatorial strategy ----
combinatorialTrackFinder.MinTrackHits = 3                       # minimum hits to form a track
combinatorialTrackFinder.MaxCombinatorialHits = 8               # cap on k in C(N,k) enumeration
combinatorialTrackFinder.NHitMaxChi2NDF = 10.0                  # max chi2/NDF of N-hit circle fit
combinatorialTrackFinder.HitIsolationCut = 0                    # min cross-region distance (cm); 0 to disable
combinatorialTrackFinder.MinLayerSpan = 3                       # min distinct detector layers per combo
combinatorialTrackFinder.MaxComboPT = 200.0                     # GeV: maximum combo pT (excludes combos above this threshold during search)
combinatorialTrackFinder.OutlierSigma = 5.0                     # outlier removal threshold (σ units)
combinatorialTrackFinder.SigmaHitDefault = 0.04                 # default hit resolution (cm = 0.4 mm)
combinatorialTrackFinder.MaxOutlierIterations = 10              # max rounds of outlier removal
# Road cuts: guard against cross-track hit mixing
combinatorialTrackFinder.MaxConsecDeltaPhi = 0.15               # max |Δφ| between consecutive layers (rad); -1 to disable
combinatorialTrackFinder.MaxComboPhiSpread = 0.5                # max total φ spread across combo (rad); -1 to disable
# Distance cuts: guard against ghost combos spanning different detector regions
combinatorialTrackFinder.MaxConsecutiveHitDistance = -1         # max 3D dist between adjacent (inner->outer) hits (cm); -1 to disable
combinatorialTrackFinder.MaxPairHitDistance = 350.0             # max 3D dist between any two hits in a combo (cm); 300 accommodates barrel->endcap
# Hit quality pre-filter and proximity-based combo selection
combinatorialTrackFinder.MaxHitEdepKeV = 10                      # reject hits with edep > 10 keV (delta-rays / hadronic secondaries)
combinatorialTrackFinder.ProximityScoreWeight = 1.0             # weight for edep-homogeneity in combo ranking (0 to disable)
combinatorialTrackFinder.EdepNormKeV = 3.0                      # normalisation: 2 keV avg deviation ~ +1 chi2 unit (for k>=4 combined score)
combinatorialTrackFinder.MaxHitsPerCompositeID = 0              # 0 = use all hits; >0 caps per layer to top-N by proximity score
# conditions for second track in the events
combinatorialTrackFinder.NeighbourTrackMaxDist = 10.0           # cm (= 50 mm)
combinatorialTrackFinder.NeighbourTrackMinHits = 4              # must have >=4 hits
combinatorialTrackFinder.NeighbourTrackMaxChi2NDF = 3.0         # must have chi2/NDF <= 3
# Guided crowded-layer hit selection
combinatorialTrackFinder.doCrowdedLayerHitSelection = False     # enable guided crowded-layer hit selection

# ---- Collections ----
combinatorialTrackFinder.InputHitCollection = ["MSTrackerHits"]
combinatorialTrackFinder.OutputTrackCollection = ["StandaloneMuonTracks"]
combinatorialTrackFinder.OutputHitCollection = ["StandaloneMuonTrackHits"]
combinatorialTrackFinder.InputRecoSimLinkCollection = ["MSTrackerHitRelations"]
combinatorialTrackFinder.OutputRecoSimLinkCollection = ["StandaloneMuonTrackHitSimLinks"]

# Input/Output service
inputOutputService = IOSvc()
import sys
inputOutputService.Input = os.environ['MS_IN']
inputOutputService.Output = os.environ['MS_OUT']

# Application manager configuration
ApplicationMgr(
    TopAlg=[combinatorialTrackFinder],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc")],
    OutputLevel=INFO
)
