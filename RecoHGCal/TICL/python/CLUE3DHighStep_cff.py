import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.TICLSeedingRegions_cff import ticlSeedingGlobal, ticlSeedingGlobalHFNose
from RecoHGCal.TICL.trackstersProducer_cfi import trackstersProducer as _trackstersProducer
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer as _filteredLayerClustersProducer

# CLUSTER FILTERING/MASKING

filteredLayerClustersCLUE3DHigh = _filteredLayerClustersProducer.clone(
    clusterFilter = "ClusterFilterByAlgoAndSize",
    min_cluster_size = 2, # inclusive
    iteration_label = "CLUE3DHigh"
)

# PATTERN RECOGNITION

ticlTrackstersCLUE3DHigh = _trackstersProducer.clone(
    filtered_mask = "filteredLayerClustersCLUE3DHigh:CLUE3DHigh",
    seeding_regions = "ticlSeedingGlobal",
    itername = "CLUE3DHigh",
    patternRecognitionBy = "CLUE3D",
    pluginPatternRecognitionByCLUE3D = dict (
        criticalDensity = [0.6, 0.6, 0.6],
        criticalEtaPhiDistance = [0.025, 0.025, 0.025],
        kernelDensityFactor = [0.2, 0.2, 0.2],
        algo_verbosity = 0,
        doPidCut = True,
        cutHadProb = 999
    ),
    inferenceAlgo = cms.string('TracksterInferenceByONNX'),
    pluginInferenceAlgoTracksterInferenceByONNX = cms.PSet(
        algo_verbosity = cms.int32(0),
        type = cms.string("TracksterInferenceByONNX"),
        onnxPIDModelPath = cms.string("RecoHGCal/TICL/data/ticlv5/onnx_models/CNN/patternrecognition/id_v0.onnx"),
        onnxEnergyModelPath = cms.string(""),
        inputNames = cms.vstring("input"),
        outputNamesPID = cms.vstring("pid_output"),
        outputNamesEnergy = cms.vstring("enreg_output"),
        outputProbabilityIndices = cms.vuint32(0, 1, 2, 3, 4, 5, 6, 7),
        inputFormat = cms.string("legacyCNN"),
        minTracksterEnergy = cms.double(1.0),
        nLayers = cms.int32(50),
        maxClusters = cms.int32(10),
        doPID = cms.bool(True),
        doRegression = cms.bool(False),
        miniBatchSize = cms.untracked.int32(64),
        reportTiming = cms.untracked.bool(False),
    ),
)

ticlCLUE3DHighStepTask = cms.Task(ticlSeedingGlobal
    ,filteredLayerClustersCLUE3DHigh
    ,ticlTrackstersCLUE3DHigh)
