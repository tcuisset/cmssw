import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.TICLSeedingRegions_cff import ticlSeedingGlobal, ticlSeedingGlobalHFNose
from RecoHGCal.TICL.trackstersProducer_cfi import trackstersProducer as _trackstersProducer
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer as _filteredLayerClustersProducer

# CLUSTER FILTERING/MASKING

filteredLayerClustersRecovery = _filteredLayerClustersProducer.clone(
    clusterFilter = "ClusterFilterByAlgoAndSize",
    min_cluster_size = 2, # inclusive
    iteration_label = "Recovery",
    LayerClustersInputMask = 'ticlTrackstersCLUE3DHigh',
    algo_number = [6, 7, 8],
)

# PATTERN RECOGNITION

ticlTrackstersRecovery = _trackstersProducer.clone(
    filtered_mask = "filteredLayerClustersRecovery:Recovery",
    original_mask = 'ticlTrackstersCLUE3DHigh',
    seeding_regions = "ticlSeedingGlobal",
    itername = "Recovery",
    patternRecognitionBy = "Recovery",
    inferenceAlgo=cms.string(''),
    pluginPatternRecognitionByRecovery = dict (
        algo_verbosity = 0
    ),
    pluginInferenceAlgoTracksterInferenceByONNX = cms.PSet(
      algo_verbosity = cms.int32(0),
      onnxPIDModelPath = cms.string('RecoHGCal/TICL/data/ticlv5/onnx_models/PFN/patternrecognition/id_v0.onnx'),
      onnxEnergyModelPath = cms.string('RecoHGCal/TICL/data/ticlv5/onnx_models/PFN/patternrecognition/energy_v0.onnx'),
      inputNames = cms.vstring(
        'input',
        'input_tr_features'
      ),
      outputNamesEnergy = cms.vstring('enreg_output'),
      outputNamesPID = cms.vstring('pid_output'),
      outputProbabilityIndices = cms.vuint32(0, 1, 2, 3, 4, 5, 6, 7),
      inputFormat = cms.string('legacyPFN'),
      minTracksterEnergy = cms.double(1),
      nLayers = cms.int32(50),
      maxClusters = cms.int32(10),
      doPID = cms.bool(False),
      doRegression = cms.bool(False),
      type = cms.string('TracksterInferenceByONNX')
    ),
)


ticlRecoveryStepTask = cms.Task(ticlSeedingGlobal
    ,filteredLayerClustersRecovery
    ,ticlTrackstersRecovery)
