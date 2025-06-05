import FWCore.ParameterSet.Config as cms

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2

from RecoHGCal.TICL.simTrackstersProducer_cfi import simTrackstersProducer as _simTrackstersProducer
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer as _filteredLayerClustersProducer
# from RecoHGCal.TICL.simTrackstersMerger_cfi import simTrackstersMerger as _simTrackstersMerger

# CA - PATTERN RECOGNITION


filteredLayerClustersSimTracksters = _filteredLayerClustersProducer.clone(
    clusterFilter = "ClusterFilterByAlgoAndSize",
    min_cluster_size = 0, # inclusive
    iteration_label = "ticlSimTracksters"
)

from RecoHGCal.TICL.simClusterMerger_cfi import simClusterMerger
premix_stage2.toModify(simClusterMerger,
    simclusters = "mixData:MergedCaloTruth",
    caloparticles = "mixData:MergedCaloTruth"
)

ticlSimTracksters = _simTrackstersProducer.clone(
    computeLocalTime = cms.bool(False),
    simClusterCollections = cms.VPSet(
        cms.PSet(
            simCluster = cms.InputTag("mix", "MergedCaloTruth"),
            simClusterToCaloParticlesMap = cms.InputTag("SimClusterToCaloParticleAssociation"),
            layerClusterToSimClusterMap = cms.InputTag("layerClusterSimClusterAssociationProducer"),
            outputProductLabel = cms.string(""),
            useForSimTICLCandidate = cms.bool(True)
        ),
        cms.PSet(
            simCluster = cms.InputTag("simClusterMerger"),
            simClusterToCaloParticlesMap = cms.InputTag("simClusterMerger"),
            layerClusterToSimClusterMap = cms.InputTag("layerClusterSimClusterMergedAssociationProducer"),
            outputProductLabel = cms.string("simClusterMerger"),
            useForSimTICLCandidate = cms.bool(False)
        ),
    )
)
from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5
ticl_v5.toModify(ticlSimTracksters,computeLocalTime = cms.bool(True))

premix_stage2.toModify(ticlSimTracksters,
    #simclusters = "mixData:MergedCaloTruth",
    simClusterCollections={
        0 : dict(simCluster="mixData:MergedCaloTruth")
    },
    caloparticles = "mixData:MergedCaloTruth",
)

# ticlSimTrackstersMerged = _simTrackstersMerger.clone()
# premix_stage2.toModify(_simTrackstersMerger,
#     simclusters = "mixData:MergedCaloTruth",
#     caloparticles = "mixData:MergedCaloTruth",
# )

ticlSimTrackstersTask = cms.Task(filteredLayerClustersSimTracksters, simClusterMerger, ticlSimTracksters,) # ticlSimTrackstersMerged

simTracksterCollections_inputTags = [
    cms.InputTag("ticlSimTracksters"), # SimTracksters built from SimCluster (ie one for every SimTrack crossing boundary)
    cms.InputTag("ticlSimTracksters", "fromCPs"), # SimTrackster built from CaloParticles (ie one for every GenParticle have any descendant SimTrack crossing boundary)
    cms.InputTag("ticlSimTrackstersMerged") # SimTrackster built from clusters of SimTrack crossing boundary, merging very close- : trying to approach what can reasonably be reconstructed given HGCAL granularity
]
#simTrackstersCollections = {"simTrackstersFromSC" : cms.InputTag("ticlSimTracksters"), "simTrackstersFromCP" : cms.InputTag("ticlSimTracksters", "fromCPs"), "simTrackstersMerged" : cms.InputTag("ticlSimTrackstersMerged")}
