import FWCore.ParameterSet.Config as cms

tracksterSimTracksterAssociationLinking = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
    associator = cms.InputTag('simTracksterHitLCAssociatorByEnergyScoreProducer'),
    label_tst = cms.InputTag("ticlTrackstersMerge"),
    label_simTst = cms.InputTag("ticlSimTracksters", "fromCPs"),
    label_lcl = cms.InputTag("hgcalMergeLayerClusters"),
    label_scl = cms.InputTag("mix", "MergedCaloTruth"),
    label_cp = cms.InputTag("mix","MergedCaloTruth"),
)

tracksterSimTracksterAssociationPR = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
    associator = cms.InputTag('simTracksterHitLCAssociatorByEnergyScoreProducer'),
    label_tst = cms.InputTag("ticlTrackstersMerge"),
    label_simTst = cms.InputTag("ticlSimTracksters"),
    label_lcl = cms.InputTag("hgcalMergeLayerClusters"),
    label_scl = cms.InputTag("mix", "MergedCaloTruth"),
    label_cp = cms.InputTag("mix","MergedCaloTruth"),
)


tracksterSimTracksterAssociationLinkingbyCLUE3D = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
    associator = cms.InputTag('simTracksterHitLCAssociatorByEnergyScoreProducer'),
    label_tst = cms.InputTag("ticlTrackstersCLUE3DHigh"),
    label_simTst = cms.InputTag("ticlSimTracksters", "fromCPs"),
    label_lcl = cms.InputTag("hgcalMergeLayerClusters"),
    label_scl = cms.InputTag("mix", "MergedCaloTruth"),
    label_cp = cms.InputTag("mix","MergedCaloTruth"),
)

tracksterSimTracksterAssociationPRbyCLUE3D = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
    associator = cms.InputTag('simTracksterHitLCAssociatorByEnergyScoreProducer'),
    label_tst = cms.InputTag("ticlTrackstersCLUE3DHigh"),
    label_simTst = cms.InputTag("ticlSimTracksters"),
    label_lcl = cms.InputTag("hgcalMergeLayerClusters"),
    label_scl = cms.InputTag("mix", "MergedCaloTruth"),
    label_cp = cms.InputTag("mix","MergedCaloTruth"),
)

# Linking CaloParticleSimTrackster to CLUE3D EM tracksters
tracksterSimTracksterAssociationLinkingbyCLUE3DEM = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
    associator = cms.InputTag('simTracksterHitLCAssociatorByEnergyScoreProducer'),
    label_tst = cms.InputTag("ticlTrackstersCLUE3DEM"),
    label_simTst = cms.InputTag("ticlSimTracksters", "fromCPs"),
    label_lcl = cms.InputTag("hgcalMergeLayerClusters"),
    label_scl = cms.InputTag("mix", "MergedCaloTruth"),
    label_cp = cms.InputTag("mix","MergedCaloTruth"),
)

tracksterSimTracksterAssociationPRbyCLUE3DEM = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
    associator = cms.InputTag('simTracksterHitLCAssociatorByEnergyScoreProducer'),
    label_tst = cms.InputTag("ticlTrackstersCLUE3DEM"),
    label_simTst = cms.InputTag("ticlSimTracksters"),
    label_lcl = cms.InputTag("hgcalMergeLayerClusters"),
    label_scl = cms.InputTag("mix", "MergedCaloTruth"),
    label_cp = cms.InputTag("mix","MergedCaloTruth"),
)

# Linking CaloParticleSimTrackster to EM superclusters
tracksterSimTracksterAssociationLinkingSuperclustering = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
    associator = cms.InputTag('simTracksterHitLCAssociatorByEnergyScoreProducer'),
    label_tst = cms.InputTag("ticlTracksterLinksSuperclustering"),
    label_simTst = cms.InputTag("ticlSimTracksters", "fromCPs"),
    label_lcl = cms.InputTag("hgcalMergeLayerClusters"),
    label_scl = cms.InputTag("mix", "MergedCaloTruth"),
    label_cp = cms.InputTag("mix","MergedCaloTruth"),
)

tracksterSimTracksterAssociationPRSuperclustering = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
    associator = cms.InputTag('simTracksterHitLCAssociatorByEnergyScoreProducer'),
    label_tst = cms.InputTag("ticlTracksterLinksSuperclustering"),
    label_simTst = cms.InputTag("ticlSimTracksters"),
    label_lcl = cms.InputTag("hgcalMergeLayerClusters"),
    label_scl = cms.InputTag("mix", "MergedCaloTruth"),
    label_cp = cms.InputTag("mix","MergedCaloTruth"),
)

tracksterSimTracksterAssociationLinkingPU = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
    associator = cms.InputTag('simTracksterHitLCAssociatorByEnergyScoreProducer'),
    label_tst = cms.InputTag("ticlTrackstersMerge"),
    label_simTst = cms.InputTag("ticlSimTracksters", "PU"),
    label_lcl = cms.InputTag("hgcalMergeLayerClusters"),
    label_scl = cms.InputTag("mix", "MergedCaloTruth"),
    label_cp = cms.InputTag("mix","MergedCaloTruth"),
)

tracksterSimTracksterAssociationPRPU = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
    associator = cms.InputTag('simTracksterHitLCAssociatorByEnergyScoreProducer'),
    label_tst = cms.InputTag("ticlTrackstersMerge"),
    label_simTst = cms.InputTag("ticlSimTracksters", "PU"),
    label_lcl = cms.InputTag("hgcalMergeLayerClusters"),
    label_scl = cms.InputTag("mix", "MergedCaloTruth"),
    label_cp = cms.InputTag("mix","MergedCaloTruth"),
)


