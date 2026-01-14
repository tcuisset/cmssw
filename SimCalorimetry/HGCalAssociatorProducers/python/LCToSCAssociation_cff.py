import FWCore.ParameterSet.Config as cms

from SimCalorimetry.HGCalAssociatorProducers.LCToSCAssociatorEDProducer_cfi import LCToSCAssociatorEDProducer as _LCToSCAssociatorEDProducer

layerClusterSimClusterAssociationProducer = _LCToSCAssociatorEDProducer.clone(
    label_scl = cms.InputTag("mix","MergedCaloTruth"),
)
layerClusterBoundaryTrackSimClusterAssociationProducer = _LCToSCAssociatorEDProducer.clone(
    label_scl = cms.InputTag("mix","MergedCaloTruthBoundaryTrackSimCluster"),
)
layerClusterMergedSimClusterAssociationProducer = _LCToSCAssociatorEDProducer.clone(
    label_scl = cms.InputTag("mix","MergedCaloTruth"),
)
# the next associator is an associator of LCs->SimCluster dataforamt but using SimCluster collection that is a 1-1 mapping to CaloParticle. 
# this way downstream code only has one dataformat (SimCluster) instead of 2 (CaloParticle & SimCluster)
# at some point layerClusterCaloParticleAssociationProducer will be removed, keeping only layerClusterBoundaryTrackSimClusterAssociationProducer (once downstream code is updated)
layerClusterCaloParticleSimClusterAssociationProducer = _LCToSCAssociatorEDProducer.clone(
    label_scl = cms.InputTag("mix","MergedCaloTruth"),
)

barrelLayerClusterSimClusterAssociation = cms.EDProducer("LCToSCAssociatorEDProducer",
    associator = cms.InputTag('barrelLCToSCAssociatorByEnergyScoreProducer'),
    label_scl = cms.InputTag("mix","MergedCaloTruth"),
    label_lcl = cms.InputTag("hgcalMergeLayerClusters")
)

layerClusterSimClusterAssociationHFNose = layerClusterSimClusterAssociation.clone(
    label_lcl = "hgcalLayerClustersHFNose"
)

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
for assoc_ in [layerClusterSimClusterAssociationProducer, layerClusterBoundaryTrackSimClusterAssociationProducer, layerClusterMergedSimClusterAssociationProducer, layerClusterCaloParticleSimClusterAssociationProducer, barrelLayerClusterSimClusterAssociation, layerClusterSimClusterAssociationHFNose]:
    premix_stage2.toModify(assoc_,
        label_scl = cms.InputTag("mixData", assoc_.label_scl.productInstanceLabel)
    )


