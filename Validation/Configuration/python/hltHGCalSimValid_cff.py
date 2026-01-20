import FWCore.ParameterSet.Config as cms

from SimCalorimetry.HGCalAssociatorProducers.LCToSCAssociation_cff import scAssocByEnergyScoreProducer as _scAssocByEnergyScoreProducer, layerClusterSimClusterAssociationProducer as _layerClusterSimClusterAssociationProducer, layerClusterBoundaryTrackSimClusterAssociationProducer as _layerClusterBoundaryTrackSimClusterAssociationProducer, layerClusterCaloParticleSimClusterAssociationProducer as _layerClusterCaloParticleSimClusterAssociationProducer

from SimCalorimetry.HGCalAssociatorProducers.SimClusterToCaloParticleAssociation_cfi import SimClusterToCaloParticleAssociation
from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import  allTrackstersToSimTrackstersAssociationsByLCs as _allTrackstersToSimTrackstersAssociationsByLCs

from Validation.HGCalValidation.HLT_TICLIterLabels_cff import hltTiclIterLabels as _hltTiclIterLabels

from RecoLocalCalo.HGCalRecProducers.recHitMapProducer_cff import recHitMapProducer as _recHitMapProducer

hits = ["hltHGCalRecHit:HGCEERecHits", "hltHGCalRecHit:HGCHEFRecHits", "hltHGCalRecHit:HGCHEBRecHits"]
hltRecHitMapProducer = _recHitMapProducer.clone(
    hits = hits,
    hgcalOnly = True,
)

hltScAssocByEnergyScoreProducer = _scAssocByEnergyScoreProducer.clone(
    hits = cms.InputTag("hltRecHitMapProducer", "RefProdVectorHGCRecHitCollection"),
    hitMapTag = cms.InputTag("hltRecHitMapProducer","hgcalRecHitMap"),
)

# Layer Cluster <-> SimCluster/CaloParticle associations
hltLayerClusterSimClusterAssociationProducer = _layerClusterSimClusterAssociationProducer.clone(
    associator = cms.InputTag("hltScAssocByEnergyScoreProducer"),
    label_lcl = cms.InputTag("hltMergeLayerClusters")
)
hltLayerClusterBoundaryTrackSimClusterAssociationProducer = _layerClusterBoundaryTrackSimClusterAssociationProducer.clone(
    associator = cms.InputTag("hltScAssocByEnergyScoreProducer"),
    label_lcl = cms.InputTag("hltMergeLayerClusters")
)
# the next associator is an associator of LCs->SimCluster dataformat but using SimCluster collection that is a 1-1 mapping to CaloParticle. 
# this way downstream code only has one dataformat (SimCluster) instead of 2 (CaloParticle & SimCluster)
hltLayerClusterCaloParticleSimClusterAssociationProducer = _layerClusterCaloParticleSimClusterAssociationProducer.clone(
    associator = cms.InputTag("hltScAssocByEnergyScoreProducer"),
    label_lcl = cms.InputTag("hltMergeLayerClusters")
)

from SimCalorimetry.HGCalAssociatorProducers.AllLayerClusterToTracksterAssociatorsProducer_cfi import AllLayerClusterToTracksterAssociatorsProducer as _AllLayerClusterToTracksterAssociatorsProducer

hltAllLayerClusterToTracksterAssociations = _AllLayerClusterToTracksterAssociatorsProducer.clone(
    layer_clusters = cms.InputTag("hltMergeLayerClusters"),
    tracksterCollections = cms.VInputTag(
        *[cms.InputTag(label) for label in _hltTiclIterLabels],
        cms.InputTag("hltTiclSimTracksters", "fromLegacySimCluster"),
        cms.InputTag("hltTiclSimTracksters", "fromBoundarySimCluster"),
        cms.InputTag("hltTiclSimTracksters", "fromCaloParticle"),
    )
)

hltAllTrackstersToSimTrackstersAssociationsByLCs = _allTrackstersToSimTrackstersAssociationsByLCs.clone(
    allLCtoTSAccoc =  cms.string("hltAllLayerClusterToTracksterAssociations"),
    layerClusters = cms.InputTag("hltMergeLayerClusters"),
    tracksterCollections = cms.VInputTag(
        *[cms.InputTag(label) for label in _hltTiclIterLabels]
    ),
    simTracksterCollections = cms.VInputTag(
        cms.InputTag("hltTiclSimTracksters", "fromLegacySimCluster"),
        cms.InputTag("hltTiclSimTracksters", "fromBoundarySimCluster"),
        cms.InputTag("hltTiclSimTracksters", "fromCaloParticle"),
    ),
)

from SimCalorimetry.HGCalAssociatorProducers.hitToSimClusterCaloParticleAssociator_cfi import hitToSimClusterCaloParticleAssociator as _hitToSimClusterCaloParticleAssociator
hltHitToLegacySimClusterAssociator = _hitToSimClusterCaloParticleAssociator.clone(
    simClusters = cms.InputTag("mix", "MergedCaloTruth"),
    hitMap = cms.InputTag("hltRecHitMapProducer","hgcalRecHitMap"),
    hits = cms.InputTag("hltRecHitMapProducer", "RefProdVectorHGCRecHitCollection"),
)
hltHitToBoundarySimClusterAssociator = _hitToSimClusterCaloParticleAssociator.clone(
    simClusters = cms.InputTag("mix", "MergedCaloTruthBoundaryTrackSimCluster"),
    hitMap = cms.InputTag("hltRecHitMapProducer","hgcalRecHitMap"),
    hits = cms.InputTag("hltRecHitMapProducer", "RefProdVectorHGCRecHitCollection"),
)
hltHitToCPSimClusterAssociator = _hitToSimClusterCaloParticleAssociator.clone(
    simClusters = cms.InputTag("mix", "MergedCaloTruthCaloParticle"), # CaloParticle but in SimCluster dataformat
    hitMap = cms.InputTag("hltRecHitMapProducer","hgcalRecHitMap"),
    hits = cms.InputTag("hltRecHitMapProducer", "RefProdVectorHGCRecHitCollection"),
)



from SimCalorimetry.HGCalAssociatorProducers.AllHitToTracksterAssociatorsProducer_cfi import AllHitToTracksterAssociatorsProducer as _AllHitToTracksterAssociatorsProducer

hltAllHitToTracksterAssociations =  _AllHitToTracksterAssociatorsProducer.clone(
    hitMapTag = cms.InputTag("hltRecHitMapProducer","hgcalRecHitMap"),
    hits = cms.InputTag("hltRecHitMapProducer", "RefProdVectorHGCRecHitCollection"),
    layerClusters = cms.InputTag("hltMergeLayerClusters"),
    tracksterCollections = cms.VInputTag(
        *[cms.InputTag(label) for label in _hltTiclIterLabels],
        cms.InputTag("hltTiclSimTracksters", "fromLegacySimCluster"),
        cms.InputTag("hltTiclSimTracksters", "fromBoundarySimCluster"),
        cms.InputTag("hltTiclSimTracksters", "fromCaloParticle"),
    )
)

from SimCalorimetry.HGCalAssociatorProducers.AllTracksterToSimTracksterAssociatorsByHitsProducer_cfi import AllTracksterToSimTracksterAssociatorsByHitsProducer as _AllTracksterToSimTracksterAssociatorsByHitsProducer
hltAllTrackstersToSimTrackstersAssociationsByHits = _AllTracksterToSimTracksterAssociatorsByHitsProducer.clone(
    allHitToTSAccoc = cms.string("hltAllHitToTracksterAssociations"),
    hits = cms.InputTag("hltRecHitMapProducer", "RefProdVectorHGCRecHitCollection"),
    tracksterCollections = cms.VInputTag(
        *[cms.InputTag(label) for label in _hltTiclIterLabels]
    ),
    simTracksters = cms.VPSet(
        cms.PSet(
            simTracksterCollection=cms.InputTag("hltTiclSimTracksters", "fromLegacySimCluster"),
            hitToSimClusterMap=cms.InputTag("hltHitToLegacySimClusterAssociator")
        ),
        cms.PSet(
            simTracksterCollection=cms.InputTag("hltTiclSimTracksters", "fromBoundarySimCluster"),
            hitToSimClusterMap=cms.InputTag("hltHitToBoundarySimClusterAssociator")
        ),
        cms.PSet(
            simTracksterCollection=cms.InputTag("hltTiclSimTracksters", "fromCaloParticle"),
            hitToSimClusterMap=cms.InputTag("hltHitToCPSimClusterAssociator")
        ),
    )
)

hltHgcalAssociatorsTask = cms.Task(hltRecHitMapProducer,
                                   hltScAssocByEnergyScoreProducer,
                                   SimClusterToCaloParticleAssociation,
                                   hltLayerClusterSimClusterAssociationProducer, hltLayerClusterBoundaryTrackSimClusterAssociationProducer, hltLayerClusterCaloParticleSimClusterAssociationProducer,
                                   hltAllLayerClusterToTracksterAssociations,
                                   hltAllTrackstersToSimTrackstersAssociationsByLCs,
                                   hltAllHitToTracksterAssociations,
                                   hltHitToLegacySimClusterAssociator, hltHitToBoundarySimClusterAssociator, hltHitToCPSimClusterAssociator,
                                   hltAllTrackstersToSimTrackstersAssociationsByHits
                                   )
