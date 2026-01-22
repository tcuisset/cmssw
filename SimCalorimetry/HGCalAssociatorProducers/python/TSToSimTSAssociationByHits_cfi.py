import FWCore.ParameterSet.Config as cms
from SimCalorimetry.HGCalAssociatorProducers.HitToTracksterAssociation_cfi import *
from SimCalorimetry.HGCalAssociatorProducers.AllTracksterToSimTracksterAssociatorsByHitsProducer_cfi import AllTracksterToSimTracksterAssociatorsByHitsProducer
from RecoHGCal.TICL.iterativeTICL_cff import ticlIterLabels

allTrackstersToSimTrackstersAssociationsByHits = AllTracksterToSimTracksterAssociatorsByHitsProducer.clone(    
    tracksterCollections = cms.VInputTag(
        *[cms.InputTag(label) for label in ticlIterLabels]
    ),
    simTracksters = cms.VPSet(
        cms.PSet(
            simTracksterCollection=cms.InputTag("ticlSimTracksters", "fromLegacySimCluster"),
            hitToSimClusterMap=cms.InputTag("hitToLegacySimClusterAssociator")
        ),
        cms.PSet(
            simTracksterCollection=cms.InputTag("ticlSimTracksters", "fromBoundarySimCluster"),
            hitToSimClusterMap=cms.InputTag("hitToBoundarySimClusterAssociator")
        ),
        cms.PSet(
            simTracksterCollection=cms.InputTag("ticlSimTracksters", "fromMergedSimCluster"),
            hitToSimClusterMap=cms.InputTag("hitToMergedSimClusterAssociator")
        ),
        cms.PSet(
            simTracksterCollection=cms.InputTag("ticlSimTracksters", "fromCaloParticle"),
            hitToSimClusterMap=cms.InputTag("hitToCPSimClusterAssociator")
        ),
    )
)

