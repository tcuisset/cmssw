import FWCore.ParameterSet.Config as cms
from SimCalorimetry.HGCalAssociatorProducers.hitToTracksterAssociator_cfi import hitToTracksterAssociator

# the "single" hitToTrackstersAssociation are not used (only the allHitToTracksterAssociations one is used)
hitToTrackstersAssociationLinking = hitToTracksterAssociator.clone(
    tracksters = cms.InputTag("ticlTrackstersMerge"),
)


hitToTrackstersAssociationPR = hitToTracksterAssociator.clone(
    tracksters = cms.InputTag("ticlTrackstersCLUE3DHigh"),
)

hitToSimTracksterAssociation = hitToTracksterAssociator.clone(
    tracksters = cms.InputTag("ticlSimTracksters", "fromLegacySimCluster"),
)

hitToSimTracksterFromCPsAssociation = hitToTracksterAssociator.clone(
    tracksters = cms.InputTag("ticlSimTracksters", "fromCaloParticle"),
)


from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5

ticl_v5.toModify(hitToTrackstersAssociationLinking, tracksters = cms.InputTag("ticlCandidate"))

from SimCalorimetry.HGCalAssociatorProducers.AllHitToTracksterAssociatorsProducer_cfi import AllHitToTracksterAssociatorsProducer
from RecoHGCal.TICL.iterativeTICL_cff import ticlIterLabels

allHitToTracksterAssociations = AllHitToTracksterAssociatorsProducer.clone(    
    tracksterCollections = cms.VInputTag(
        *[cms.InputTag(label) for label in ticlIterLabels],
        cms.InputTag("ticlSimTracksters", "fromLegacySimCluster"),
        cms.InputTag("ticlSimTracksters", "fromBoundarySimCluster"),
        cms.InputTag("ticlSimTracksters", "fromMergedSimCluster"),
        cms.InputTag("ticlSimTracksters", "fromCaloParticle"),
    )
)


