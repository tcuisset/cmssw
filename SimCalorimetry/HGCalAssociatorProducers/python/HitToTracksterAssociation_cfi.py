import FWCore.ParameterSet.Config as cms
from SimCalorimetry.HGCalAssociatorProducers.hitToTracksterAssociator_cfi import hitToTracksterAssociator as _hitToTracksterAssociator

# the "single" hitToTrackstersAssociation are not used (only the allHitToTracksterAssociations one is used)
hitToTrackstersAssociationLinking = _hitToTracksterAssociator.clone(
    tracksters = cms.InputTag("ticlTrackstersMerge"),
)


hitToTrackstersAssociationPR = _hitToTracksterAssociator.clone(
    tracksters = cms.InputTag("ticlTrackstersCLUE3DHigh"),
)

hitToSimTracksterAssociation = _hitToTracksterAssociator.clone(
    tracksters = cms.InputTag("ticlSimTracksters", "fromLegacySimCluster"),
)

hitToSimTracksterFromCPsAssociation = _hitToTracksterAssociator.clone(
    tracksters = cms.InputTag("ticlSimTracksters", "fromCaloParticle"),
)


from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5

ticl_v5.toModify(hitToTrackstersAssociationLinking, tracksters = cms.InputTag("ticlCandidate"))

from SimCalorimetry.HGCalAssociatorProducers.AllHitToTracksterAssociatorsProducer_cfi import AllHitToTracksterAssociatorsProducer as _AllHitToTracksterAssociatorsProducer
from RecoHGCal.TICL.iterativeTICL_cff import ticlIterLabels

allHitToTracksterAssociations = _AllHitToTracksterAssociatorsProducer.clone(    
    tracksterCollections = cms.VInputTag(
        *[cms.InputTag(label) for label in ticlIterLabels],
        cms.InputTag("ticlSimTracksters", "fromLegacySimCluster"),
        cms.InputTag("ticlSimTracksters", "fromBoundarySimCluster"),
        cms.InputTag("ticlSimTracksters", "fromMergedSimCluster"),
        cms.InputTag("ticlSimTracksters", "fromCaloParticle"),
    )
)


