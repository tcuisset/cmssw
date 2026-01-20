import FWCore.ParameterSet.Config as cms

from Validation.HGCalValidation.hgcalValidator_cfi import hgcalValidator as _hgcalValidator
from Validation.HGCalValidation.HLT_TICLIterLabels_cff import hltTiclIterLabels as _hltTiclIterLabels

hltAssociatorInstances = []

for labelts in _hltTiclIterLabels:
    for labelsts in ["hltTiclSimTrackstersfromBoundarySimCluster", 'hltTiclSimTrackstersfromCaloParticle']:
        hltAssociatorInstances.append(labelts+'To'+labelsts)
        hltAssociatorInstances.append(labelsts+'To'+labelts)

hltHgcalValidator = _hgcalValidator.clone(
    LayerClustersInputMask = cms.VInputTag("hltTiclTrackstersCLUE3DHigh", "hltTiclSimTracksters:fromCaloParticle", "hltTiclSimTracksters:fromLegacySimCluster", "hltTiclSimTracksters:fromBoundarySimCluster"),
    label_tst = cms.VInputTag(*[cms.InputTag(label) for label in _hltTiclIterLabels] + [cms.InputTag("hltTiclSimTracksters", "fromCaloParticle"), cms.InputTag("hltTiclSimTracksters", "fromLegacySimCluster"), cms.InputTag("hltTiclSimTracksters", "fromBoundarySimCluster")]),
    allTracksterTracksterAssociatorsLabels = cms.VInputTag( *[cms.InputTag('hltAllTrackstersToSimTrackstersAssociationsByLCs:'+associator) for associator in hltAssociatorInstances] ),
    allTracksterTracksterByHitsAssociatorsLabels = cms.VInputTag( *[cms.InputTag('hltAllTrackstersToSimTrackstersAssociationsByHits:'+associator) for associator in hltAssociatorInstances] ),
    associator = cms.untracked.InputTag("hltLayerClusterCaloParticleSimClusterAssociationProducer"),
    associatorSim = cms.untracked.InputTag("hltLayerClusterBoundaryTrackSimClusterAssociationProducer"),
    dirName = cms.string('HLT/HGCAL/HGCalValidator/'),
    hits = cms.InputTag("hltRecHitMapProducer", "RefProdVectorHGCRecHitCollection"),
    hitMap = cms.InputTag("hltRecHitMapProducer","hgcalRecHitMap"),
    simTrackstersMap = cms.InputTag("hltTiclSimTracksters", "fromBoundarySimCluster"),
    label_layerClustersPlots = cms.string("hltHgcalMergeLayerClusters"),
    label_lcl = cms.InputTag("hltMergeLayerClusters"),
    label_simTS = cms.InputTag("hltTiclSimTracksters", "fromBoundarySimCluster"),
    label_simTSFromCP = cms.InputTag("hltTiclSimTracksters","fromCaloParticle"),
    recoTracks = cms.InputTag("hltGeneralTracks"),
    simTiclCandidates = cms.InputTag("hltTiclSimTICLCandidatesFromBoundary"),
    ticlCandidates = cms.string('hltTiclCandidate'),
    ticlTrackstersMerge = cms.InputTag("hltTiclTrackstersMerge"),
    mergeRecoToSimAssociator = cms.InputTag("hltAllTrackstersToSimTrackstersAssociationsByLCs","hltTiclTrackstersMergeTohltTiclSimTrackstersfromCaloParticle"),
    mergeSimToRecoAssociator = cms.InputTag("hltAllTrackstersToSimTrackstersAssociationsByLCs","hltTiclSimTrackstersfromCaloParticleTohltTiclTrackstersMerge"),
)

from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5

lcInputMask_v5  = ["hltTiclTrackstersCLUE3DHigh"]
lcInputMask_v5.extend([cms.InputTag("hltTiclSimTracksters", "fromCaloParticle"), cms.InputTag("hltTiclSimTracksters", "fromBoundarySimCluster")])

ticl_v5.toModify(hltHgcalValidator,
                 LayerClustersInputMask = cms.VInputTag(lcInputMask_v5),
                 ticlTrackstersMerge = cms.InputTag("hltTiclCandidate"),
                 isticlv5 = cms.untracked.bool(True),
                 mergeSimToRecoAssociator = cms.InputTag("hltAllTrackstersToSimTrackstersAssociationsByLCs:hltTiclSimTrackstersfromCaloParticleTohltTiclCandidate"),
                 mergeRecoToSimAssociator = cms.InputTag("hltAllTrackstersToSimTrackstersAssociationsByLCs:hltTiclCandidateTohltTiclSimTrackstersfromCaloParticle"),
                 )

