import FWCore.ParameterSet.Config as cms


superclusterNtuplizer = cms.EDAnalyzer(
    "SuperclusterNtuplizer",
    superclusters=cms.InputTag("ticlTracksterLinksSuperclusteringDNN"),
    linkedTracksters=cms.InputTag(
        "ticlTracksterLinksSuperclusteringDNN",
        "linkedTracksterIdToInputTracksterId",
    ),
    recoSuperclusters=cms.InputTag("ticlEGammaSuperClusterProducer"),
    simTracksters=cms.InputTag("ticlSimTracksters", "fromCPs"),
    recoToSimAssociator=cms.InputTag(
        "allTrackstersToSimTrackstersAssociationsByLCs",
        "ticlTracksterLinksSuperclusteringDNNToticlSimTrackstersfromCPs",
    ),
    simToRecoAssociator=cms.InputTag(
        "allTrackstersToSimTrackstersAssociationsByLCs",
        "ticlSimTrackstersfromCPsToticlTracksterLinksSuperclusteringDNN",
    ),
    caloParticles=cms.InputTag("mix", "MergedCaloTruth"),
    recoSuperclusterEtThreshold=cms.double(4.0),
)
