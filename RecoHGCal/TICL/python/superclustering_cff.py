import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.tracksterLinksProducer_cfi import tracksterLinksProducer as _tracksterLinksProducer

ticlTracksterLinksSuperclustering = _tracksterLinksProducer.clone(
    linkingPSet = cms.PSet(
        type=cms.string("SuperClustering"),
        algo_verbosity=cms.int32(0)
    ),
    tracksters_collections = [cms.InputTag("ticlTrackstersCLUE3DEM")],
)
