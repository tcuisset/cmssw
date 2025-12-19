import FWCore.ParameterSet.Config as cms
from RecoHGCal.TICL.ticlDumper_cfi import ticlDumper as ticlDumper_

from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5
from Configuration.ProcessModifiers.ticl_superclustering_dnn_cff import ticl_superclustering_dnn
from Configuration.ProcessModifiers.ticl_superclustering_mustache_pf_cff import ticl_superclustering_mustache_pf
from Configuration.ProcessModifiers.ticl_superclustering_mustache_ticl_cff import ticl_superclustering_mustache_ticl


from RecoHGCal.TICL.iterativeTICL_cff import ticlIterLabels, associatorsInstances


dumperAssociators = []

for tracksterIteration in ticlIterLabels:
    simTrackstersCollection = "ticlSimTrackstersfromCaloParticle"
    dumperAssociators.append(
        cms.PSet(
            branchName=cms.string(tracksterIteration),
            suffix=cms.string("CP"),
            associatorRecoToSimInputTag=cms.InputTag(f"allTrackstersToSimTrackstersAssociationsByLCs:{tracksterIteration}To{simTrackstersCollection}"),
            associatorSimToRecoInputTag=cms.InputTag(f"allTrackstersToSimTrackstersAssociationsByLCs:{simTrackstersCollection}To{tracksterIteration}")
        )
    )

    simTrackstersCollection = "ticlSimTrackstersfromBoundarySimCluster"
    dumperAssociators.append(
        cms.PSet(
            branchName=cms.string(tracksterIteration),
            suffix=cms.string("SC"),
            associatorRecoToSimInputTag=cms.InputTag(f"allTrackstersToSimTrackstersAssociationsByLCs:{tracksterIteration}To{simTrackstersCollection}"),
            associatorSimToRecoInputTag=cms.InputTag(f"allTrackstersToSimTrackstersAssociationsByLCs:{simTrackstersCollection}To{tracksterIteration}")
        )
    )


ticlDumper = ticlDumper_.clone(
    tracksterCollections = [*[cms.PSet(treeName=cms.string(label), inputTag=cms.InputTag(label)) for label in ticlIterLabels],
        cms.PSet(
            treeName=cms.string("simtrackstersSC"),
            inputTag=cms.InputTag("ticlSimTracksters", "fromBoundarySimCluster"),
            tracksterType=cms.string("SimTracksterSC")
        ),
        cms.PSet(
            treeName=cms.string("simtrackstersCP"),
            inputTag=cms.InputTag("ticlSimTracksters", "fromCaloParticle"),
            tracksterType=cms.string("SimTracksterCP")
        ),
    ],

    associators=dumperAssociators.copy(),
    saveSuperclustering = cms.bool(False)
)

ticl_v5.toModify(ticlDumper, 
                 ticlcandidates = cms.InputTag("ticlCandidate"), 
                 recoSuperClusters_sourceTracksterCollection=cms.InputTag("ticlTrackstersCLUE3DHigh"), 
                 saveSuperclustering = cms.bool(True), 
                 trackstersInCand=cms.InputTag("ticlCandidate"))

(ticl_v5 & ticl_superclustering_mustache_pf).toModify(ticlDumper, saveSuperclustering=False, recoSuperClusters_sourceTracksterCollection=cms.InputTag("ticlTrackstersCLUE3DHigh"))
