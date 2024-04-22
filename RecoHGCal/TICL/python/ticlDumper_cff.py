import FWCore.ParameterSet.Config as cms
from RecoHGCal.TICL.ticlDumper_cfi import ticlDumper as ticlDumper_

ticlDumper = ticlDumper_.clone(
    tracksterCollections = [
        cms.PSet(
            treeName=cms.string("trackstersclue3d"),
            inputTag=cms.InputTag("ticlTrackstersCLUE3DHigh")
        ),
        cms.PSet(
            treeName=cms.string("trackstersmerged"),
            inputTag=cms.InputTag("ticlTrackstersMerge")
        ),

        cms.PSet(
            treeName=cms.string("simtrackstersSC"),
            inputTag=cms.InputTag("ticlSimTracksters"),
            tracksterType=cms.string("SimTracksterSC")
        ),
        cms.PSet(
            treeName=cms.string("simtrackstersCP"),
            inputTag=cms.InputTag("ticlSimTracksters", "fromCPs"),
            tracksterType=cms.string("SimTracksterCP")
        ),
    ],
    
    associators=[
        cms.PSet(
            branchName=cms.string("tsCLUE3D"),
            suffix=cms.string("SC"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationPRbyCLUE3D"),
            tracksterCollection=cms.InputTag("ticlTrackstersCLUE3DHigh"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters")
        ),
        cms.PSet(
            branchName=cms.string("tsCLUE3D"),
            suffix=cms.string("CP"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationLinkingbyCLUE3D"),
            tracksterCollection=cms.InputTag("ticlTrackstersCLUE3DHigh"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters", "fromCPs")
        ),
        
        cms.PSet(
            branchName=cms.string("Mergetstracksters"),
            suffix=cms.string("SC"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationPR"),
            tracksterCollection=cms.InputTag("ticlTrackstersMerge"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters")
        ),
        cms.PSet(
            branchName=cms.string("Mergetracksters"),
            suffix=cms.string("CP"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationLinking"),
            tracksterCollection=cms.InputTag("ticlTrackstersMerge"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters", "fromCPs")
        ),

        cms.PSet(
            branchName=cms.string("Mergetracksters"),
            suffix=cms.string("PU"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationLinkingPU"),
            tracksterCollection=cms.InputTag("ticlTrackstersMerge"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters", "PU")
        ),
    ]
)