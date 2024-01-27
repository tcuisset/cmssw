from RecoHGCal.TICL.TICLDumper_cfi import ticlDumper as ticlDumper_

ticlDumper = ticlDumper_.clone(
    tracksterCollections = [
        cms.PSet(
            treeName="trackstersclue3d",
            inputTag=cms.InputTag("ticlTrackstersCLUE3DHigh")
        ),
        cms.PSet(
            treeName="trackstersmerged",
            inputTag=cms.InputTag("ticlTrackstersMerge")
        ),
    ],
    
    associators=[
        cms.PSet(
            branchName="tsCLUE3D",
            suffix="SC",
            associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationPRbyCLUE3D"),
            tracksterCollection=cms.InputTag("ticlTrackstersCLUE3DHigh"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters")
        ),
        cms.PSet(
            branchName="tsCLUE3D",
            suffix="CP",
            associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationLinkingbyCLUE3D"),
            tracksterCollection=cms.InputTag("ticlTrackstersCLUE3DHigh"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters", "fromCPs")
        ),
        
        cms.PSet(
            branchName="Mergetstracksters",
            suffix="SC",
            associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationPR"),
            tracksterCollection=cms.InputTag("ticlTrackstersMerge"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters")
        ),
        cms.PSet(
            branchName="Mergetracksters",
            suffix="CP",
            associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationLinking"),
            tracksterCollection=cms.InputTag("ticlTrackstersMerge"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters", "fromCPs")
        ),

        cms.PSet(
            branchName="Mergetracksters",
            suffix="PU",
            associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationLinkingPU"),
            tracksterCollection=cms.InputTag("ticlTrackstersMerge"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters", "PU")
        ),
    ]
)