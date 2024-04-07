import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.ticlLayerTileProducer_cfi import ticlLayerTileProducer

from RecoHGCal.TICL.CLUE3DEM_cff import *
from RecoHGCal.TICL.CLUE3DHAD_cff import *
from RecoHGCal.TICL.pfTICLProducer_cfi import pfTICLProducer as _pfTICLProducer

from RecoHGCal.TICL.ticlLayerTileProducer_cfi import ticlLayerTileProducer
from RecoHGCal.TICL.pfTICLProducer_cfi import pfTICLProducer as _pfTICLProducer
from RecoHGCal.TICL.tracksterSelectionTf_cfi import *

from RecoHGCal.TICL.tracksterLinksProducer_cfi import tracksterLinksProducer as _tracksterLinksProducer
from RecoHGCal.TICL.ticlEGammaSuperClusterProducer_cfi import ticlEGammaSuperClusterProducer
from RecoHGCal.TICL.ticlCandidateProducer_cfi import ticlCandidateProducer as _ticlCandidateProducer
from RecoHGCal.Configuration.RecoHGCal_EventContent_cff import customiseForTICLv5EventContent
from RecoHGCal.TICL.iterativeTICL_cff import ticlIterLabels, ticlIterLabelsMerge
from RecoHGCal.TICL.ticlDumper_cff import ticlDumper
from RecoHGCal.TICL.mergedTrackstersProducer_cfi import mergedTrackstersProducer as _mergedTrackstersProducer
from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import tracksterSimTracksterAssociationLinkingbyCLUE3D as _tracksterSimTracksterAssociationLinkingbyCLUE3D
from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import tracksterSimTracksterAssociationPRbyCLUE3D  as _tracksterSimTracksterAssociationPRbyCLUE3D 
from Validation.HGCalValidation.HGCalValidator_cff import hgcalValidatorv5 
from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import HGCalUncalibRecHit
from RecoHGCal.TICL.SimTracksters_cff import ticlSimTracksters

from RecoHGCal.TICL.FastJetStep_cff import ticlTrackstersFastJet
from RecoHGCal.TICL.EMStep_cff import ticlTrackstersEM, ticlTrackstersHFNoseEM
from RecoHGCal.TICL.TrkStep_cff import ticlTrackstersTrk, ticlTrackstersHFNoseTrk
from RecoHGCal.TICL.MIPStep_cff import ticlTrackstersMIP, ticlTrackstersHFNoseMIP
from RecoHGCal.TICL.HADStep_cff import ticlTrackstersHAD, ticlTrackstersHFNoseHAD
from RecoHGCal.TICL.CLUE3DEM_cff import ticlTrackstersCLUE3DEM
from RecoHGCal.TICL.CLUE3DHAD_cff import ticlTrackstersCLUE3DHAD
from RecoHGCal.TICL.CLUE3DHighStep_cff import ticlTrackstersCLUE3DHigh
from RecoHGCal.TICL.TrkEMStep_cff import ticlTrackstersTrkEM, filteredLayerClustersHFNoseTrkEM


def customiseForTICLv5(process, enableDumper = False, enableSuperclusteringDNN=True):
    """ enableSuperclusteringDNN : if True, use the superlcustering DNN for electron seeds. If False, use Mustache (fed from CLUE3D EM tracksters) """

    process.HGCalUncalibRecHit.computeLocalTime = cms.bool(True)
    process.ticlSimTracksters.computeLocalTime = cms.bool(True)

    process.ticlTrackstersFastJet.pluginPatternRecognitionByFastJet.computeLocalTime = cms.bool(True)

    process.ticlTrackstersEM.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)
    process.ticlTrackstersHFNoseEM.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)

    process.ticlTrackstersTrk.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)
    process.ticlTrackstersHFNoseTrk.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)

    process.ticlTrackstersMIP.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)
    process.ticlTrackstersHFNoseMIP.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)

    process.ticlTrackstersHAD.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)
    process.ticlTrackstersHFNoseHAD.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)

    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.computeLocalTime = cms.bool(True)
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.computeLocalTime = cms.bool(True)
    process.ticlTrackstersCLUE3DHigh.pluginPatternRecognitionByCLUE3D.computeLocalTime = cms.bool(True)

    process.ticlTrackstersTrkEM.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)
    process.ticlTrackstersHFNoseTrkEM.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)

    process.ticlLayerTileTask = cms.Task(ticlLayerTileProducer)

    process.ticlIterationsTask = cms.Task(
        ticlCLUE3DEMStepTask,
        ticlCLUE3DHADStepTask,
    )

    process.ticlTracksterLinks = _tracksterLinksProducer.clone()
    process.ticlTracksterLinksTask = cms.Task(process.ticlTracksterLinks, process.ticlTracksterLinksSuperclustering)

    process.ticlCandidate = _ticlCandidateProducer.clone()
    process.ticlCandidateTask = cms.Task(process.ticlCandidate)

    process.tracksterSimTracksterAssociationLinkingbyCLUE3DEM = _tracksterSimTracksterAssociationLinkingbyCLUE3D.clone(
        label_tst = cms.InputTag("ticlTrackstersCLUE3DEM")
        )
    process.tracksterSimTracksterAssociationPRbyCLUE3DEM = _tracksterSimTracksterAssociationPRbyCLUE3D.clone(
        label_tst = cms.InputTag("ticlTrackstersCLUE3DEM")
        )
    process.tracksterSimTracksterAssociationLinkingbyCLUE3DHAD = _tracksterSimTracksterAssociationLinkingbyCLUE3D.clone(
        label_tst = cms.InputTag("ticlTrackstersCLUE3DHAD")
        )
    process.tracksterSimTracksterAssociationPRbyCLUE3DHAD = _tracksterSimTracksterAssociationPRbyCLUE3D.clone(
        label_tst = cms.InputTag("ticlTrackstersCLUE3DHAD")
        )
    # Linking CaloParticleSimTrackster to EM superclusters
    process.tracksterSimTracksterAssociationLinkingSuperclustering = _tracksterSimTracksterAssociationLinkingbyCLUE3D.clone(
        label_tst = cms.InputTag("ticlTracksterLinksSuperclustering"),
        )
    process.tracksterSimTracksterAssociationPRSuperclustering = _tracksterSimTracksterAssociationPRbyCLUE3D.clone(
        label_tst = cms.InputTag("ticlTracksterLinksSuperclustering"),
        )

    process.mergedTrackstersProducer = _mergedTrackstersProducer.clone()    

    process.tracksterSimTracksterAssociationLinkingbyCLUE3D = _tracksterSimTracksterAssociationLinkingbyCLUE3D.clone(
        label_tst = cms.InputTag("mergedTrackstersProducer")
        )
    process.tracksterSimTracksterAssociationPRbyCLUE3D = _tracksterSimTracksterAssociationPRbyCLUE3D.clone(
        label_tst = cms.InputTag("mergedTrackstersProducer")
        )
    process.iterTICLTask = cms.Task(process.ticlLayerTileTask,
                                     process.ticlIterationsTask,
                                     process.ticlTracksterLinksTask,
                                     process.ticlCandidateTask)
    process.particleFlowClusterHGCal.initialClusteringStep.tracksterSrc = "ticlCandidate"
    process.globalrecoTask.remove(process.ticlTrackstersMerge)

    # modifying superclustering
    if enableSuperclusteringDNN:
        process.ticlEGammaSuperClusterProducer = ticlEGammaSuperClusterProducer
        process.FEVTDEBUGHLToutput.outputCommands.extend(["keep *_ticlEGammaSuperClusterProducer_*_*"])
        process.particleFlowSuperClusteringTask.replace(process.particleFlowSuperClusterHGCal, process.ticlEGammaSuperClusterProducer)
        try:
            process.mergedSuperClustersHGC.src[1] = "ticlEGammaSuperClusterProducer" # original config : cms.EDProducer("SuperClusterMerger", src = cms.VInputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel", "particleFlowSuperClusterHGCal"))
        except AttributeError: pass # in case we run without DQM
        process.ecalDrivenElectronSeeds.endcapSuperClusters = cms.InputTag("ticlEGammaSuperClusterProducer")
        process.photonCoreHGC.scIslandEndcapProducer = cms.InputTag("ticlEGammaSuperClusterProducer")
    else:
        # make CLUE3D EM tracksters flow into Mustache (instead of ticlCandidate)

        # particleFlowClusterHGCal is used for particleFlowSuperClusterHGCal as well as for hgcalHitCalibration (which we don't want to change)
        # so duplicate particleFlowClusterHGCal to make it consume CLUE3D EM tracksters, keeping the old config in in parallel for hgcalHitCalibration
        process.particleFlowClusterHGCalCLUE3DEM = process.particleFlowClusterHGCal.clone()
        process.particleFlowClusterHGCalCLUE3DEM.initialClusteringStep.tracksterSrc = cms.InputTag("ticlTrackstersCLUE3DEM")
        process.particleFlowClusterHGCalCLUE3DEM.filterByTracksterPID = cms.bool(False)
        process.hgcalLocalRecoTask.add(process.particleFlowClusterHGCalCLUE3DEM)
        process.particleFlowSuperClusterHGCal.PFClusters = cms.InputTag("particleFlowClusterHGCalCLUE3DEM")

    process.tracksterSimTracksterAssociationLinking.label_tst = cms.InputTag("ticlCandidate")
    process.tracksterSimTracksterAssociationPR.label_tst = cms.InputTag("ticlCandidate")

    process.tracksterSimTracksterAssociationLinkingPU.label_tst = cms.InputTag("ticlTracksterLinks")
    process.tracksterSimTracksterAssociationPRPU.label_tst = cms.InputTag("ticlTracksterLinks")
    process.mergeTICLTask = cms.Task()
    process.pfTICL.ticlCandidateSrc = cms.InputTag("ticlCandidate") 
    process.hgcalAssociators = cms.Task(process.mergedTrackstersProducer, process.lcAssocByEnergyScoreProducer, process.layerClusterCaloParticleAssociationProducer,
                            process.scAssocByEnergyScoreProducer, process.layerClusterSimClusterAssociationProducer,
                            process.lcSimTSAssocByEnergyScoreProducer, process.layerClusterSimTracksterAssociationProducer,
                            process.simTsAssocByEnergyScoreProducer,  process.simTracksterHitLCAssociatorByEnergyScoreProducer,
                            process.tracksterSimTracksterAssociationLinking, process.tracksterSimTracksterAssociationPR,
                            process.tracksterSimTracksterAssociationLinkingbyCLUE3D, process.tracksterSimTracksterAssociationPRbyCLUE3D,
                            process.tracksterSimTracksterAssociationLinkingbyCLUE3DEM, process.tracksterSimTracksterAssociationPRbyCLUE3DEM,
                            process.tracksterSimTracksterAssociationLinkingbyCLUE3DHAD, process.tracksterSimTracksterAssociationPRbyCLUE3DHAD,
                            process.tracksterSimTracksterAssociationLinkingSuperclustering, process.tracksterSimTracksterAssociationPRSuperclustering,
                            process.tracksterSimTracksterAssociationLinkingPU, process.tracksterSimTracksterAssociationPRPU
                            )

    process.hgcalValidatorv5 = hgcalValidatorv5.clone(
        ticlTrackstersMerge = cms.InputTag("ticlCandidate"),
        trackstersclue3d = cms.InputTag("mergedTrackstersProducer")
    )
    process.hgcalValidatorSequence = cms.Sequence(process.hgcalValidatorv5)
    process.hgcalValidation = cms.Sequence(process.hgcalSimHitValidationEE+process.hgcalSimHitValidationHEF+process.hgcalSimHitValidationHEB+process.hgcalDigiValidationEE+process.hgcalDigiValidationHEF+process.hgcalDigiValidationHEB+process.hgcalRecHitValidationEE+process.hgcalRecHitValidationHEF+process.hgcalRecHitValidationHEB+process.hgcalHitValidationSequence+process.hgcalValidatorSequence+process.hgcalTiclPFValidation+process.hgcalPFJetValidation)
    process.globalValidationHGCal = cms.Sequence(process.hgcalValidation)
    process.validation_step9 = cms.EndPath(process.globalValidationHGCal)
    if(enableDumper):
        process.ticlDumper = ticlDumper.clone(
            tracksterCollections=[
                cms.PSet(
                    treeName=cms.string("trackstersCLUE3DEM"),
                    inputTag=cms.InputTag("ticlTrackstersCLUE3DEM")
                ),
                cms.PSet(
                    treeName=cms.string("trackstersSuperclustering"),
                    inputTag=cms.InputTag("ticlTracksterLinksSuperclustering")
                ),
                cms.PSet(
                    treeName=cms.string("trackstersCLUE3DHAD"),
                    inputTag=cms.InputTag("ticlTrackstersCLUE3DHAD")
                ),
                cms.PSet(
                    treeName=cms.string("trackstersMerged"),
                    inputTag=cms.InputTag("mergedTrackstersProducer")
                ),
                cms.PSet(
                    treeName=cms.string("trackstersTiclCandidate"),
                    inputTag=cms.InputTag("ticlCandidate")
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
            recoSuperClusters = cms.InputTag("ticlEGammaSuperClusterProducer" if enableSuperclusteringDNN else "particleFlowSuperClusterHGCal"),
            recoSuperClusters_sourceTracksterCollection = cms.InputTag("ticlTrackstersCLUE3DEM"),
            ticlcandidates = cms.InputTag("ticlCandidate"),
                associators=[
                    cms.PSet(
                        branchName=cms.string("tsCLUE3DEM"),
                        suffix=cms.string("SC"),
                        associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationPRbyCLUE3DEM"),
                        tracksterCollection=cms.InputTag("ticlTrackstersCLUE3DEM"),
                        simTracksterCollection=cms.InputTag("ticlSimTracksters")
                    ),
                    cms.PSet(
                        branchName=cms.string("tsCLUE3DEM"),
                        suffix=cms.string("CP"),
                        associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationLinkingbyCLUE3DEM"),
                        tracksterCollection=cms.InputTag("ticlTrackstersCLUE3DEM"),
                        simTracksterCollection=cms.InputTag("ticlSimTracksters", "fromCPs")
                    ),

                    cms.PSet(
                        branchName=cms.string("tsCLUE3DHAD"),
                        suffix=cms.string("SC"),
                        associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationPRbyCLUE3DHAD"),
                        tracksterCollection=cms.InputTag("ticlTrackstersCLUE3DHAD"),
                        simTracksterCollection=cms.InputTag("ticlSimTracksters")
                    ),
                    cms.PSet(
                        branchName=cms.string("tsCLUE3DHAD"),
                        suffix=cms.string("CP"),
                        associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationLinkingbyCLUE3DHAD"),
                        tracksterCollection=cms.InputTag("ticlTrackstersCLUE3DHAD"),
                        simTracksterCollection=cms.InputTag("ticlSimTracksters", "fromCPs")
                    ),

                    cms.PSet(
                        branchName=cms.string("tsSuperclusters"),
                        suffix=cms.string("SC"),
                        associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationPRSuperclustering"),
                        tracksterCollection=cms.InputTag("ticlTracksterLinksSuperclustering"),
                        simTracksterCollection=cms.InputTag("ticlSimTracksters")
                    ),
                    cms.PSet(
                        branchName=cms.string("tsSuperclusters"),
                        suffix=cms.string("CP"),
                        associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationLinkingSuperclustering"),
                        tracksterCollection=cms.InputTag("ticlTracksterLinksSuperclustering"),
                        simTracksterCollection=cms.InputTag("ticlSimTracksters", "fromCPs")
                    ),

                    cms.PSet(
                        branchName=cms.string("Mergetracksters"),
                        suffix=cms.string("SC"),
                        associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationPRbyCLUE3D"),
                        tracksterCollection=cms.InputTag("mergedTrackstersProducer"),
                        simTracksterCollection=cms.InputTag("ticlSimTracksters")
                    ),
                    cms.PSet(
                        branchName=cms.string("Mergetracksters"),
                        suffix=cms.string("CP"),
                        associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationLinkingbyCLUE3D"),
                        tracksterCollection=cms.InputTag("mergedTrackstersProducer"),
                        simTracksterCollection=cms.InputTag("ticlSimTracksters", "fromCPs")
                    ),

                    cms.PSet(
                        branchName=cms.string("Mergetracksters"),
                        suffix=cms.string("PU"),
                        associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationLinkingPU"),
                        tracksterCollection=cms.InputTag("ticlTracksterLinks"),
                        simTracksterCollection=cms.InputTag("ticlSimTracksters", "PU")
                    ),
                ]
        )
        process.TFileService = cms.Service("TFileService",
                                           fileName=cms.string("histo.root")
                                           )
        process.FEVTDEBUGHLToutput_step = cms.EndPath(
            process.FEVTDEBUGHLToutput + process.ticlDumper)


    process = customiseForTICLv5EventContent(process)

    return process
