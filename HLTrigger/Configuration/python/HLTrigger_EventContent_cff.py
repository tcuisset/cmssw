import FWCore.ParameterSet.Config as cms

# EventContent for HLT related products.

# This file exports the following EventContent blocks:
#   HLTriggerRAW  HLTriggerRECO  HLTriggerAOD (without DEBUG products)
#   HLTDebugRAW   HLTDebugFEVT                (with    DEBUG products)
#   HLTScouting                               (with Scouting products)
#
# as these are used in Configuration/EventContent
#
HLTriggerRAW  = cms.PSet(
    outputCommands = cms.vstring( *(
        'drop *_hlt*_*_*',
        'keep FEDRawDataCollection_rawDataCollector_*_*',
        'keep FEDRawDataCollection_source_*_*',
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*',
        'keep edmTriggerResults_*_*_*',
        'keep triggerTriggerEvent_*_*_*',
        'keep *_hltFEDSelectorL1_*_*',
        'keep *_hltScoutingEgammaPacker_*_*',
        'keep *_hltScoutingMuonPacker_*_*',
        'keep *_hltScoutingPFPacker_*_*',
        'keep *_hltScoutingPrimaryVertexPacker_*_*',
        'keep *_hltScoutingTrackPacker_*_*',
        'keep edmTriggerResults_*_*_*'
    ) )
)

HLTriggerRECO = cms.PSet(
    outputCommands = cms.vstring( *(
        'drop *_hlt*_*_*',
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*',
        'keep edmTriggerResults_*_*_*',
        'keep triggerTriggerEvent_*_*_*',
        'keep *_hltFEDSelectorL1_*_*',
        'keep *_hltScoutingEgammaPacker_*_*',
        'keep *_hltScoutingMuonPacker_*_*',
        'keep *_hltScoutingPFPacker_*_*',
        'keep *_hltScoutingPrimaryVertexPacker_*_*',
        'keep *_hltScoutingTrackPacker_*_*',
        'keep edmTriggerResults_*_*_*'
    ) )
)

HLTriggerAOD  = cms.PSet(
    outputCommands = cms.vstring( *(
        'drop *_hlt*_*_*',
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*',
        'keep edmTriggerResults_*_*_*',
        'keep triggerTriggerEvent_*_*_*',
        'keep *_hltFEDSelectorL1_*_*',
        'keep *_hltScoutingEgammaPacker_*_*',
        'keep *_hltScoutingMuonPacker_*_*',
        'keep *_hltScoutingPFPacker_*_*',
        'keep *_hltScoutingPrimaryVertexPacker_*_*',
        'keep *_hltScoutingTrackPacker_*_*',
        'keep edmTriggerResults_*_*_*'
    ) )
)

HLTDebugRAW   = cms.PSet(
    outputCommands = cms.vstring( *(
        'drop *_hlt*_*_*',
        'keep *_hltAK4CaloJetsCorrectedIDPassed_*_*',
        'keep *_hltAK4CaloJetsIDPassed_*_*',
        'keep *_hltAK4CaloJets_*_*',
        'keep *_hltAK4PFJetsCorrected_*_*',
        'keep *_hltAK4PFJetsForTaus_*_*',
        'keep *_hltAK4PFJets_*_*',
        'keep *_hltAlCaEtaEBRechitsToDigis_*_*',
        'keep *_hltAlCaEtaEERechitsToDigis_*_*',
        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegional_etaEcalRecHitsES_*',
        'keep *_hltAlCaPi0EBRechitsToDigis_*_*',
        'keep *_hltAlCaPi0EERechitsToDigis_*_*',
        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegional_pi0EcalRecHitsES_*',
        'keep *_hltAlcaPixelClusterCounts_*_*',
        'keep *_hltBSoftMuonMu5L3_*_*',
        'keep *_hltCsc2DRecHits_*_*',
        'keep *_hltCscSegments_*_*',
        'keep *_hltDeepBLifetimeTagInfosPF_*_*',
        'keep *_hltDeepCombinedSecondaryVertexBJetTagsCalo_*_*',
        'keep *_hltDeepCombinedSecondaryVertexBJetTagsInfosCalo_*_*',
        'keep *_hltDeepCombinedSecondaryVertexBJetTagsInfos_*_*',
        'keep *_hltDeepCombinedSecondaryVertexBJetTagsPF_*_*',
        'keep *_hltDeepSecondaryVertexTagInfosPF_*_*',
        'keep *_hltDisplacedhltIter4PFlowTrackSelectionHighPurity_*_*',
        'keep *_hltDoubletRecoveryPFlowTrackSelectionHighPurity_*_*',
        'keep *_hltDt4DSegments_*_*',
        'keep *_hltEcalPhiSymFilter_*_*',
        'keep *_hltEcalRecHit_*_*',
        'keep *_hltEgammaCandidates_*_*',
        'keep *_hltEgammaGsfElectrons_*_*',
        'keep *_hltEgammaGsfTracks_*_*',
        'keep *_hltFEDSelectorGEM_*_*',
        'keep *_hltFEDSelectorTCDS_*_*',
        'keep *_hltFastPVPixelTracksMerger_*_*',
        'keep *_hltFastPVPixelTracksRecover_*_*',
        'keep *_hltFastPVPixelTracks_*_*',
        'keep *_hltFastPVPixelVertices_*_*',
        'keep *_hltFastPixelBLifetimeL3Associator_*_*',
        'keep *_hltFastPrimaryVertex_*_*',
        'keep *_hltGtStage2Digis_*_*',
        'keep *_hltGtStage2ObjectMap_*_*',
        'keep *_hltHbhereco_*_*',
        'keep *_hltHfreco_*_*',
        'keep *_hltHoreco_*_*',
        'keep *_hltImpactParameterTagInfos_*_*',
        'keep *_hltInclusiveSecondaryVertexFinderTagInfos_*_*',
        'keep *_hltIsolPixelTrackProdHB_*_*',
        'keep *_hltIsolPixelTrackProdHE_*_*',
        'keep *_hltIter0HighPtTkMuTrackSelectionHighPurity_*_*',
        'keep *_hltIter2MergedForDisplaced_*_*',
        'keep *_hltIterL3GlbMuon_*_*',
        'keep *_hltIterL3MuonAndMuonFromL1Merged_*_*',
        'keep *_hltIterL3MuonMerged_*_*',
        'keep *_hltIterL3MuonsNoID_*_*',
        'keep *_hltIterL3Muons_*_*',
        'keep *_hltIterL3OIMuonTrackSelectionHighPurity_*_*',
        'keep *_hltL2MuonCandidatesNoVtx_*_*',
        'keep *_hltL2MuonCandidates_*_*',
        'keep *_hltL2MuonSeeds_*_*',
        'keep *_hltL2Muons_*_*',
        'keep *_hltL2TauJets_*_*',
        'keep *_hltL3MuonsIOHit_*_*',
        'keep *_hltL3MuonsLinksCombination_*_*',
        'keep *_hltL3MuonsOIHit_*_*',
        'keep *_hltL3MuonsOIState_*_*',
        'keep *_hltL3Muons_*_*',
        'keep *_hltL3NoFiltersNoVtxMuonCandidates_*_*',
        'keep *_hltL3NoFiltersNoVtxMuons_*_*',
        'keep *_hltL3TkFromL2OICombination_*_*',
        'keep *_hltL3TkTracksFromL2IOHit_*_*',
        'keep *_hltL3TkTracksFromL2OIHit_*_*',
        'keep *_hltL3TkTracksFromL2OIState_*_*',
        'keep *_hltL3TkTracksFromL2_*_*',
        'keep *_hltL3TrackCandidateFromL2IOHit_*_*',
        'keep *_hltL3TrackCandidateFromL2OIHit_*_*',
        'keep *_hltL3TrackCandidateFromL2OIState_*_*',
        'keep *_hltL3TrajSeedIOHit_*_*',
        'keep *_hltL3TrajSeedOIHit_*_*',
        'keep *_hltL3TrajSeedOIState_*_*',
        'keep *_hltL3TrajectorySeed_*_*',
        'keep *_hltMergedTracksForBTag_*_*',
        'keep *_hltMergedTracks_*_*',
        'keep *_hltMet_*_*',
        'keep *_hltMuonCSCDigis_*_*',
        'keep *_hltMuonCSCDigis_MuonCSCStripDigi_*',
        'keep *_hltMuonCSCDigis_MuonCSCWireDigi_*',
        'keep *_hltMuonDTDigis_*_*',
        'keep *_hltMuonRPCDigis_*_*',
        'keep *_hltOnlineBeamSpot_*_*',
        'keep *_hltPFJetForBtag_*_*',
        'keep *_hltPFJetForPNetAK8_*_*',
        'keep *_hltPFMETNoMuProducer_*_*',
        'keep *_hltPFMETProducer_*_*',
        'keep *_hltPFMETTypeOne_*_*',
        'keep *_hltPFMuonMerging_*_*',
        'keep *_hltPFTau35Track_*_*',
        'keep *_hltPFTau35_*_*',
        'keep *_hltPPSCalibrationRaw_*_*',
        'keep *_hltParticleFlowForTaus_*_*',
        'keep *_hltParticleFlow_*_*',
        'keep *_hltParticleNetDiscriminatorsJetTagsAK8_*_*',
        'keep *_hltParticleNetDiscriminatorsJetTags_*_*',
        'keep *_hltPixelTracks_*_*',
        'keep *_hltPixelVertices_*_*',
        'keep *_hltRpcRecHits_*_*',
        'keep *_hltSelector4CentralJetsL1FastJet_*_*',
        'keep *_hltSelector8CentralJetsL1FastJet_*_*',
        'keep *_hltSelectorJets20L1FastJet_*_*',
        'keep *_hltSiPixelClustersCache_*_*',
        'keep *_hltSiPixelClusters_*_*',
        'keep *_hltSiStripClusterizerForRawPrime_*_*',
        'keep *_hltSiStripClusters2ApproxClusters_*_*',
        'keep *_hltSiStripRawToClustersFacility_*_*',
        'keep *_hltTowerMakerForAll_*_*',
        'keep *_hltTriggerSummaryAOD_*_*',
        'keep *_hltTriggerSummaryRAW_*_*',
        'keep *_hltTrimmedPixelVertices_*_*',
        'keep *_hltVerticesL3_*_*',
        'keep *_hltVerticesPFFilter_*_*',
        'keep *_hltVerticesPFSelector_*_*',
        'keep DetIds_hltSiStripRawToDigi_*_*',
        'keep FEDRawDataCollection_rawDataCollector_*_*',
        'keep FEDRawDataCollection_rawDataRepacker_*_*',
        'keep FEDRawDataCollection_source_*_*',
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*',
        'keep L2MuonTrajectorySeeds_hltL2MuonSeeds_*_*',
        'keep L3MuonTrajectorySeeds_hltL3TrajSeedOIState_*_*',
        'keep SiPixelClusteredmNewDetSetVector_hltSiPixelClusters_*_*',
        'keep TrackingRecHitsOwned_hltL3Muons_*_*',
        'keep edmTriggerResults_*_*_*',
        'keep recoCaloJets_*_*_*',
        'keep recoCaloMETs_*_*_*',
        'keep recoCaloMETs_hltMet_*_*',
        'keep recoCompositeCandidates_*_*_*',
        'keep recoElectrons_*_*_*',
        'keep recoIsolatedPixelTrackCandidates_*_*_*',
        'keep recoMETs_*_*_*',
        'keep recoPFJets_*_*_*',
        'keep recoPFTaus_*_*_*',
        'keep recoRecoChargedCandidates_*_*_*',
        'keep recoRecoChargedCandidates_hltL2MuonCandidates_*_*',
        'keep recoRecoEcalCandidates_*_*_*',
        'keep triggerTriggerEventWithRefs_*_*_*',
        'keep triggerTriggerEvent_*_*_*',
        'keep triggerTriggerFilterObjectWithRefs_*_*_*',
        'keep *_hltFEDSelectorL1_*_*',
        'keep *_hltScoutingEgammaPacker_*_*',
        'keep *_hltScoutingMuonPacker_*_*',
        'keep *_hltScoutingPFPacker_*_*',
        'keep *_hltScoutingPrimaryVertexPacker_*_*',
        'keep *_hltScoutingTrackPacker_*_*',
        'keep edmTriggerResults_*_*_*'
    ) )
)

HLTDebugFEVT  = cms.PSet(
    outputCommands = cms.vstring( *(
        'drop *_hlt*_*_*',
        'keep *_hltAK4CaloJetsCorrectedIDPassed_*_*',
        'keep *_hltAK4CaloJetsIDPassed_*_*',
        'keep *_hltAK4CaloJets_*_*',
        'keep *_hltAK4PFJetsCorrected_*_*',
        'keep *_hltAK4PFJetsForTaus_*_*',
        'keep *_hltAK4PFJets_*_*',
        'keep *_hltAlCaEtaEBRechitsToDigis_*_*',
        'keep *_hltAlCaEtaEERechitsToDigis_*_*',
        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegional_etaEcalRecHitsES_*',
        'keep *_hltAlCaPi0EBRechitsToDigis_*_*',
        'keep *_hltAlCaPi0EERechitsToDigis_*_*',
        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegional_pi0EcalRecHitsES_*',
        'keep *_hltAlcaPixelClusterCounts_*_*',
        'keep *_hltBSoftMuonMu5L3_*_*',
        'keep *_hltCsc2DRecHits_*_*',
        'keep *_hltCscSegments_*_*',
        'keep *_hltDeepBLifetimeTagInfosPF_*_*',
        'keep *_hltDeepCombinedSecondaryVertexBJetTagsCalo_*_*',
        'keep *_hltDeepCombinedSecondaryVertexBJetTagsInfosCalo_*_*',
        'keep *_hltDeepCombinedSecondaryVertexBJetTagsInfos_*_*',
        'keep *_hltDeepCombinedSecondaryVertexBJetTagsPF_*_*',
        'keep *_hltDeepSecondaryVertexTagInfosPF_*_*',
        'keep *_hltDisplacedhltIter4PFlowTrackSelectionHighPurity_*_*',
        'keep *_hltDoubletRecoveryPFlowTrackSelectionHighPurity_*_*',
        'keep *_hltDt4DSegments_*_*',
        'keep *_hltEcalPhiSymFilter_*_*',
        'keep *_hltEcalRecHit_*_*',
        'keep *_hltEgammaCandidates_*_*',
        'keep *_hltEgammaGsfElectrons_*_*',
        'keep *_hltEgammaGsfTracks_*_*',
        'keep *_hltFEDSelectorGEM_*_*',
        'keep *_hltFEDSelectorTCDS_*_*',
        'keep *_hltFastPVPixelTracksMerger_*_*',
        'keep *_hltFastPVPixelTracksRecover_*_*',
        'keep *_hltFastPVPixelTracks_*_*',
        'keep *_hltFastPVPixelVertices_*_*',
        'keep *_hltFastPixelBLifetimeL3Associator_*_*',
        'keep *_hltFastPrimaryVertex_*_*',
        'keep *_hltGtStage2Digis_*_*',
        'keep *_hltGtStage2ObjectMap_*_*',
        'keep *_hltHbhereco_*_*',
        'keep *_hltHfreco_*_*',
        'keep *_hltHoreco_*_*',
        'keep *_hltImpactParameterTagInfos_*_*',
        'keep *_hltInclusiveSecondaryVertexFinderTagInfos_*_*',
        'keep *_hltIsolPixelTrackProdHB_*_*',
        'keep *_hltIsolPixelTrackProdHE_*_*',
        'keep *_hltIter0HighPtTkMuTrackSelectionHighPurity_*_*',
        'keep *_hltIter2MergedForDisplaced_*_*',
        'keep *_hltIterL3GlbMuon_*_*',
        'keep *_hltIterL3MuonAndMuonFromL1Merged_*_*',
        'keep *_hltIterL3MuonMerged_*_*',
        'keep *_hltIterL3MuonsNoID_*_*',
        'keep *_hltIterL3Muons_*_*',
        'keep *_hltIterL3OIMuonTrackSelectionHighPurity_*_*',
        'keep *_hltL2MuonCandidatesNoVtx_*_*',
        'keep *_hltL2MuonCandidates_*_*',
        'keep *_hltL2MuonSeeds_*_*',
        'keep *_hltL2Muons_*_*',
        'keep *_hltL2TauJets_*_*',
        'keep *_hltL3MuonsIOHit_*_*',
        'keep *_hltL3MuonsLinksCombination_*_*',
        'keep *_hltL3MuonsOIHit_*_*',
        'keep *_hltL3MuonsOIState_*_*',
        'keep *_hltL3Muons_*_*',
        'keep *_hltL3NoFiltersNoVtxMuonCandidates_*_*',
        'keep *_hltL3NoFiltersNoVtxMuons_*_*',
        'keep *_hltL3TkFromL2OICombination_*_*',
        'keep *_hltL3TkTracksFromL2IOHit_*_*',
        'keep *_hltL3TkTracksFromL2OIHit_*_*',
        'keep *_hltL3TkTracksFromL2OIState_*_*',
        'keep *_hltL3TkTracksFromL2_*_*',
        'keep *_hltL3TrackCandidateFromL2IOHit_*_*',
        'keep *_hltL3TrackCandidateFromL2OIHit_*_*',
        'keep *_hltL3TrackCandidateFromL2OIState_*_*',
        'keep *_hltL3TrajSeedIOHit_*_*',
        'keep *_hltL3TrajSeedOIHit_*_*',
        'keep *_hltL3TrajSeedOIState_*_*',
        'keep *_hltL3TrajectorySeed_*_*',
        'keep *_hltMergedTracksForBTag_*_*',
        'keep *_hltMergedTracks_*_*',
        'keep *_hltMet_*_*',
        'keep *_hltMuonCSCDigis_*_*',
        'keep *_hltMuonCSCDigis_MuonCSCStripDigi_*',
        'keep *_hltMuonCSCDigis_MuonCSCWireDigi_*',
        'keep *_hltMuonDTDigis_*_*',
        'keep *_hltMuonRPCDigis_*_*',
        'keep *_hltOnlineBeamSpot_*_*',
        'keep *_hltPFJetForBtag_*_*',
        'keep *_hltPFJetForPNetAK8_*_*',
        'keep *_hltPFMETNoMuProducer_*_*',
        'keep *_hltPFMETProducer_*_*',
        'keep *_hltPFMETTypeOne_*_*',
        'keep *_hltPFMuonMerging_*_*',
        'keep *_hltPFTau35Track_*_*',
        'keep *_hltPFTau35_*_*',
        'keep *_hltPPSCalibrationRaw_*_*',
        'keep *_hltParticleFlowForTaus_*_*',
        'keep *_hltParticleFlow_*_*',
        'keep *_hltParticleNetDiscriminatorsJetTagsAK8_*_*',
        'keep *_hltParticleNetDiscriminatorsJetTags_*_*',
        'keep *_hltPixelTracks_*_*',
        'keep *_hltPixelVertices_*_*',
        'keep *_hltRpcRecHits_*_*',
        'keep *_hltSelector4CentralJetsL1FastJet_*_*',
        'keep *_hltSelector8CentralJetsL1FastJet_*_*',
        'keep *_hltSelectorJets20L1FastJet_*_*',
        'keep *_hltSiPixelClustersCache_*_*',
        'keep *_hltSiPixelClusters_*_*',
        'keep *_hltSiStripClusterizerForRawPrime_*_*',
        'keep *_hltSiStripClusters2ApproxClusters_*_*',
        'keep *_hltSiStripRawToClustersFacility_*_*',
        'keep *_hltTowerMakerForAll_*_*',
        'keep *_hltTriggerSummaryAOD_*_*',
        'keep *_hltTriggerSummaryRAW_*_*',
        'keep *_hltTrimmedPixelVertices_*_*',
        'keep *_hltVerticesL3_*_*',
        'keep *_hltVerticesPFFilter_*_*',
        'keep *_hltVerticesPFSelector_*_*',
        'keep DetIds_hltSiStripRawToDigi_*_*',
        'keep FEDRawDataCollection_rawDataCollector_*_*',
        'keep FEDRawDataCollection_rawDataRepacker_*_*',
        'keep FEDRawDataCollection_source_*_*',
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*',
        'keep L2MuonTrajectorySeeds_hltL2MuonSeeds_*_*',
        'keep L3MuonTrajectorySeeds_hltL3TrajSeedOIState_*_*',
        'keep SiPixelClusteredmNewDetSetVector_hltSiPixelClusters_*_*',
        'keep TrackingRecHitsOwned_hltL3Muons_*_*',
        'keep edmTriggerResults_*_*_*',
        'keep recoCaloJets_*_*_*',
        'keep recoCaloMETs_*_*_*',
        'keep recoCaloMETs_hltMet_*_*',
        'keep recoCompositeCandidates_*_*_*',
        'keep recoElectrons_*_*_*',
        'keep recoIsolatedPixelTrackCandidates_*_*_*',
        'keep recoMETs_*_*_*',
        'keep recoPFJets_*_*_*',
        'keep recoPFTaus_*_*_*',
        'keep recoRecoChargedCandidates_*_*_*',
        'keep recoRecoChargedCandidates_hltL2MuonCandidates_*_*',
        'keep recoRecoEcalCandidates_*_*_*',
        'keep triggerTriggerEventWithRefs_*_*_*',
        'keep triggerTriggerEvent_*_*_*',
        'keep triggerTriggerFilterObjectWithRefs_*_*_*',
        'keep *_hltFEDSelectorL1_*_*',
        'keep *_hltScoutingEgammaPacker_*_*',
        'keep *_hltScoutingMuonPacker_*_*',
        'keep *_hltScoutingPFPacker_*_*',
        'keep *_hltScoutingPrimaryVertexPacker_*_*',
        'keep *_hltScoutingTrackPacker_*_*',
        'keep edmTriggerResults_*_*_*'
    ) )
)

HLTScouting   = cms.PSet(
    outputCommands = cms.vstring( *(
        'keep *_hltFEDSelectorL1_*_*',
        'keep *_hltScoutingEgammaPacker_*_*',
        'keep *_hltScoutingMuonPacker_*_*',
        'keep *_hltScoutingPFPacker_*_*',
        'keep *_hltScoutingPrimaryVertexPacker_*_*',
        'keep *_hltScoutingTrackPacker_*_*',
        'keep edmTriggerResults_*_*_*'
    ) )
)

