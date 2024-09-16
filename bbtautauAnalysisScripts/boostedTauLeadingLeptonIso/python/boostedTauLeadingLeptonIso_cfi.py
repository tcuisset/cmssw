import FWCore.ParameterSet.Config as cms

boostedTauLeadingLeptonIso = cms.EDProducer('boostedTauLeadingLeptonIso',
                                            boostedTauCollection = cms.InputTag('slimmedTausBoosted'),
                                            electronCollection = cms.InputTag('slimmedElectrons'),
                                            muonCollection = cms.InputTag('slimmedMuons'),
                                            verboseDebug=cms.bool(False),
                                            rhoSrc = cms.InputTag("fixedGridRhoFastjetAll"),
                                            EAConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
)
