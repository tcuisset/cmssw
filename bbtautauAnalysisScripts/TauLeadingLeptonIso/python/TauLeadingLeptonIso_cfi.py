import FWCore.ParameterSet.Config as cms

TauLeadingLeptonIso = cms.EDProducer('TauLeadingLeptonIso',
                                            TauCollection = cms.InputTag('slimmedTaus'),
                                            electronCollection = cms.InputTag('slimmedElectrons'),
                                            muonCollection = cms.InputTag('slimmedMuons'),
                                            verboseDebug=cms.bool(False),
                                            rhoSrc = cms.InputTag("fixedGridRhoFastjetAll"),
                                            EAConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
)
