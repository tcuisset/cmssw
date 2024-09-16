import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.nano_eras_cff import *
from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer


#Ganesh Changes for the private production of NanoAOD
from bbtautauAnalysisScripts.boostedTauLeadingLeptonIso.boostedTauLeadingLeptonIso_cfi import *

boostedTauLeadingLeptonIso.boostedTauCollection  = cms.InputTag("slimmedTausBoostedNewID")
boostedTauLeadingLeptonIso.muonCollection  = cms.InputTag("slimmedMuonsUpdated")
from RecoTauTag.RecoTau.tauIdWPsDefs import WORKING_POINTS_v2p7
#Import from the user defined plugins


##################### Import reusable funtions and objects from std taus ########
from PhysicsTools.NanoAOD.taus_cff import _tauIdWPMask, tausMCMatchLepTauForTable, tausMCMatchHadTauForTable,tauMCTable

##################### User floats producers, selectors ##########################


finalBoostedTaus = cms.EDFilter("PATTauRefSelector",
    #Ganesh
    src = cms.InputTag("slimmedboostedTauWithUserData"),
    #cut = cms.string("pt > 40 && tauID('decayModeFindingNewDMs') && (tauID('byVVLooseIsolationMVArun2DBoldDMwLT') || tauID('byVVLooseIsolationMVArun2DBnewDMwLT'))")
    cut = cms.string("pt > 18")
)
run2_nanoAOD_106Xv2.toModify(
    finalBoostedTaus,
    src = cms.InputTag("slimmedboostedTauWithUserData"),
    #cut = "pt > 40 && tauID('decayModeFindingNewDMs') && (tauID('byVVLooseIsolationMVArun2017v2DBoldDMwLT2017') || tauID('byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017') || tauID('byVVLooseIsolationMVArun2017v2DBnewDMwLT2017'))"
    cut = "pt > 18"
)

#Changes by Ganesh for Custom nanoAOD production

slimmedboostedTauWithUserData = cms.EDProducer("PATTauUserDataEmbedder",
     src = cms.InputTag("slimmedTausBoostedNewID"),
     userFloats = cms.PSet(

        LeadingElectrondelR = cms.InputTag("boostedTauLeadingLeptonIso:LeadingElectrondelR"),
        SubLeadingElectrondelR = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingElectrondelR"),
        SubSubLeadingElectrondelR = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingElectrondelR"),
        LeadingMuondelR = cms.InputTag("boostedTauLeadingLeptonIso:LeadingMuondelR"),
        SubLeadingMuondelR = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingMuondelR"),
        SubSubLeadingMuondelR = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingMuondelR"),



        LeadingElectronPt = cms.InputTag("boostedTauLeadingLeptonIso:LeadingElectronPt"),
        LeadingElectronEta = cms.InputTag("boostedTauLeadingLeptonIso:LeadingElectronEta"),
        LeadingElectronPhi = cms.InputTag("boostedTauLeadingLeptonIso:LeadingElectronPhi"),
        LeadingElectronM = cms.InputTag("boostedTauLeadingLeptonIso:LeadingElectronM"),
        LeadingElectronCorrIso = cms.InputTag("boostedTauLeadingLeptonIso:LeadingElectronCorrIso"),
        LeadingElectronsumPFChargedHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:LeadingElectronsumPFChargedHadronPt"),
        LeadingElectronsumPFNeutralHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:LeadingElectronsumPFNeutralHadronPt"),
        LeadingElectronsumPFPhotonPt = cms.InputTag("boostedTauLeadingLeptonIso:LeadingElectronsumPFPhotonPt"),
        LeadingElectronea = cms.InputTag("boostedTauLeadingLeptonIso:LeadingElectronea"),
        LeadingElectronrho = cms.InputTag("boostedTauLeadingLeptonIso:LeadingElectronrho"),
        LeadingElectrontausumPFChargedHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:LeadingElectrontausumPFChargedHadronPt"),
        LeadingElectrontausumPFNeutralHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:LeadingElectrontausumPFNeutralHadronPt"),
        LeadingElectrontausumPFPhotonPt = cms.InputTag("boostedTauLeadingLeptonIso:LeadingElectrontausumPFPhotonPt"),


        SubLeadingElectronPt = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingElectronPt"),
        SubLeadingElectronEta = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingElectronEta"),
        SubLeadingElectronPhi = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingElectronPhi"),
        SubLeadingElectronM = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingElectronM"),
        SubLeadingElectronCorrIso = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingElectronCorrIso"),
        SubLeadingElectronsumPFChargedHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingElectronsumPFChargedHadronPt"),
        SubLeadingElectronsumPFNeutralHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingElectronsumPFNeutralHadronPt"),
        SubLeadingElectronsumPFPhotonPt = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingElectronsumPFPhotonPt"),
        SubLeadingElectronea = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingElectronea"),
        SubLeadingElectronrho = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingElectronrho"),
        SubLeadingElectrontausumPFChargedHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingElectrontausumPFChargedHadronPt"),
        SubLeadingElectrontausumPFNeutralHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingElectrontausumPFNeutralHadronPt"),
        SubLeadingElectrontausumPFPhotonPt = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingElectrontausumPFPhotonPt"),  

        SubSubLeadingElectronPt = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingElectronPt"),
        SubSubLeadingElectronEta = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingElectronEta"),
        SubSubLeadingElectronPhi = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingElectronPhi"),
        SubSubLeadingElectronM = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingElectronM"),
        SubSubLeadingElectronCorrIso = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingElectronCorrIso"),
        SubSubLeadingElectronsumPFChargedHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingElectronsumPFChargedHadronPt"),
        SubSubLeadingElectronsumPFNeutralHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingElectronsumPFNeutralHadronPt"),
        SubSubLeadingElectronsumPFPhotonPt = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingElectronsumPFPhotonPt"),
        SubSubLeadingElectronea = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingElectronea"),
        SubSubLeadingElectronrho = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingElectronrho"),
        SubSubLeadingElectrontausumPFChargedHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingElectrontausumPFChargedHadronPt"),
        SubSubLeadingElectrontausumPFNeutralHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingElectrontausumPFNeutralHadronPt"),
        SubSubLeadingElectrontausumPFPhotonPt = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingElectrontausumPFPhotonPt"),  
        ##############################################################################


        LeadingMuonPt = cms.InputTag("boostedTauLeadingLeptonIso:LeadingMuonPt"),
        LeadingMuonEta = cms.InputTag("boostedTauLeadingLeptonIso:LeadingMuonEta"),
        LeadingMuonPhi = cms.InputTag("boostedTauLeadingLeptonIso:LeadingMuonPhi"),
        LeadingMuonM = cms.InputTag("boostedTauLeadingLeptonIso:LeadingMuonM"),
        LeadingMuonCorrIso = cms.InputTag("boostedTauLeadingLeptonIso:LeadingMuonCorrIso"),
        LeadingMuonsumPFChargedHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:LeadingMuonsumPFChargedHadronPt"),
        LeadingMuonsumPFNeutralHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:LeadingMuonsumPFNeutralHadronPt"),
        LeadingMuonsumPFPhotonPt = cms.InputTag("boostedTauLeadingLeptonIso:LeadingMuonsumPFPhotonPt"),
        LeadingMuonsumPUPt = cms.InputTag("boostedTauLeadingLeptonIso:LeadingMuonsumPUPt"),
        LeadingMuontausumPFChargedHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:LeadingMuontausumPFChargedHadronPt"),
        LeadingMuontausumPFNeutralHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:LeadingMuontausumPFNeutralHadronPt"),
        LeadingMuontausumPFPhotonPt = cms.InputTag("boostedTauLeadingLeptonIso:LeadingMuontausumPFPhotonPt"),

        SubLeadingMuonPt = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingMuonPt"),
        SubLeadingMuonEta = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingMuonEta"),
        SubLeadingMuonPhi = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingMuonPhi"),
        SubLeadingMuonM = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingMuonM"),
        SubLeadingMuonCorrIso = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingMuonCorrIso"),
        SubLeadingMuonsumPFChargedHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingMuonsumPFChargedHadronPt"),
        SubLeadingMuonsumPFNeutralHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingMuonsumPFNeutralHadronPt"),
        SubLeadingMuonsumPFPhotonPt = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingMuonsumPFPhotonPt"),
        SubLeadingMuonsumPUPt = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingMuonsumPUPt"),
        SubLeadingMuontausumPFChargedHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingMuontausumPFChargedHadronPt"),
        SubLeadingMuontausumPFNeutralHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingMuontausumPFNeutralHadronPt"),
        SubLeadingMuontausumPFPhotonPt = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingMuontausumPFPhotonPt"),        

        SubSubLeadingMuonPt = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingMuonPt"),
        SubSubLeadingMuonEta = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingMuonEta"),
        SubSubLeadingMuonPhi = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingMuonPhi"),
        SubSubLeadingMuonM = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingMuonM"),
        SubSubLeadingMuonCorrIso = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingMuonCorrIso"),
        SubSubLeadingMuonsumPFChargedHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingMuonsumPFChargedHadronPt"),
        SubSubLeadingMuonsumPFNeutralHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingMuonsumPFNeutralHadronPt"),
        SubSubLeadingMuonsumPFPhotonPt = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingMuonsumPFPhotonPt"),
        SubSubLeadingMuonsumPUPt = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingMuonsumPUPt"),
        SubSubLeadingMuontausumPFChargedHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingMuontausumPFChargedHadronPt"),
        SubSubLeadingMuontausumPFNeutralHadronPt = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingMuontausumPFNeutralHadronPt"),
        SubSubLeadingMuontausumPFPhotonPt = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingMuontausumPFPhotonPt"),         
  
     ),
     userInts = cms.PSet(
      Ecounter = cms.InputTag("boostedTauLeadingLeptonIso:Ecounter"),
      Mcounter = cms.InputTag("boostedTauLeadingLeptonIso:Mcounter"),

      LeadingElectron_electronIdx = cms.InputTag("boostedTauLeadingLeptonIso:LeadingElectronelectronIdx"),
     SubLeadingElectron_electronIdx = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingElectronelectronIdx"),      
     SubSubLeadingElectron_electronIdx = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingElectronelectronIdx"),

      LeadingMuon_muonIdx = cms.InputTag("boostedTauLeadingLeptonIso:LeadingMuonmuonIdx"),
     SubLeadingMuon_muonIdx = cms.InputTag("boostedTauLeadingLeptonIso:SubLeadingMuonmuonIdx"),      
     SubSubLeadingMuon_muonIdx = cms.InputTag("boostedTauLeadingLeptonIso:SubSubLeadingMuonmuonIdx"),
     ),

)

boostedTauTable = simpleCandidateFlatTableProducer.clone(
    src = cms.InputTag("linkedObjects", "boostedTaus"),
    name= cms.string("boostedTau"),
    doc = cms.string("slimmedBoostedTaus after basic selection (" + finalBoostedTaus.cut.value()+")"),
    variables = cms.PSet() # PSet defined below in era dependent way
)
_boostedTauVarsBase = cms.PSet(P4Vars,
       charge = Var("charge", int, doc="electric charge"),
       jetIdx = Var("?hasUserCand('jet')?userCand('jet').key():-1", "int16", doc="index of the associated jet (-1 if none)"),
       decayMode = Var("decayMode()",int),
       idDecayModeOldDMs = Var("(?isTauIDAvailable('decayModeFinding')?tauID('decayModeFinding'):-1) > 0", bool),
       idDecayModeNewDMs = Var("(?isTauIDAvailable('decayModeFindingNewDMs')?tauID('decayModeFindingNewDMs'):-1) > 0", bool),
       leadTkPtOverTauPt = Var("leadChargedHadrCand.pt/pt ",float, doc="pt of the leading track divided by tau pt",precision=10),
       leadTkDeltaEta = Var("leadChargedHadrCand.eta - eta ",float, doc="eta of the leading track, minus tau eta",precision=8),
       leadTkDeltaPhi = Var("deltaPhi(leadChargedHadrCand.phi, phi) ",float, doc="phi of the leading track, minus tau phi",precision=8),


       #variables added by Ganesh for Electron-Match

       #add DelR
       LeadingElectrondelR = Var("userFloat('LeadingElectrondelR')",float,doc="Leading Matched Electron delR"),
       SubLeadingElectrondelR = Var("userFloat('SubLeadingElectrondelR')",float,doc="Sub Leading Matched Electron delR"),      
       SubSubLeadingElectrondelR =  Var("userFloat('SubSubLeadingElectrondelR')",float,doc="Sub Sub Leading Matched Electron delR"),
       LeadingMuondelR = Var("userFloat('LeadingMuondelR')",float,doc="Leading Matched Muon delR"),
       SubLeadingMuondelR = Var("userFloat('SubLeadingMuondelR')",float,doc="Sub Leading Matched Muon delR"),
       SubSubLeadingMuondelR = Var("userFloat('SubSubLeadingMuondelR')",float,doc="Sub-Sub Leading Matched Muon delR"),

       Ecounter = Var("userInt('Ecounter')",int,doc="Number of electrons that passed & matched with Tau and has the Loose ID and the delta R < 0.4 and > 0.02 requirements. Sadly we only store 3 leading"),

       LeadingElectronPt = Var("userFloat('LeadingElectronPt')",float,doc="Leading Matched Electron Pt"),
       LeadingElectronEta = Var("userFloat('LeadingElectronEta')",float,doc="Leading Matched Electron eta"),
       LeadingElectronPhi = Var("userFloat('LeadingElectronPhi')",float,doc="Leading Matched Electron phi"),
       LeadingElectronM = Var("userFloat('LeadingElectronM')",float,doc="Leading Matched Electron mass"),
       LeadingElectronCorrIso = Var("userFloat('LeadingElectronCorrIso')",float,doc="Corrected isolation for the Leading Electron matched"),
       LeadingElectronsumPFChargedHadronPt = Var("userFloat('LeadingElectronsumPFChargedHadronPt')",float,doc="sumPFChargedHadronPt for the Leading Electron matched"),
       LeadingElectronsumPFNeutralHadronPt = Var("userFloat('LeadingElectronsumPFNeutralHadronPt')",float,doc="sumPFNeutralHadronPt for the Leading Electron matched"),
       LeadingElectronsumPFPhotonPt = Var("userFloat('LeadingElectronsumPFPhotonPt')",float,doc="sumPFPhotonPt for the Leading Electron matched"),
       LeadingElectronea = Var("userFloat('LeadingElectronea')",float,doc="Electronea for the Leading Electron matched"),
       LeadingElectronrho = Var("userFloat('LeadingElectronrho')",float,doc="Electronrho for the Leading Electron matched"),
       LeadingElectrontausumPFChargedHadronPt = Var("userFloat('LeadingElectrontausumPFChargedHadronPt')",float,doc="tausumPFChargedHadronPt for the Leading Electron matched"),
       LeadingElectrontausumPFNeutralHadronPt = Var("userFloat('LeadingElectrontausumPFNeutralHadronPt')",float,doc="tausumPFNeutralHadronPt for the Leading Electron matched"),
       LeadingElectrontausumPFPhotonPt = Var("userFloat('LeadingElectrontausumPFPhotonPt')",float,doc="tausumPFPhotonPt for the Leading Electron matched"),
       LeadingElectron_electronIdx = Var("userInt('LeadingElectron_electronIdx')","int16",doc="index of the associated electron (-1 if none)"), 

       SubLeadingElectronPt = Var("userFloat('SubLeadingElectronPt')",float,doc="Sub Leading Matched Electron Pt"),
       SubLeadingElectronEta = Var("userFloat('SubLeadingElectronEta')",float,doc="Sub Leading Matched Electron eta"),
       SubLeadingElectronPhi = Var("userFloat('SubLeadingElectronPhi')",int,doc="Sub Leading Matched Electron phi"),
       SubLeadingElectronM = Var("userFloat('SubLeadingElectronM')",float,doc="Sub Leading Matched Electron mass"),
       SubLeadingElectronCorrIso = Var("userFloat('SubLeadingElectronCorrIso')",float,doc="Corrected isolation for the sub Leading Electron matched"),
       SubLeadingElectronsumPFChargedHadronPt = Var("userFloat('SubLeadingElectronsumPFChargedHadronPt')",float,doc="sumPFChargedHadronPt for the Sub Leading Electron matched"),
       SubLeadingElectronsumPFNeutralHadronPt = Var("userFloat('SubLeadingElectronsumPFNeutralHadronPt')",float,doc="sumPFNeutralHadronPt for the Sub Leading Electron matched"),
       SubLeadingElectronsumPFPhotonPt = Var("userFloat('SubLeadingElectronsumPFPhotonPt')",float,doc="sumPFPhotonPt for the Sub Leading Electron matched"),
       SubLeadingElectronea = Var("userFloat('SubLeadingElectronea')",float,doc="Electronea for the Sub Leading Electron matched"),
       SubLeadingElectronrho = Var("userFloat('SubLeadingElectronrho')",float,doc="Electronrho for the Sub Leading Electron matched"),
       SubLeadingElectrontausumPFChargedHadronPt = Var("userFloat('SubLeadingElectrontausumPFChargedHadronPt')",float,doc="tausumPFChargedHadronPt for the Sub Leading Electron matched"),
       SubLeadingElectrontausumPFNeutralHadronPt = Var("userFloat('SubLeadingElectrontausumPFNeutralHadronPt')",float,doc="tausumPFNeutralHadronPt for the Sub Leading Electron matched"),
       SubLeadingElectrontausumPFPhotonPt = Var("userFloat('SubLeadingElectrontausumPFPhotonPt')",float,doc="tausumPFPhotonPt for the Sub Leading Electron matched"),
       SubLeadingElectron_electronIdx = Var("userInt('SubLeadingElectron_electronIdx')","int16",doc="index of the associated electron (-1 if none)"), 

       SubSubLeadingElectronPt =  Var("userFloat('SubSubLeadingElectronPt')",float,doc="Sub Sub Leading Matched Electron Pt"),
       SubSubLeadingElectronEta = Var("userFloat('SubSubLeadingElectronEta')",float,doc="Sub Sub Leading Matched Electron eta"), 
       SubSubLeadingElectronPhi = Var("userFloat('SubSubLeadingElectronPhi')",int,doc="Sub Sub Leading Matched Electron phi"),
       SubSubLeadingElectronM = Var("userFloat('SubSubLeadingElectronM')",float,doc="Sub Sub Leading Matched Electron mass"),
       SubSubLeadingElectronCorrIso = Var("userFloat('SubSubLeadingElectronCorrIso')",float,doc="Corrected isolation for the sub sub Leading Electron matched"),
       SubSubLeadingElectronsumPFChargedHadronPt = Var("userFloat('SubSubLeadingElectronsumPFChargedHadronPt')",float,doc="sumPFChargedHadronPt for the Sub Sub Leading Electron matched"),
       SubSubLeadingElectronsumPFNeutralHadronPt = Var("userFloat('SubSubLeadingElectronsumPFNeutralHadronPt')",float,doc="sumPFNeutralHadronPt for the Sub Sub Leading Electron matched"),
       SubSubLeadingElectronsumPFPhotonPt = Var("userFloat('SubSubLeadingElectronsumPFPhotonPt')",float,doc="sumPFPhotonPt for the Sub Sub Leading Electron matched"),
       SubSubLeadingElectronea = Var("userFloat('SubSubLeadingElectronea')",float,doc="Electronea for the Sub Sub Leading Electron matched"),
       SubSubLeadingElectronrho = Var("userFloat('SubSubLeadingElectronrho')",float,doc="Electronrho for the Sub Sub Leading Electron matched"),
       SubSubLeadingElectrontausumPFChargedHadronPt = Var("userFloat('SubSubLeadingElectrontausumPFChargedHadronPt')",float,doc="tausumPFChargedHadronPt for the Sub Sub Leading Electron matched"),
       SubSubLeadingElectrontausumPFNeutralHadronPt = Var("userFloat('SubSubLeadingElectrontausumPFNeutralHadronPt')",float,doc="tausumPFNeutralHadronPt for the Sub Sub Leading Electron matched"),
       SubSubLeadingElectrontausumPFPhotonPt = Var("userFloat('SubSubLeadingElectrontausumPFPhotonPt')",float,doc="tausumPFPhotonPt for the Sub Sub Leading Electron matched"),      
       SubSubLeadingElectron_electronIdx = Var("userInt('SubSubLeadingElectron_electronIdx')","int16",doc="index of the associated electron (-1 if none)"),         

       #variables added by Ganesh for Muon Match
       Mcounter = Var("userInt('Mcounter')",int,doc="Number of muons that passed and matched with the Tau and has the Loose ID and the delta R < 0.4 and > 0.02 requirements. But sadly we only store 3"),

       LeadingMuonPt = Var("userFloat('LeadingMuonPt')",float,doc="Leading Matched Muon Pt"),
       LeadingMuonEta = Var("userFloat('LeadingMuonEta')",float,doc="Leading Matched Muon eta"),
       LeadingMuonPhi = Var("userFloat('LeadingMuonPhi')",float,doc="Leading Matched Muon phi"),
       LeadingMuonM = Var("userFloat('LeadingMuonM')",float,doc="Leading Matched Muon mass"),
       LeadingMuonCorrIso = Var("userFloat('LeadingMuonCorrIso')",float,doc="Corrected isolation for the Muon Leading matched"),
       LeadingMuonsumPFChargedHadronPt = Var("userFloat('LeadingMuonsumPFChargedHadronPt')",float,doc="sumPFChargedHadronPt for the Muon Leading matched"),
       LeadingMuonsumPFNeutralHadronPt = Var("userFloat('LeadingMuonsumPFNeutralHadronPt')",float,doc="sumPFNeutralHadronPt for the Muon Leading matched"),
       LeadingMuonsumPFPhotonPt = Var("userFloat('LeadingMuonsumPFPhotonPt')",float,doc="sumPFPhotonPt for the Muon Leading matched"),
       LeadingMuonsumPUPt = Var("userFloat('LeadingMuonsumPUPt')",float,doc="sumPUPt for the Muon Leading matched"),
       LeadingMuontausumPFChargedHadronPt = Var("userFloat('LeadingMuontausumPFChargedHadronPt')",float,doc="tausumPFChargedHadronPt for the Muon Leading matched"),
       LeadingMuontausumPFNeutralHadronPt = Var("userFloat('LeadingMuontausumPFNeutralHadronPt')",float,doc="tausumPFNeutralHadronPt for the Muon Leading matched"),
       LeadingMuontausumPFPhotonPt = Var("userFloat('LeadingMuontausumPFPhotonPt')",float,doc="tausumPFPhotonPt for the Muon Leading matched"),
       LeadingMuon_muonIdx = Var("userInt('LeadingMuon_muonIdx')","int16",doc="index of the associated muon (-1 if none)"), 

       SubLeadingMuonPt = Var("userFloat('SubLeadingMuonPt')",float,doc="Sub Leading Matched Muon Pt"),
       SubLeadingMuonEta = Var("userFloat('SubLeadingMuonEta')",float,doc="Sub Leading Matched Muon eta"),
       SubLeadingMuonPhi = Var("userFloat('SubLeadingMuonPhi')",float,doc="Sub Leading Matched Muon phi"),
       SubLeadingMuonM = Var("userFloat('SubLeadingMuonM')",int,doc="Sub Leading Matched Muon mass"),
       SubLeadingMuonCorrIso = Var("userFloat('SubSubLeadingMuonCorrIso')",float,doc="Corrected isolation for the sub Leading Muon matched"),
       SubLeadingMuonsumPFChargedHadronPt = Var("userFloat('SubLeadingMuonsumPFChargedHadronPt')",float,doc="sumPFChargedHadronPt for the Muon Sub Leading matched"),
       SubLeadingMuonsumPFNeutralHadronPt = Var("userFloat('SubLeadingMuonsumPFNeutralHadronPt')",float,doc="sumPFNeutralHadronPt for the Muon Sub Leading matched"),
       SubLeadingMuonsumPFPhotonPt = Var("userFloat('SubLeadingMuonsumPFPhotonPt')",float,doc="sumPFPhotonPt for the Muon Sub Leading matched"),
       SubLeadingMuonsumPUPt = Var("userFloat('SubLeadingMuonsumPUPt')",float,doc="sumPUPt for the Muon Sub Leading matched"),
       SubLeadingMuontausumPFChargedHadronPt = Var("userFloat('SubLeadingMuontausumPFChargedHadronPt')",float,doc="tausumPFChargedHadronPt for the Muon Sub Leading matched"),
       SubLeadingMuontausumPFNeutralHadronPt = Var("userFloat('SubLeadingMuontausumPFNeutralHadronPt')",float,doc="tausumPFNeutralHadronPt for the Muon Sub Leading matched"),
       SubLeadingMuontausumPFPhotonPt = Var("userFloat('SubLeadingMuontausumPFPhotonPt')",float,doc="tausumPFPhotonPt for the Muon Sub Leading matched"),
       SubLeadingMuon_muonIdx = Var("userInt('SubLeadingMuon_muonIdx')","int16",doc="index of the associated muon (-1 if none)"), 

       SubSubLeadingMuonPt = Var("userFloat('SubSubLeadingMuonPt')",float,doc="Sub-Sub Leading Matched Muon Pt"),
       SubSubLeadingMuonEta = Var("userFloat('SubSubLeadingMuonEta')",float,doc="Sub-Sub Leading Matched Muon eta"),
       SubSubLeadingMuonPhi = Var("userFloat('SubSubLeadingMuonPhi')",float,doc="Sub-Sub Leading Matched Muon phi"),
       SubSubLeadingMuonM = Var("userFloat('SubSubLeadingMuonM')",int,doc="Sub-Sub Leading Matched Muon mass"),
       SubSubLeadingMuonCorrIso = Var("userFloat('SubSubLeadingMuonCorrIso')",float,doc="Corrected isolation for the sub-sub Leading Muon matched"),
       SubSubLeadingMuonsumPFChargedHadronPt = Var("userFloat('SubSubLeadingMuonsumPFChargedHadronPt')",float,doc="sumPFChargedHadronPt for the Muon Sub Sub Leading matched"),
       SubSubLeadingMuonsumPFNeutralHadronPt = Var("userFloat('SubSubLeadingMuonsumPFNeutralHadronPt')",float,doc="sumPFNeutralHadronPt for the Muon Sub Sub Leading matched"),
       SubSubLeadingMuonsumPFPhotonPt = Var("userFloat('SubSubLeadingMuonsumPFPhotonPt')",float,doc="sumPFPhotonPt for the Muon Sub Sub Leading matched"),
       SubSubLeadingMuonsumPUPt = Var("userFloat('SubSubLeadingMuonsumPUPt')",float,doc="sumPUPt for the Muon Sub Sub Leading matched"),
       SubSubLeadingMuontausumPFChargedHadronPt = Var("userFloat('SubSubLeadingMuontausumPFChargedHadronPt')",float,doc="tausumPFChargedHadronPt for the Muon Sub Sub Leading matched"),
       SubSubLeadingMuontausumPFNeutralHadronPt = Var("userFloat('SubSubLeadingMuontausumPFNeutralHadronPt')",float,doc="tausumPFNeutralHadronPt for the Muon Sub Sub Leading matched"),
       SubSubLeadingMuontausumPFPhotonPt = Var("userFloat('SubSubLeadingMuontausumPFPhotonPt')",float,doc="tausumPFPhotonPt for the Muon Sub Sub Leading matched"),
       SubSubLeadingMuon_muonIdx = Var("userInt('SubSubLeadingMuon_muonIdx')","int16",doc="index of the associated muon (-1 if none)"),        


       rawIso = Var( "tauID('byCombinedIsolationDeltaBetaCorrRaw3Hits')", float, doc = "combined isolation (deltaBeta corrections)", precision=10),
       rawIsodR03 = Var( "(tauID('chargedIsoPtSumdR03')+max(0.,tauID('neutralIsoPtSumdR03')-0.072*tauID('puCorrPtSum')))", float, doc = "combined isolation (deltaBeta corrections, dR=0.3)", precision=10),
       chargedIso = Var( "tauID('chargedIsoPtSum')", float, doc = "charged isolation", precision=10),
       neutralIso = Var( "tauID('neutralIsoPtSum')", float, doc = "neutral (photon) isolation", precision=10),
       puCorr = Var( "tauID('puCorrPtSum')", float, doc = "pileup correction", precision=10),
       photonsOutsideSignalCone = Var( "tauID('photonPtSumOutsideSignalCone')", float, doc = "sum of photons outside signal cone", precision=10),
       idAntiMu = _tauIdWPMask("againstMuon%s3", choices=("Loose","Tight"), doc= "Anti-muon discriminator V3: ")
)
#MVA 2017 v2 variables
_boostedTauVarsMVAIso = cms.PSet(
       rawMVAoldDM2017v2 = Var("tauID('byIsolationMVArun2DBoldDMwLTraw')",float, doc="byIsolationMVArun2DBoldDMwLT raw output discriminator (2017v2)",precision=10),
       rawMVAnewDM2017v2 = Var("tauID('byIsolationMVArun2DBnewDMwLTraw')",float,doc='byIsolationMVArun2DBnewDMwLT raw output discriminator (2017v2)',precision=10),
       idMVAnewDM2017v2 = _tauIdWPMask("by%sIsolationMVArun2DBnewDMwLT", choices=("VVLoose","VLoose","Loose","Medium","Tight","VTight","VVTight"), doc="IsolationMVArun2DBnewDMwLT ID working point (2017v2)"),
       idMVAoldDM2017v2 = _tauIdWPMask("by%sIsolationMVArun2DBoldDMwLT",choices=("VVLoose","VLoose","Loose","Medium","Tight","VTight","VVTight"), doc="IsolationMVArun2DBoldDMwLT ID working point (2017v2)"),
)
#MVA 2017 v2 dR<0.3 variables
_boostedTauVarsMVAIsoDr03 = cms.PSet(
       rawMVAoldDMdR032017v2 = Var("tauID('byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017')",float,doc='byIsolationMVArun2DBoldDMdR0p3wLT raw output discriminator (2017v2)'),
       idMVAoldDMdR032017v2 = _tauIdWPMask("by%sIsolationMVArun2017v2DBoldDMdR0p3wLT2017",choices=("VVLoose","VLoose","Loose","Medium","Tight","VTight","VVTight"),doc="IsolationMVArun2DBoldDMdR0p3wLT ID working point (2017v2)")
)
#AntiEle MVA 2018 variables
_boostedTauVarsAntiEleMVA = cms.PSet(
       rawAntiEle2018 = Var("tauID('againstElectronMVA6Raw')", float, doc= "Anti-electron MVA discriminator V6 raw output discriminator (2018)", precision=10),
       rawAntiEleCat2018 = Var("tauID('againstElectronMVA6category')", "int16", doc="Anti-electron MVA discriminator V6 category (2018)"),
       idAntiEle2018 = _tauIdWPMask("againstElectron%sMVA6", choices=("VLoose","Loose","Medium","Tight","VTight"), doc= "Anti-electron MVA discriminator V6 (2018)")
)
#DeepBoostedTau ID variables added by Ganesh
_deepTauVars2018v2p7 = cms.PSet(
    rawDeepTau2018v2p7VSe = Var("tauID('byDeepTau2018v2p7VSeraw')", float, doc="byDeepTau2018v2p7VSe for boostedTaus raw output discriminator (deepTau2018v2p7)", precision=10),
    rawDeepTau2018v2p7VSmu = Var("tauID('byDeepTau2018v2p7VSmuraw')", float, doc="byDeepTau2018v2p7VSmu for boostedTaus raw output discriminator (deepTau2018v2p7)", precision=10),
    rawDeepTau2018v2p7VSjet = Var("tauID('byDeepTau2018v2p7VSjetraw')", float, doc="byDeepTau2018v2p7VSjet for boostedTaus raw output discriminator (deepTau2018v2p7)", precision=10),
    idDeepTau2018v2p7VSe = _tauIdWPMask("by%sDeepTau2018v2p7VSe",
                                            choices=("VVVLoose","VVLoose","VLoose","Loose","Medium","Tight","VTight","VVTight"),
                                            doc="byDeepTau2018v2p7VSe for boostedTaus ID working points (deepTau2018v2p7)"),
    idDeepTau2018v2p7VSmu = _tauIdWPMask("by%sDeepTau2018v2p7VSmu",
                                            choices=("VLoose", "Loose", "Medium", "Tight"),
                                            doc="byDeepTau2018v2p7VSmu for boostedTaus ID working points (deepTau2018v2p7)"),
    idDeepTau2018v2p7VSjet = _tauIdWPMask("by%sDeepTau2018v2p7VSjet",
                                            choices=("VLoose","Loose","Medium","Tight"),
                                            doc="byDeepTau2018v2p7VSjet for boostedTaus ID working points (deepTau2018v2p7)"),
)

#Ganesh Changes
boostedTauTable.variables = cms.PSet(
    _boostedTauVarsBase,
    _boostedTauVarsMVAIso,
    _boostedTauVarsAntiEleMVA,
    _deepTauVars2018v2p7
)

#Ganesh Changes
_boostedTauVarsWithDr03 = cms.PSet(
    _boostedTauVarsBase,
    _boostedTauVarsMVAIso,
    _boostedTauVarsMVAIsoDr03,
    _boostedTauVarsAntiEleMVA,
    _deepTauVars2018v2p7
)
run2_nanoAOD_106Xv2.toModify(
    boostedTauTable,
    variables = _boostedTauVarsWithDr03
).toModify(
    boostedTauTable.variables,
    rawMVAoldDM2017v2 = Var("tauID('byIsolationMVArun2017v2DBoldDMwLTraw2017')",float, doc="byIsolationMVArun2DBoldDMwLT raw output discriminator (2017v2)",precision=10),
    rawMVAnewDM2017v2 = Var("tauID('byIsolationMVArun2017v2DBnewDMwLTraw2017')",float,doc='byIsolationMVArun2DBnewDMwLT raw output discriminator (2017v2)',precision=10),
    idMVAnewDM2017v2 = _tauIdWPMask("by%sIsolationMVArun2017v2DBnewDMwLT2017", choices=("VVLoose","VLoose","Loose","Medium","Tight","VTight","VVTight"),doc="IsolationMVArun2DBnewDMwLT ID working point (2017v2)"),
    idMVAoldDM2017v2 = _tauIdWPMask("by%sIsolationMVArun2017v2DBoldDMwLT2017",choices=("VVLoose","VLoose","Loose","Medium","Tight","VTight","VVTight"),doc="IsolationMVArun2DBoldDMwLT ID working point (2017v2)"),
    rawAntiEle2018 = Var("tauID('againstElectronMVA6Raw2018')", float, doc= "Anti-electron MVA discriminator V6 raw output discriminator (2018)", precision=10),
    rawAntiEleCat2018 = Var("tauID('againstElectronMVA6category2018')", "int16", doc="Anti-electron MVA discriminator V6 category (2018)"),
    idAntiEle2018 = _tauIdWPMask("againstElectron%sMVA62018", choices=("VLoose","Loose","Medium","Tight","VTight"), doc= "Anti-electron MVA discriminator V6 (2018)")
)

boostedTausMCMatchLepTauForTable = tausMCMatchLepTauForTable.clone(
    src = boostedTauTable.src
)

#This requires genVisTaus in taus_cff.py
boostedTausMCMatchHadTauForTable = tausMCMatchHadTauForTable.clone(
    src = boostedTauTable.src
)

boostedTauMCTable = tauMCTable.clone(
    src = boostedTauTable.src,
    mcMap = cms.InputTag("boostedTausMCMatchLepTauForTable"),
    mcMapVisTau = cms.InputTag("boostedTausMCMatchHadTauForTable"),
    objName = boostedTauTable.name,
)

#Ganesh Changes
#boostedTauTask = cms.Task(finalBoostedTaus)
boostedTauTask = cms.Task(boostedTauLeadingLeptonIso,slimmedboostedTauWithUserData,finalBoostedTaus)
boostedTauTablesTask = cms.Task(boostedTauTable)
boostedTauMCTask = cms.Task(boostedTausMCMatchLepTauForTable,boostedTausMCMatchHadTauForTable,boostedTauMCTable)

#remove boosted tau from previous eras
(run3_nanoAOD_122).toReplaceWith(
    boostedTauTask,cms.Task()
).toReplaceWith(
    boostedTauTablesTask,cms.Task()
).toReplaceWith(
    boostedTauMCTask,cms.Task()
)