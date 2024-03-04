import FWCore.ParameterSet.Config as cms

from ..sequences.HLTBeginSequence_cfi import *
from ..sequences.HLTEndSequence_cfi import *
from ..modules.hltPreDoublePFTauHPS_cfi import *
from ..sequences.HLTDoLocalPixelSequence_cfi import *
from ..sequences.HLTParticleFlowSequence_cfi import *

#HLTL2TauTagNNSequence = cms.Sequence(HLTDoLocalPixelSequence + fragment.HLTRecoPixelTracksSequence + fragment.HLTRecopixelvertexingSequence + fragment.HLTDoCaloSequence + cms.ignore(fragment.hltL1sDoubleTauBigOR) + cms.ignore(fragment.hltL1sSingleTau) + cms.ignore(fragment.hltL1sBigOrMuXXerIsoTauYYer) + cms.ignore(fragment.hltL1sMu22erIsoTau40er) + cms.ignore(fragment.hltL1sBigORDoubleTauJet) + cms.ignore(fragment.hltL1VBFDiJetIsoTau) + cms.ignore(fragment.hltL1sVeryBigORMu18erTauXXer2p1) + cms.ignore(fragment.hltL1sDoubleTauBigORWithLowMass) + fragment.hltL2TauTagNNProducer)

hltL2TauTagNNProducer = cms.EDProducer( "L2TauNNProducer",
    debugLevel = cms.int32( 0 ),
    L1Taus = cms.VPSet( 
      cms.PSet(  L1CollectionName = cms.string( "DoubleTau" ),
        L1TauTrigger = cms.InputTag( "hltL1sDoubleTauBigOR" )
      ),
      cms.PSet(  L1CollectionName = cms.string( "SingleTau" ),
        L1TauTrigger = cms.InputTag( "hltL1sSingleTau" )
      ),
      cms.PSet(  L1CollectionName = cms.string( "MuXXTauYY" ),
        L1TauTrigger = cms.InputTag( "hltL1sBigOrMuXXerIsoTauYYer" )
      ),
      cms.PSet(  L1CollectionName = cms.string( "Mu22Tau40" ),
        L1TauTrigger = cms.InputTag( "hltL1sMu22erIsoTau40er" )
      ),
      cms.PSet(  L1CollectionName = cms.string( "DoubleTauJet" ),
        L1TauTrigger = cms.InputTag( "hltL1sBigORDoubleTauJet" )
      ),
      cms.PSet(  L1CollectionName = cms.string( "VBFIsoTau" ),
        L1TauTrigger = cms.InputTag( "hltL1VBFDiJetIsoTau" )
      ),
      cms.PSet(  L1CollectionName = cms.string( "Mu18TauXX" ),
        L1TauTrigger = cms.InputTag( "hltL1sVeryBigORMu18erTauXXer2p1" )
      ),
      cms.PSet(  L1CollectionName = cms.string( "DoubleTauLowMass" ),
        L1TauTrigger = cms.InputTag( "hltL1sDoubleTauBigORWithLowMass" )
      )
    ),
    hbheInput = cms.InputTag( "hltHbhereco" ),
    hoInput = cms.InputTag( "hltHoreco" ),
    ebInput = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEB' ),
    eeInput = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEE' ),
    pataVertices = cms.InputTag( "hltPixelVerticesSoA" ),
    pataTracks = cms.InputTag( "hltPixelTracksSoA" ),
    BeamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    maxVtx = cms.uint32( 100 ),
    fractionSumPt2 = cms.double( 0.3 ),
    minSumPt2 = cms.double( 0.0 ),
    track_pt_min = cms.double( 1.0 ),
    track_pt_max = cms.double( 10.0 ),
    track_chi2_max = cms.double( 99999.0 ),
    graphPath = cms.string( "RecoTauTag/TrainingFiles/data/L2TauNNTag/L2TauTag_Run3v1.pb" ),
    normalizationDict = cms.string( "RecoTauTag/TrainingFiles/data/L2TauNNTag/NormalizationDict.json" )
)

hltL2DoubleTauTagNNFilter = cms.EDFilter( "L2TauTagFilter",
    saveTags = cms.bool( True ),
    nExpected = cms.int32( 2 ),
    L1TauSrc = cms.InputTag( "hltL1sDoubleTauBigOR" ),
    L2Outcomes = cms.InputTag( 'hltL2TauTagNNProducer','DoubleTau' ),
    DiscrWP = cms.double( 0.386 ),
    l1TauPtThreshold = cms.double( 250.0 )
)

HLT_DoublePNetTauh = cms.Path(
    HLTBeginSequence +
    hltPreDoublePFTauHPS + 
    HLTParticleFlowSequence +
    hltL2TauTagNNProducer + 
    hltL2DoubleTauTagNNFilter + 
    HLTEndSequence
)