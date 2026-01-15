import FWCore.ParameterSet.Config as cms
from SimCalorimetry.HGCalAssociatorProducers.hgCalLCToSCAssociatorByEnergyScoreProducer_cfi import hgCalLCToSCAssociatorByEnergyScoreProducer as _scAssocByEnergyScoreProducer

lcAssocByEnergyScoreProducer = _scAssocByEnergyScoreProducer.clone(
    label_scl = cms.InputTag("mix", "MergedCaloTruthCaloParticle"), # CaloParticle as SimCluster dataformat
    hardScatterOnly = cms.bool(True)
)
scAssocByEnergyScoreProducer = _scAssocByEnergyScoreProducer.clone(
    label_scl = cms.InputTag("mix", "MergedCaloTruthCaloParticle"), # CaloParticle as SimCluster dataformat
    hardScatterOnly = cms.bool(True)
)

from Configuration.ProcessModifiers.enableCPfromPU_cff import enableCPfromPU

enableCPfromPU.toModify(lcAssocByEnergyScoreProducer, hardScatterOnly = cms.bool(False))
enableCPfromPU.toModify(scAssocByEnergyScoreProducer, hardScatterOnly = cms.bool(False))
