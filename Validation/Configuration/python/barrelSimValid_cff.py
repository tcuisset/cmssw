import FWCore.ParameterSet.Config as cms

from SimCalorimetry.HGCalAssociatorProducers.LCToSCAssociation_cff import barrelLCToSCAssociatorByEnergyScoreProducer, barrelLayerClusterCaloParticleAssociation, barrelLayerClusterSimClusterAssociation

from Validation.HGCalValidation.BarrelValidator_cff import barrelValidator

barrelValidatorSequence = cms.Sequence(barrelValidator)

barrelAssociators = cms.Task(
    barrelLCToSCAssociatorByEnergyScoreProducer,
    barrelLayerClusterCaloParticleAssociation,
    barrelLayerClusterSimClusterAssociation
)

barrelValidation = cms.Sequence(barrelValidatorSequence)