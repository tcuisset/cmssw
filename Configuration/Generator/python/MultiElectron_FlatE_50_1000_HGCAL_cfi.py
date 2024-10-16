import FWCore.ParameterSet.Config as cms

generator = cms.EDProducer("MultiParticleProducer",
    PGunParameters = cms.PSet(
        PartID = cms.vint32(11, -11),
        MaxEta = cms.double(2.9),
        MinEta = cms.double(1.6),
        
        MinE = cms.double(50.),
        MaxE = cms.double(1000),

        NParticlesInEta = cms.uint32(7),
        EtaMinSeparation = cms.double(0.3),
        EtaRandomSpread = cms.double(0.2),

        NParticlesInPhi = cms.uint32(10),
        PhiRandomSpread = cms.double(0.4),

        # ignored : 
        MinPhi = cms.double(-3.14159265359), ## in radians
        MaxPhi = cms.double(3.14159265359)

        
    ),
    Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts

    psethack = cms.string('multi electron E50to1000 HGCAL'),
    AddAntiParticle = cms.bool(True),
    firstRun = cms.untracked.uint32(1)
)
