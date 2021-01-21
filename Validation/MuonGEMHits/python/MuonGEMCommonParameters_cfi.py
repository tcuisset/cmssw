import FWCore.ParameterSet.Config as cms


GEMValidationCommonParameters = cms.PSet(
    detailPlot = cms.bool(False),
    # ZR occupancy plots
    ZROccRange = cms.untracked.vdouble(
    #   xlow, xup, ylow, yup
        525, 555, 66, 160, # station 0
        564, 574, 110, 290, # station 1
        792, 802, 120, 390), # station 2
    ZROccNumBins = cms.untracked.vint32(
    # nbinsx, nbinsy
        30, 100, # station0
        200, 150, # station1
        200, 250), # station2
    XYOccNumBins = cms.untracked.int32(720),
    EtaOccRange = cms.untracked.vdouble(1.55, 2.85),
)
