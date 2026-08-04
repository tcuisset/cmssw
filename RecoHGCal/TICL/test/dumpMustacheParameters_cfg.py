import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.register("seedEta", 2.0, VarParsing.multiplicity.singleton, VarParsing.varType.float,
                 "Seed trackster eta")
options.register("seedPhi", 0.0, VarParsing.multiplicity.singleton, VarParsing.varType.float,
                 "Seed trackster phi")
options.register("clusterEt", 10.0, VarParsing.multiplicity.singleton, VarParsing.varType.float,
                 "Candidate trackster ET [GeV]; set <=0 when using clusterEnergy")
options.register("clusterEnergy", -1.0, VarParsing.multiplicity.singleton, VarParsing.varType.float,
                 "Fixed candidate raw energy [GeV]; set <=0 when using clusterEt")
options.register("etaHalfRange", 0.3, VarParsing.multiplicity.singleton, VarParsing.varType.float,
                 "Half-width of eta scan")
options.register("phiHalfRange", 0.6, VarParsing.multiplicity.singleton, VarParsing.varType.float,
                 "Half-width of phi scan")
options.register("nEta", 301, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                 "Number of eta grid points")
options.register("nPhi", 401, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                 "Number of phi grid points")
options.register("output", "mustache_scan.csv", VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Output CSV file")
options.parseArguments()

process = cms.Process("DUMPMUSTACHE")
process.source = cms.Source("EmptySource", firstRun=cms.untracked.uint32(1))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))

# These records are populated by GlobalTag in normal configurations.
# However there are also configurable ESProducers in that serve to prepare conditions DB upload and to force conditions locally:
# process.load("RecoEcal.EgammaCoreTools.EcalMustacheSCParametersESProducer_cff")
# process.load("RecoEcal.EgammaCoreTools.EcalSCDynamicDPhiParametersESProducer_cff")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T35', '')

process.dumpMustache = cms.EDAnalyzer(
    "MustacheESProductDumper",
    seedEta=cms.double(options.seedEta),
    seedPhi=cms.double(options.seedPhi),
    clusterEt=cms.double(options.clusterEt),
    clusterEnergy=cms.double(options.clusterEnergy),
    etaHalfRange=cms.double(options.etaHalfRange),
    phiHalfRange=cms.double(options.phiHalfRange),
    nEta=cms.uint32(options.nEta),
    nPhi=cms.uint32(options.nPhi),
    output=cms.string(options.output),
)

process.path = cms.Path(process.dumpMustache)
