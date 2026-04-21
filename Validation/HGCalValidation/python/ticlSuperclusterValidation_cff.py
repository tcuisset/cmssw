# Validation for the superclustering sequence in TICL
# monitors efficiency of the PID cut and of supercluster building (in trackster dataformat)
# to be used on samples containing prompt electrons (or photons)
import FWCore.ParameterSet.Config as cms

from Validation.HGCalValidation.tracksterAssociationMaskProducer_cfi import tracksterAssociationMaskProducer as _tracksterAssociationMaskProducer
from Validation.HGCalValidation.tracksterSuperclusteringValidCandidateMaskProducer_cfi import tracksterSuperclusteringValidCandidateMaskProducer as _tracksterSuperclusteringValidCandidateMaskProducer


sourceTracksterIteration = "ticlTrackstersCLUE3DHigh" # trackster collection used as input to superclustering (CLUE3D tracksters)
superclusterTracksterIteration = "ticlTracksterLinksSuperclusteringDNN" # trackster collection output by superclustering (superclusters in trackster dataformat)
simTrackstersCollection = "ticlSimTrackstersfromCPs" # simtrackster collection used for genmatching (CaloParticle)
sim2recoscore_forEfficiency = 0.3 # cut on sim2reco score to consider a supercluster genmatched for efficiency plots
sim2recoscore_forResolution = 0.99 # cut for resolution plots (looser as otherwise the scale is biased to always have Ereco>=Etrue)

################ Validation of PID cut for superclustering

## Building masks to select tracksters for efficiency and fake rate "baseline"
# finding the reco trackster that is the "seed" cluster for each electron caloparticle
ticlValidSuperclusteringSeedMask = _tracksterAssociationMaskProducer.clone(
    tracksters=cms.InputTag(sourceTracksterIteration),
    associatorRecoToSim=cms.InputTag(f"allTrackstersToSimTrackstersAssociationsByLCs:{sourceTracksterIteration}To{simTrackstersCollection}"),
    associatorSimToReco=cms.InputTag(f"allTrackstersToSimTrackstersAssociationsByLCs:{simTrackstersCollection}To{sourceTracksterIteration}"),
    recoToSimScoreCut=cms.double(0.1), # for the seed cluster of an electron, we want a tight cut (will remove low energy electron clusters that are heavily pileup contaminated)
    recoToSimScoreCut_forFakes=cms.double(0.8), # we don't want to study fake rates of PID in tracksters that have some electron/photon in it
    particleTypesSignal=cms.vint32(0, 1) # ele+photon (SimTrackster ParticleType enum, photon=0, ele=1, muon=2, pi0=3, h+-=4, h0=5)
)

# finding all reco tracksters that are a "brem"/"superclustering candidate" cluster of an electron
tracksterSuperclusteringValidCandidateMaskProducer = _tracksterSuperclusteringValidCandidateMaskProducer.clone(
    tracksters=cms.InputTag(sourceTracksterIteration),
    associatorRecoToSim=cms.InputTag(f"allTrackstersToSimTrackstersAssociationsByLCs:{sourceTracksterIteration}To{simTrackstersCollection}"),
    associatorSimToReco=cms.InputTag(f"allTrackstersToSimTrackstersAssociationsByLCs:{simTrackstersCollection}To{sourceTracksterIteration}"),
    recoToSimScoreCut=cms.double(0.15), # for the candidate brem clusters of an electron, we want a slightly looser cut
    particleTypesSignal=cms.vint32(0, 1) # ele+photon (SimTrackster ParticleType enum, photon=0, ele=1, muon=2, pi0=3, h+-=4, h0=5)
)

## actual validation DQM
from Validation.HGCalValidation.ticlTracksterPIDValidation_cfi import ticlTracksterPIDValidation as _ticlTracksterPIDValidation
ticlValidSuperclusteringSeedPID = _ticlTracksterPIDValidation.clone(
    tracksters = cms.InputTag(sourceTracksterIteration),
    tracksterMask = cms.InputTag("ticlValidSuperclusteringSeedMask"),
    tracksterMaskFakes = cms.InputTag("ticlValidSuperclusteringSeedMask", "fakes"),
    folder = cms.string("HGCAL/TICLTracksterPIDValidation/superclusteringSeedTrackster/"),
    pidCut = cms.double(0.5)
)
ticlValidSuperclusteringCandidatePID = _ticlTracksterPIDValidation.clone(
    tracksters = cms.InputTag(sourceTracksterIteration),
    tracksterMask = cms.InputTag("tracksterSuperclusteringValidCandidateMaskProducer"),
    # for fakes the selection is the same as for seeds, only the PID cut is potentially different
    tracksterMaskFakes = cms.InputTag("ticlValidSuperclusteringSeedMask", "fakes"),
    folder = cms.string("HGCAL/TICLTracksterPIDValidation/superclusteringCandidateTrackster/"),
    pidCut = cms.double(0.2),
)

ticlSuperclusterPIDValidation = cms.Sequence(
    ticlValidSuperclusteringSeedMask + tracksterSuperclusteringValidCandidateMaskProducer +
    ticlValidSuperclusteringSeedPID + ticlValidSuperclusteringCandidatePID
)



#### post-processing : computing efficiencies of PID cut as ratios of histogram
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester
postProcessorTICLPIDValid = DQMEDHarvester('DQMGenericClient',
    subDirs = cms.untracked.vstring("HGCAL/TICLTracksterPIDValidation/superclusteringSeedTrackster", "HGCAL/TICLTracksterPIDValidation/superclusteringCandidateTrackster"),
    efficiencySets = cms.untracked.VPSet(
        cms.untracked.PSet( # validating the efficiency of the gen-matching selections on electron caloparticles
            name=cms.untracked.string("pt_eta_reco2SimSelection_eff"),
            title=cms.untracked.string("Efficiency of reco2Sim cut as a function of pt and abs(eta) on superclusteringSeedTrackster"),
            numerator=cms.untracked.string("pt_eta_reco2SimSelected"),
            denominator=cms.untracked.string("pt_eta_noReco2SimSelection")
        ),

        cms.untracked.PSet(
            name=cms.untracked.string("pt_eta_PID_eff"),
            title=cms.untracked.string("Efficiency of PID cut as a function of pt and abs(eta) on superclusteringSeedTrackster"),
            numerator=cms.untracked.string("pt_eta_pidNum"),
            denominator=cms.untracked.string("pt_eta_reco2SimSelected")
        ),
        cms.untracked.PSet(
            name=cms.untracked.string("pt_eta_PID_fakeRate"),
            title=cms.untracked.string("Fake rate of PID cut as a function of pt and abs(eta) on superclusteringSeedTrackster (computed on tracksters failing sim-matching to electron)"),
            numerator=cms.untracked.string("pt_eta_fakes_pid_Num"),
            denominator=cms.untracked.string("pt_eta_fakes")
        ),

        
    ),
    efficiency = cms.vstring(),
    resolution = cms.vstring(),
    verbose = cms.untracked.uint32(4))


################ Validation of superclusters (in trackster dataformat)
from Validation.HGCalValidation.hgcalSuperClusterValidator_cfi import hgcalSuperClusterValidator as _hgcalSuperClusterValidator
hgcalSuperClusterValidator = _hgcalSuperClusterValidator.clone(
    sc_tracksters = cms.InputTag(superclusterTracksterIteration),
    sc_tracksters_before_linking = cms.InputTag(sourceTracksterIteration),
    linkedTracksterIdToInputTracksterId = cms.InputTag(superclusterTracksterIteration, "linkedTracksterIdToInputTracksterId"),
    associatorsimToReco = cms.InputTag(f"allTrackstersToSimTrackstersAssociationsByLCs:{simTrackstersCollection}To{superclusterTracksterIteration}"),
    simToRecoScoreThreshold_forEfficiency = sim2recoscore_forEfficiency,
    simToRecoScoreThreshold_forResolution = sim2recoscore_forResolution
)

from Configuration.ProcessModifiers.ticl_superclustering_mustache_ticl_cff import ticl_superclustering_mustache_ticl
ticl_superclustering_mustache_ticl.toModify(hgcalSuperClusterValidator,
    sc_tracksters = cms.InputTag("ticlTracksterLinksSuperclusteringMustache"),
    linkedTracksterIdToInputTracksterId = cms.InputTag("ticlTracksterLinksSuperclusteringMustache", "linkedTracksterIdToInputTracksterId"))


## post-processing 
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester
postProcessorHGCalSuperClusterValidator = DQMEDHarvester('DQMGenericClient',
    subDirs = cms.untracked.vstring("HGCal/SuperClusters"),
    efficiencySets = cms.untracked.VPSet(
        cms.untracked.PSet(
            name=cms.untracked.string("sc_ET_eff"),
            title=cms.untracked.string(f"Efficiency of SuperCluster reconstruction (denom=e/g CaloParticle, num=denom+supercluster with sim2reco score<{sim2recoscore_forEfficiency})"),
            numerator=cms.untracked.string("sim_ET_vs_Eta_num"),
            denominator=cms.untracked.string("sim_ET_vs_Eta_denom")
        ),
    ),
    resolutionProfileSets = cms.untracked.VPSet( # resolutionProfileSets
        cms.untracked.PSet(
            namePrefix=cms.untracked.string("sc_resolution_vs_pt"),
            titlePrefix=cms.untracked.string("HGCAL SuperCluster energy resolution (Esc / Etrue) vs ET (RMS, Etrue=simTrackster->regressed_energy)"),
            srcName=cms.untracked.string("sc_EoverEtruth_vs_ET"),
            typeName=cms.untracked.string("rms"),
        ),

        cms.untracked.PSet(
            namePrefix=cms.untracked.string("sc_resolution_vs_pt_lowEta"),
            titlePrefix=cms.untracked.string("HGCAL SuperCluster energy resolution (Esc / Etrue) vs ET (abs(eta)<2.1 region) (RMS, Etrue=simTrackster->regressed_energy)"),
            srcName=cms.untracked.string("sc_EoverEtruth_vs_ET_lowEta"),
            typeName=cms.untracked.string("rms"),
        ),
        cms.untracked.PSet(
            namePrefix=cms.untracked.string("sc_resolution_vs_pt_midEta"),
            titlePrefix=cms.untracked.string("HGCAL SuperCluster energy resolution (Esc / Etrue) vs ET (2.1<abs(eta)<2.6 region) (RMS, Etrue=simTrackster->regressed_energy)"),
            srcName=cms.untracked.string("sc_EoverEtruth_vs_ET_midEta"),
            typeName=cms.untracked.string("rms"),
        ),
        cms.untracked.PSet(
            namePrefix=cms.untracked.string("sc_resolution_vs_pt_highEta"),
            titlePrefix=cms.untracked.string("HGCAL SuperCluster energy resolution (Esc / Etrue) vs ET (abs(eta)>2.6 region) (RMS, Etrue=simTrackster->regressed_energy)"),
            srcName=cms.untracked.string("sc_EoverEtruth_vs_ET_highEta"),
            typeName=cms.untracked.string("rms"),
        ),

    ),
    efficiency = cms.vstring(),
    resolution = cms.vstring(),
    verbose = cms.untracked.uint32(4))


########################### Sequences
ticlSuperclusterValidation = cms.Sequence(
    ticlSuperclusterPIDValidation + hgcalSuperClusterValidator
)
postProcessorTiclSupercluster = cms.Sequence(
    postProcessorTICLPIDValid + postProcessorHGCalSuperClusterValidator
)
