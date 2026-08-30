// Flat ntuple for TICL supercluster energy-response studies.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "TTree.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"

namespace {
  using TracksterCollection = std::vector<ticl::Trackster>;
  using TracksterToTracksterMap =
      ticl::AssociationMap<ticl::mapWithSharedEnergyAndScore, TracksterCollection, TracksterCollection>;

  constexpr float kMissingFloat = -999.f;
}  // namespace

class SuperclusterNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SuperclusterNtuplizer(edm::ParameterSet const&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void beginJob() override;
  void analyze(edm::Event const&, edm::EventSetup const&) override;
  void resetRow();

  edm::EDGetTokenT<TracksterCollection> superclustersToken_;
  edm::EDGetTokenT<std::vector<std::vector<unsigned int>>> linksToken_;
  edm::EDGetTokenT<reco::SuperClusterCollection> recoSuperclustersToken_;
  edm::EDGetTokenT<TracksterCollection> simTrackstersToken_;
  edm::EDGetTokenT<TracksterToTracksterMap> recoToSimToken_;
  edm::EDGetTokenT<TracksterToTracksterMap> simToRecoToken_;
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  float recoSuperclusterEtThreshold_;

  TTree* tree_ = nullptr;

  uint32_t run_ = 0;
  uint32_t lumi_ = 0;
  uint64_t event_ = 0;
  uint32_t scIndex_ = 0;
  int32_t recoScIndex_ = -1;
  float scRawEnergy_ = kMissingFloat;
  float scRegressedEnergy_ = kMissingFloat;
  float scRawPt_ = kMissingFloat;
  float scEta_ = kMissingFloat;
  float scPhi_ = kMissingFloat;
  uint32_t scNTracksters_ = 0;
  int32_t scSeedTracksterIndex_ = -1;
  bool hasRecoSupercluster_ = false;
  float recoScRawEnergy_ = kMissingFloat;
  float recoScEnergy_ = kMissingFloat;
  float recoScCorrectedEnergy_ = kMissingFloat;

  uint32_t nRecoToSimAssociations_ = 0;
  bool hasBestSimTrackster_ = false;
  int32_t bestSimTracksterIndex_ = -1;
  float bestSimSharedEnergy_ = kMissingFloat;
  float bestSimRecoToSimScore_ = kMissingFloat;
  float bestSimRawEnergy_ = kMissingFloat;
  float bestSimRegressedEnergy_ = kMissingFloat;
  float bestSimEta_ = kMissingFloat;
  float bestSimPhi_ = kMissingFloat;
  bool bestSimIsHadronic_ = false;
  int32_t bestSimCaloParticleIndex_ = -1;
  int32_t bestSimPdgId_ = 0;
  float bestSimCaloParticleEnergy_ = kMissingFloat;
  float bestSimCaloParticlePt_ = kMissingFloat;
  float bestSimCaloParticleEta_ = kMissingFloat;
  float bestSimCaloParticlePhi_ = kMissingFloat;

  bool isBestRecoForSim_ = false;
  int32_t bestRecoIndexForSim_ = -1;
  float bestRecoForSimSharedEnergy_ = kMissingFloat;
  float bestRecoForSimScore_ = kMissingFloat;
};

SuperclusterNtuplizer::SuperclusterNtuplizer(edm::ParameterSet const& config)
    : superclustersToken_(consumes<TracksterCollection>(config.getParameter<edm::InputTag>("superclusters"))),
      linksToken_(
          consumes<std::vector<std::vector<unsigned int>>>(config.getParameter<edm::InputTag>("linkedTracksters"))),
      recoSuperclustersToken_(
          consumes<reco::SuperClusterCollection>(config.getParameter<edm::InputTag>("recoSuperclusters"))),
      simTrackstersToken_(consumes<TracksterCollection>(config.getParameter<edm::InputTag>("simTracksters"))),
      recoToSimToken_(consumes<TracksterToTracksterMap>(config.getParameter<edm::InputTag>("recoToSimAssociator"))),
      simToRecoToken_(consumes<TracksterToTracksterMap>(config.getParameter<edm::InputTag>("simToRecoAssociator"))),
      caloParticlesToken_(consumes<std::vector<CaloParticle>>(config.getParameter<edm::InputTag>("caloParticles"))),
      recoSuperclusterEtThreshold_(config.getParameter<double>("recoSuperclusterEtThreshold")) {
  usesResource("TFileService");
}

void SuperclusterNtuplizer::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("superclusters", "One row per TICL supercluster");

  tree_->Branch("run", &run_);
  tree_->Branch("lumi", &lumi_);
  tree_->Branch("event", &event_);
  tree_->Branch("sc_index", &scIndex_);
  tree_->Branch("reco_sc_index", &recoScIndex_);
  tree_->Branch("sc_raw_energy", &scRawEnergy_);
  tree_->Branch("sc_regressed_energy", &scRegressedEnergy_);
  tree_->Branch("sc_raw_pt", &scRawPt_);
  tree_->Branch("sc_eta", &scEta_);
  tree_->Branch("sc_phi", &scPhi_);
  tree_->Branch("sc_num_tracksters", &scNTracksters_);
  tree_->Branch("sc_seed_trackster_index", &scSeedTracksterIndex_);
  tree_->Branch("has_reco_supercluster", &hasRecoSupercluster_);
  tree_->Branch("reco_sc_raw_energy", &recoScRawEnergy_);
  tree_->Branch("reco_sc_energy", &recoScEnergy_);
  tree_->Branch("reco_sc_corrected_energy", &recoScCorrectedEnergy_);

  tree_->Branch("num_reco_to_sim_associations", &nRecoToSimAssociations_);
  tree_->Branch("has_best_simtrackster", &hasBestSimTrackster_);
  tree_->Branch("best_simtrackster_index", &bestSimTracksterIndex_);
  tree_->Branch("best_sim_shared_energy", &bestSimSharedEnergy_);
  tree_->Branch("best_sim_reco_to_sim_score", &bestSimRecoToSimScore_);
  tree_->Branch("best_sim_raw_energy", &bestSimRawEnergy_);
  tree_->Branch("best_sim_regressed_energy", &bestSimRegressedEnergy_);
  tree_->Branch("best_sim_eta", &bestSimEta_);
  tree_->Branch("best_sim_phi", &bestSimPhi_);
  tree_->Branch("best_sim_is_hadronic", &bestSimIsHadronic_);
  tree_->Branch("best_sim_caloparticle_index", &bestSimCaloParticleIndex_);
  tree_->Branch("best_sim_pdgid", &bestSimPdgId_);
  tree_->Branch("best_sim_caloparticle_energy", &bestSimCaloParticleEnergy_);
  tree_->Branch("best_sim_caloparticle_pt", &bestSimCaloParticlePt_);
  tree_->Branch("best_sim_caloparticle_eta", &bestSimCaloParticleEta_);
  tree_->Branch("best_sim_caloparticle_phi", &bestSimCaloParticlePhi_);

  tree_->Branch("is_best_reco_for_sim", &isBestRecoForSim_);
  tree_->Branch("best_reco_index_for_sim", &bestRecoIndexForSim_);
  tree_->Branch("best_reco_for_sim_shared_energy", &bestRecoForSimSharedEnergy_);
  tree_->Branch("best_reco_for_sim_score", &bestRecoForSimScore_);
}

void SuperclusterNtuplizer::resetRow() {
  recoScIndex_ = -1;
  scRawEnergy_ = kMissingFloat;
  scRegressedEnergy_ = kMissingFloat;
  scRawPt_ = kMissingFloat;
  scEta_ = kMissingFloat;
  scPhi_ = kMissingFloat;
  scNTracksters_ = 0;
  scSeedTracksterIndex_ = -1;
  hasRecoSupercluster_ = false;
  recoScRawEnergy_ = kMissingFloat;
  recoScEnergy_ = kMissingFloat;
  recoScCorrectedEnergy_ = kMissingFloat;

  nRecoToSimAssociations_ = 0;
  hasBestSimTrackster_ = false;
  bestSimTracksterIndex_ = -1;
  bestSimSharedEnergy_ = kMissingFloat;
  bestSimRecoToSimScore_ = kMissingFloat;
  bestSimRawEnergy_ = kMissingFloat;
  bestSimRegressedEnergy_ = kMissingFloat;
  bestSimEta_ = kMissingFloat;
  bestSimPhi_ = kMissingFloat;
  bestSimIsHadronic_ = false;
  bestSimCaloParticleIndex_ = -1;
  bestSimPdgId_ = 0;
  bestSimCaloParticleEnergy_ = kMissingFloat;
  bestSimCaloParticlePt_ = kMissingFloat;
  bestSimCaloParticleEta_ = kMissingFloat;
  bestSimCaloParticlePhi_ = kMissingFloat;

  isBestRecoForSim_ = false;
  bestRecoIndexForSim_ = -1;
  bestRecoForSimSharedEnergy_ = kMissingFloat;
  bestRecoForSimScore_ = kMissingFloat;
}

void SuperclusterNtuplizer::analyze(edm::Event const& event, edm::EventSetup const&) {
  auto const superclustersHandle = event.getHandle(superclustersToken_);
  auto const linksHandle = event.getHandle(linksToken_);
  auto const recoSuperclustersHandle = event.getHandle(recoSuperclustersToken_);
  auto const simTrackstersHandle = event.getHandle(simTrackstersToken_);
  auto const caloParticlesHandle = event.getHandle(caloParticlesToken_);
  auto const& recoToSim = event.get(recoToSimToken_);
  auto const& simToReco = event.get(simToRecoToken_);

  auto const& superclusters = *superclustersHandle;
  auto const& links = *linksHandle;
  auto const& recoSuperclusters = *recoSuperclustersHandle;
  auto const& simTracksters = *simTrackstersHandle;
  auto const& caloParticles = *caloParticlesHandle;

  if (superclusters.size() != links.size() || recoToSim.size() != superclusters.size() ||
      simToReco.size() != simTracksters.size()) {
    throw cms::Exception("DataCorruption")
        << "SuperclusterNtuplizer found inconsistent collection/map sizes: superclusters=" << superclusters.size()
        << ", links=" << links.size() << ", recoToSim=" << recoToSim.size()
        << ", simTracksters=" << simTracksters.size() << ", simToReco=" << simToReco.size();
  }
  if (recoToSim.getCollectionIDs().first.id() != superclustersHandle.id() ||
      recoToSim.getCollectionIDs().second.id() != simTrackstersHandle.id() ||
      simToReco.getCollectionIDs().first.id() != simTrackstersHandle.id() ||
      simToReco.getCollectionIDs().second.id() != superclustersHandle.id()) {
    throw cms::Exception("ProductMismatch")
        << "The configured association maps do not reference the configured reco and SimTrackster collections.";
  }

  run_ = event.id().run();
  lumi_ = event.id().luminosityBlock();
  event_ = event.id().event();

  std::size_t nextRecoScIndex = 0;
  for (std::size_t index = 0; index < superclusters.size(); ++index) {
    resetRow();
    scIndex_ = index;

    auto const& sc = superclusters[index];
    scRawEnergy_ = sc.raw_energy();
    scRegressedEnergy_ = sc.regressed_energy();
    scRawPt_ = sc.raw_pt();
    scEta_ = sc.barycenter().eta();
    scPhi_ = sc.barycenter().phi();
    scNTracksters_ = links[index].size();
    if (!links[index].empty())
      scSeedTracksterIndex_ = links[index].front();

    if (sc.raw_pt() >= recoSuperclusterEtThreshold_) {
      if (nextRecoScIndex >= recoSuperclusters.size()) {
        throw cms::Exception("ProductMismatch")
            << "Fewer reco::SuperClusters than expected from the configured ET threshold "
            << recoSuperclusterEtThreshold_;
      }
      auto const& recoSc = recoSuperclusters[nextRecoScIndex];
      const float energyScale = std::max(1.f, std::abs(sc.raw_energy()));
      if (std::abs(recoSc.rawEnergy() - sc.raw_energy()) > 1.e-5f * energyScale ||
          std::abs(recoSc.eta() - sc.barycenter().eta()) > 1.e-5f) {
        throw cms::Exception("ProductMismatch")
            << "reco::SuperCluster index " << nextRecoScIndex << " does not match TICL supercluster index " << index
            << ". Check recoSuperclusterEtThreshold against EGammaSuperclusterProducer.";
      }
      hasRecoSupercluster_ = true;
      recoScIndex_ = nextRecoScIndex;
      recoScRawEnergy_ = recoSc.rawEnergy();
      recoScEnergy_ = recoSc.energy();
      recoScCorrectedEnergy_ = recoSc.correctedEnergy();
      ++nextRecoScIndex;
    }

    auto const& associations = recoToSim[index];
    nRecoToSimAssociations_ = associations.size();
    auto const bestSim = std::max_element(
        associations.begin(),
        associations.end(),
        [](TracksterToTracksterMap::AssociationElementType const& a,
           TracksterToTracksterMap::AssociationElementType const& b) { return a.sharedEnergy() < b.sharedEnergy(); });

    if (bestSim != associations.end()) {
      hasBestSimTrackster_ = true;
      bestSimTracksterIndex_ = bestSim->index();
      bestSimSharedEnergy_ = bestSim->sharedEnergy();
      bestSimRecoToSimScore_ = bestSim->score();

      auto const& sim = simTracksters[bestSim->index()];
      bestSimRawEnergy_ = sim.raw_energy();
      bestSimRegressedEnergy_ = sim.regressed_energy();
      bestSimEta_ = sim.barycenter().eta();
      bestSimPhi_ = sim.barycenter().phi();
      bestSimIsHadronic_ = sim.isHadronic();

      if (sim.seedID() == caloParticlesHandle.id() && sim.seedIndex() >= 0 &&
          static_cast<std::size_t>(sim.seedIndex()) < caloParticles.size()) {
        bestSimCaloParticleIndex_ = sim.seedIndex();
        auto const& cp = caloParticles[sim.seedIndex()];
        bestSimPdgId_ = cp.pdgId();
        bestSimCaloParticleEnergy_ = cp.energy();
        bestSimCaloParticlePt_ = cp.pt();
        bestSimCaloParticleEta_ = cp.eta();
        bestSimCaloParticlePhi_ = cp.phi();
      }

      auto const& reverseAssociations = simToReco[bestSim->index()];
      auto const bestReco = std::max_element(
          reverseAssociations.begin(),
          reverseAssociations.end(),
          [](TracksterToTracksterMap::AssociationElementType const& a,
             TracksterToTracksterMap::AssociationElementType const& b) { return a.sharedEnergy() < b.sharedEnergy(); });
      if (bestReco != reverseAssociations.end()) {
        bestRecoIndexForSim_ = bestReco->index();
        bestRecoForSimSharedEnergy_ = bestReco->sharedEnergy();
        bestRecoForSimScore_ = bestReco->score();
        isBestRecoForSim_ = bestReco->index() == index;
      }
    }

    tree_->Fill();
  }

  if (nextRecoScIndex != recoSuperclusters.size()) {
    throw cms::Exception("ProductMismatch")
        << "Found " << recoSuperclusters.size() << " reco::SuperClusters but matched " << nextRecoScIndex
        << " using recoSuperclusterEtThreshold=" << recoSuperclusterEtThreshold_;
  }
}

void SuperclusterNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription description;
  description.add<edm::InputTag>("superclusters", edm::InputTag("ticlTracksterLinksSuperclusteringDNN"));
  description.add<edm::InputTag>(
      "linkedTracksters", edm::InputTag("ticlTracksterLinksSuperclusteringDNN", "linkedTracksterIdToInputTracksterId"));
  description.add<edm::InputTag>("recoSuperclusters", edm::InputTag("ticlEGammaSuperClusterProducer"));
  description.add<edm::InputTag>("simTracksters", edm::InputTag("ticlSimTracksters", "fromCPs"));
  description.add<edm::InputTag>("recoToSimAssociator",
                                 edm::InputTag("allTrackstersToSimTrackstersAssociationsByLCs",
                                               "ticlTracksterLinksSuperclusteringDNNToticlSimTrackstersfromCPs"));
  description.add<edm::InputTag>("simToRecoAssociator",
                                 edm::InputTag("allTrackstersToSimTrackstersAssociationsByLCs",
                                               "ticlSimTrackstersfromCPsToticlTracksterLinksSuperclusteringDNN"));
  description.add<edm::InputTag>("caloParticles", edm::InputTag("mix", "MergedCaloTruth"));
  description.add<double>("recoSuperclusterEtThreshold", 4.)
      ->setComment("Must match ticlEGammaSuperClusterProducer.superclusterEtThreshold.");
  descriptions.add("superclusterNtuplizer", description);
}

DEFINE_FWK_MODULE(SuperclusterNtuplizer);
