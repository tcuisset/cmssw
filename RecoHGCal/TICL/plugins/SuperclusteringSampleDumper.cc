// Original Author:  Theo Cuisset
//         Created:  Nov 2023
/** Produce samples for electron superclustering DNN training in TICL

 Pairs of seed-candidate tracksters (in compatible geometric windows) are iterated over, in similar manner as in TracksterLinkingBySuperclustering.
 For each of these pairs, the DNN features are computed and saved to a TTree. 
 Also saved is the best (=lowest) association score of the seed trackster with CaloParticles. The association score of the candidate trackster 
 with the same CaloParticle is also saved.
*/
#include <cmath>
#include <vector>
#include <memory>
#include <algorithm>
#include <numeric>

#include <TTree.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/allowedValues.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"

#include "RecoHGCal/TICL/plugins/TracksterLinkingbySuperClusteringDNN.h"
#include "RecoHGCal/TICL/interface/SuperclusteringDNNInputs.h"

using namespace ticl;

class SuperclusteringSampleDumper : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SuperclusteringSampleDumper(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  float explainedVarianceRatio(ticl::Trackster const& ts) const;
  bool checkExplainedVarianceRatioCut(ticl::Trackster const& ts) const;

  const edm::EDGetTokenT<std::vector<Trackster>> tracksters_clue3d_token_;
  const edm::EDGetTokenT<ticl::AssociationMap<ticl::mapWithSharedEnergyAndScore, std::vector<ticl::Trackster>, std::vector<ticl::Trackster>>> tsRecoToSimCP_token_;
  float deltaEtaWindow_;
  float deltaPhiWindow_;
  float seedPtThreshold_;
  float candidateEnergyThreshold_;
  bool enableExplainedVarianceRatioCut_;
  float explVarRatioCut_energyBoundary_;  // Boundary energy between low and high energy explVarRatio cut threshold
  float explVarRatioMinimum_lowEnergy_;  // Cut on explained variance ratio of tracksters to be considered as candidate, for trackster raw_energy < explVarRatioCut_energyBoundary
  float explVarRatioMinimum_highEnergy_;  // Cut on explained variance ratio of tracksters to be considered as candidate, for trackster raw_energy > explVarRatioCut_energyBoundary

  TTree* output_tree_;
  edm::EventID eventId_;
  std::unique_ptr<AbstractSuperclusteringDNNInput> dnnInput_;
  std::vector<std::vector<float>>
      features_;  // Outer index : feature number (split into branches), inner index : inference pair index
  std::vector<unsigned int> seedTracksterIdx_;       // ID of seed trackster used for inference pair
  std::vector<unsigned int> candidateTracksterIdx_;  // ID of candidate trackster used for inference pair
  std::vector<float> seedTracksterExplainedVarianceRatio_;
  std::vector<float> candidateTracksterExplainedVarianceRatio_;
  std::vector<float> seedTracksterPIDScore_;
  std::vector<float> candidateTracksterPIDScore_;

  std::vector<float>
      seedTracksterBestAssociationScore_;  // Best association score of seed trackster (seedTracksterIdx) with CaloParticle
  std::vector<long>
      seedTracksterBestAssociation_simTsIdx_;  // Index of SimTrackster that has the best association score to the seedTrackster
  std::vector<float> seedTracksterBestAssociation_caloParticleEnergy_;  // Energy of best associated CaloParticle to seed

  std::vector<float>
      candidateTracksterBestAssociationScore_;  // Best association score of candidate trackster (seedTracksterIdx) with CaloParticle
  std::vector<long>
      candidateTracksterBestAssociation_simTsIdx_;  // Index of SimTrackster that has the best association score to the candidate

  std::vector<float>
      candidateTracksterAssociationWithSeed_score_;  // Association score of candidate trackster with the CaloParticle of seedTracksterBestAssociation_simTsIdx_
};

SuperclusteringSampleDumper::SuperclusteringSampleDumper(const edm::ParameterSet& ps)
    : tracksters_clue3d_token_(consumes<std::vector<Trackster>>(ps.getParameter<edm::InputTag>("tracksters"))),
      tsRecoToSimCP_token_(
          consumes<ticl::AssociationMap<ticl::mapWithSharedEnergyAndScore, std::vector<ticl::Trackster>, std::vector<ticl::Trackster>>>(ps.getParameter<edm::InputTag>("recoToSimAssociatorCP"))),
      deltaEtaWindow_(ps.getParameter<double>("deltaEtaWindow")),
      deltaPhiWindow_(ps.getParameter<double>("deltaPhiWindow")),
      seedPtThreshold_(ps.getParameter<double>("seedPtThreshold")),
      candidateEnergyThreshold_(ps.getParameter<double>("candidateEnergyThreshold")),
      enableExplainedVarianceRatioCut_(ps.getParameter<bool>("enableExplainedVarianceRatioCut")),
      explVarRatioCut_energyBoundary_(ps.getParameter<double>("candidateEnergyThreshold")),
      explVarRatioMinimum_lowEnergy_(ps.getParameter<double>("explVarRatioMinimum_lowEnergy")),
      explVarRatioMinimum_highEnergy_(ps.getParameter<double>("explVarRatioMinimum_highEnergy")),
      eventId_(),
      dnnInput_(makeSuperclusteringDNNInputFromString(ps.getParameter<std::string>("dnnInputsVersion"))),
      features_(dnnInput_->featureCount()) {
  usesResource("TFileService");
}

void SuperclusteringSampleDumper::beginJob() {
  edm::Service<TFileService> fs;
  output_tree_ = fs->make<TTree>("superclusteringTraining", "Superclustering training samples");
  output_tree_->Branch("Event", &eventId_);
  output_tree_->Branch("seedTracksterIdx", &seedTracksterIdx_);
  output_tree_->Branch("candidateTracksterIdx", &candidateTracksterIdx_);
  output_tree_->Branch("seedTracksterExplainedVarianceRatio", &seedTracksterExplainedVarianceRatio_);
  output_tree_->Branch("candidateTracksterExplainedVarianceRatio", &candidateTracksterExplainedVarianceRatio_);
  output_tree_->Branch("seedTracksterPIDScore", &seedTracksterPIDScore_);
  output_tree_->Branch("candidateTracksterPIDScore", &candidateTracksterPIDScore_);
  output_tree_->Branch("seedTracksterBestAssociationScore", &seedTracksterBestAssociationScore_);
  output_tree_->Branch("seedTracksterBestAssociation_simTsIdx", &seedTracksterBestAssociation_simTsIdx_);
  output_tree_->Branch("seedTracksterBestAssociation_caloParticleEnergy",
                       &seedTracksterBestAssociation_caloParticleEnergy_);
  output_tree_->Branch("candidateTracksterBestAssociationScore", &candidateTracksterBestAssociationScore_);
  output_tree_->Branch("candidateTracksterBestAssociation_simTsIdx", &candidateTracksterBestAssociation_simTsIdx_);
  output_tree_->Branch("candidateTracksterAssociationWithSeed_score", &candidateTracksterAssociationWithSeed_score_);
  std::vector<std::string> featureNames = dnnInput_->featureNames();
  assert(featureNames.size() == dnnInput_->featureCount());
  for (unsigned int i = 0; i < dnnInput_->featureCount(); i++) {
    output_tree_->Branch(("feature_" + featureNames[i]).c_str(), &features_[i]);
  }
}

/** 
 * Check if trackster passes cut on explained variance ratio. The DNN is trained only on pairs where both seed and candidate pass this cut
 * Explained variance ratio is (largest PCA eigenvalue) / (sum of PCA eigenvalues)
*/
float SuperclusteringSampleDumper::explainedVarianceRatio(ticl::Trackster const& ts) const {
  float explVar_denominator =
      std::accumulate(std::begin(ts.eigenvalues()), std::end(ts.eigenvalues()), 0.f, std::plus<float>());
  if (explVar_denominator != 0.) {
    return ts.eigenvalues()[0] / explVar_denominator;
  }
  return -1.f;
}

bool SuperclusteringSampleDumper::checkExplainedVarianceRatioCut(ticl::Trackster const& ts) const {
  if (!enableExplainedVarianceRatioCut_) {
    return true;
  }

  const float explVarRatio = explainedVarianceRatio(ts);
  if (ts.raw_energy() > explVarRatioCut_energyBoundary_)
    return explVarRatio >= explVarRatioMinimum_highEnergy_;
  else
    return explVarRatio >= explVarRatioMinimum_lowEnergy_;
}

void SuperclusteringSampleDumper::analyze(const edm::Event& evt, const edm::EventSetup& iSetup) {
  eventId_ = evt.id();

  edm::Handle<std::vector<Trackster>> inputTracksters;
  evt.getByToken(tracksters_clue3d_token_, inputTracksters);

  auto const& assoc_CP_recoToSim = evt.get(tsRecoToSimCP_token_);
  assert(assoc_CP_recoToSim.getCollectionIDs().first.id() == inputTracksters.id() && "Mismatch between trackster collection and recoToSim associator");

  auto const& tracksters = *inputTracksters;
  const auto nTs = static_cast<unsigned int>(tracksters.size());
  if (nTs == 0u) {
    output_tree_->Fill();
    return;
  }

  // ---- pT-sorted indices (out-of-place)
  std::vector<unsigned int> trackstersIndicesPt(nTs);
  std::iota(trackstersIndicesPt.begin(), trackstersIndicesPt.end(), 0u);
  std::stable_sort(
      trackstersIndicesPt.begin(), trackstersIndicesPt.end(), [&tracksters](unsigned int a, unsigned int b) {
        return tracksters[a].raw_pt() > tracksters[b].raw_pt();
      });

  const auto nFeatures = static_cast<unsigned int>(features_.size());
  // Sanity: this module expects the configured feature set to match the DNN input helper.
  // (If you prefer, replace with a cms::Exception instead of assert.)
  assert(nFeatures == dnnInput_->featureCount());

  // ---- optional: reserve to reduce reallocs (cheap heuristic)
  // Worst-case is O(N^2); don't do that. Use something modest.
  // You can tune: e.g. nTs*8 tends to be safe-ish without huge memory.
  const size_t reservePairs = static_cast<size_t>(nTs) * 8u;
  for (auto& col : features_) {
    col.reserve(reservePairs);
  }
  seedTracksterIdx_.reserve(reservePairs);
  candidateTracksterIdx_.reserve(reservePairs);
  seedTracksterExplainedVarianceRatio_.reserve(reservePairs);
  candidateTracksterExplainedVarianceRatio_.reserve(reservePairs);
  seedTracksterPIDScore_.reserve(reservePairs);
  candidateTracksterPIDScore_.reserve(reservePairs);
  seedTracksterBestAssociationScore_.reserve(reservePairs);
  seedTracksterBestAssociation_simTsIdx_.reserve(reservePairs);
  seedTracksterBestAssociation_caloParticleEnergy_.reserve(reservePairs);
  candidateTracksterBestAssociationScore_.reserve(reservePairs);
  candidateTracksterBestAssociation_simTsIdx_.reserve(reservePairs);
  candidateTracksterAssociationWithSeed_score_.reserve(reservePairs);

  // ---- scratch buffer for features (no per-pair allocation)
  std::vector<float> featScratch;
  featScratch.resize(nFeatures);

  // Find the best association score among all associations to one trackster
  auto bestAssoc = [](std::vector<AssociationElement<std::pair<SharedEnergyType, float>>> const& val) {
    return *std::min_element(val.begin(), val.end(), [](auto const& a, auto const& b) {
      // pair<Ref, pair<sharedEnergy, associationScore>>; best is smallest score
      return a.score() < b.score();
    });
  };

  // Outer: seed
  for (unsigned int seed_pt = 0; seed_pt < nTs; ++seed_pt) {
    const unsigned int seed_idx = trackstersIndicesPt[seed_pt];
    auto const& ts_seed = tracksters[seed_idx];

    if (ts_seed.raw_pt() < seedPtThreshold_) {
      break;  // remaining seeds are lower-pT due to sorting
    }
    const float seedExplainedVarianceRatio = explainedVarianceRatio(ts_seed);
    const float seedPIDScore = ts_seed.id_probability(Trackster::ParticleType::electron);
    if (!checkExplainedVarianceRatioCut(ts_seed)) {
      continue;
    }

    // Find best associated CaloParticle to the seed
    if (assoc_CP_recoToSim[seed_idx].empty()) continue;
    auto const& seed_best = bestAssoc(assoc_CP_recoToSim[seed_idx]);

    // Inner: candidate (only lower-pT than seed)
    for (unsigned int cand_pt = seed_pt + 1; cand_pt < nTs; ++cand_pt) {
      const unsigned int cand_idx = trackstersIndicesPt[cand_pt];
      auto const& ts_cand = tracksters[cand_idx];

      if (ts_cand.raw_energy() < candidateEnergyThreshold_) {
        continue;
      }
      const float candidateExplainedVarianceRatio = explainedVarianceRatio(ts_cand);
      const float candidatePIDScore = ts_cand.id_probability(Trackster::ParticleType::electron);
      if (!checkExplainedVarianceRatioCut(ts_cand)) {
        continue;
      }
      if (std::abs(ts_seed.barycenter().Eta() - ts_cand.barycenter().Eta()) >= deltaEtaWindow_) {
        continue;
      }
      if (std::abs(deltaPhi(ts_seed.barycenter().Phi(), ts_cand.barycenter().Phi())) >= deltaPhiWindow_) {
        continue;
      }

      // ---- compute features in-place
      dnnInput_->computeInto(ts_seed, ts_cand, std::span<float>(featScratch.data(), featScratch.size()));

      // ---- store features (columnar vectors for the tree)
      // Unrolled-ish simple loop, good locality on featScratch
      for (unsigned int f = 0; f < nFeatures; ++f) {
        features_[f].push_back(featScratch[f]);
      }

      seedTracksterIdx_.push_back(seed_idx);
      candidateTracksterIdx_.push_back(cand_idx);
      seedTracksterExplainedVarianceRatio_.push_back(seedExplainedVarianceRatio);
      candidateTracksterExplainedVarianceRatio_.push_back(candidateExplainedVarianceRatio);
      seedTracksterPIDScore_.push_back(seedPIDScore);
      candidateTracksterPIDScore_.push_back(candidatePIDScore);

      float candBestScore = 1.f;
      long candBestSimIdx = -1;
      float candScoreWithSeed = 1.f;

      auto& cand_assocs = assoc_CP_recoToSim[cand_idx];
      if (!cand_assocs.empty()) {
        auto const& cand_best = bestAssoc(cand_assocs);
        candBestScore = cand_best.score();
        candBestSimIdx = cand_best.index();

        auto itSeed = std::find_if(cand_assocs.begin(), cand_assocs.end(), [&seed_best](auto const& assoc) {
          return assoc.index() == seed_best.index();
        });
        if (itSeed != cand_assocs.end()) {
          candScoreWithSeed = itSeed->score();
        }
      }

      seedTracksterBestAssociationScore_.push_back(seed_best.score());
      seedTracksterBestAssociation_simTsIdx_.push_back(seed_best.index());
      seedTracksterBestAssociation_caloParticleEnergy_.push_back(assoc_CP_recoToSim.getRefSecond(seed_best.index())->regressed_energy());

      candidateTracksterBestAssociationScore_.push_back(candBestScore);
      candidateTracksterBestAssociation_simTsIdx_.push_back(candBestSimIdx);
      candidateTracksterAssociationWithSeed_score_.push_back(candScoreWithSeed);
    }
  }

  output_tree_->Fill();

  // Clear but keep capacity (important for perf across events)
  for (auto& col : features_) {
    col.clear();
  }
  seedTracksterIdx_.clear();
  candidateTracksterIdx_.clear();
  seedTracksterExplainedVarianceRatio_.clear();
  candidateTracksterExplainedVarianceRatio_.clear();
  seedTracksterPIDScore_.clear();
  candidateTracksterPIDScore_.clear();
  seedTracksterBestAssociationScore_.clear();
  seedTracksterBestAssociation_simTsIdx_.clear();
  seedTracksterBestAssociation_caloParticleEnergy_.clear();
  candidateTracksterBestAssociationScore_.clear();
  candidateTracksterBestAssociation_simTsIdx_.clear();
  candidateTracksterAssociationWithSeed_score_.clear();
}

void SuperclusteringSampleDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracksters", edm::InputTag("ticlTrackstersCLUE3DHigh"))
      ->setComment("Input trackster collection, same as what is used for superclustering inference.");
  desc.add<edm::InputTag>("recoToSimAssociatorCP",
                          edm::InputTag("allTrackstersToSimTrackstersAssociationsByLCs", "ticlTrackstersCLUE3DHighToticlSimTrackstersfromCPs"));
  desc.ifValue(edm::ParameterDescription<std::string>("dnnInputsVersion", "v3", true),
               edm::allowedValues<std::string>("v1", "v2", "v3"))
      ->setComment(
          "DNN inputs version tag. Defines which set of features is fed to the DNN. Must match with the actual DNN.");
  // Cuts are intentionally looser than those used for inference in TracksterLinkingBySuperClustering.cpp
  desc.add<double>("deltaEtaWindow", 0.2)
      ->setComment(
          "Size of delta eta window to consider for superclustering. Seed-candidate pairs outside this window "
          "are not considered for DNN inference.");
  desc.add<double>("deltaPhiWindow", 0.7)
      ->setComment(
          "Size of delta phi window to consider for superclustering. Seed-candidate pairs outside this window "
          "are not considered for DNN inference.");
  desc.add<double>("seedPtThreshold", 3.)
      ->setComment("Minimum transverse momentum of trackster to be considered as seed of a supercluster");
  desc.add<double>("candidateEnergyThreshold", 1.5)
      ->setComment("Minimum energy of trackster to be considered as candidate for superclustering");
  desc.add<bool>("enableExplainedVarianceRatioCut", true)
      ->setComment(
          "Whether to apply the explained variance ratio selection to seed and candidate tracksters. If false, "
          "tracksters are not selected on explained variance ratio.");
  desc.add<double>("explVarRatioCut_energyBoundary", 50.)
      ->setComment("Boundary energy between low and high energy explVarRatio cut threshold");
  desc.add<double>("explVarRatioMinimum_lowEnergy", 0.85)
      ->setComment(
          "Cut on explained variance ratio of tracksters to be considered as candidate, "
          "for trackster raw_energy < explVarRatioCut_energyBoundary");
  desc.add<double>("explVarRatioMinimum_highEnergy", 0.9)
      ->setComment(
          "Cut on explained variance ratio of tracksters to be considered as candidate, "
          "for trackster raw_energy > explVarRatioCut_energyBoundary");
  descriptions.add("superclusteringSampleDumper", desc);
}

DEFINE_FWK_MODULE(SuperclusteringSampleDumper);
