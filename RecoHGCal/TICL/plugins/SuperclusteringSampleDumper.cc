// Original Author:  Theo Cuisset
//         Created:  Nov 2023
/**\class SuperclusteringSampleDumper SuperclusteringSampleDumper.cc RecoHGCal/TICL/plugins/SuperclusteringSampleDumper.cc

 Description: Produce samples for electron superclustering DNN training in HGCAL

 Pairs of seed-candidate tracksters (in compatible geometric windows) are iterated over, in similar manner as in SuperclusteringProducer.
 For each of these pairs, the DNN features are computed and saved to a TTree. 
 Also saved is the best (=lowest) association score of the seed trackster with CaloParticles. The association score of the candidate trackster 
 with the same CaloParticle is also saved.
*/
#include <memory>
#include <algorithm>

#include <TTree.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/Associations/interface/TracksterToSimTracksterAssociator.h"

#include "RecoHGCal/TICL/plugins/SuperclusteringDNNInputs.h"


using namespace ticl;

class SuperclusteringSampleDumper : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SuperclusteringSampleDumper(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  const edm::EDGetTokenT<std::vector<Trackster>> tracksters_clue3d_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimCP_token_;
  //const edm::EDGetTokenT<hgcal::SimToRecoCollectionSimTracksters> tsSimToRecoCP_token_;
  float deltaEtaWindow_;
  float deltaPhiWindow_;
  float seedPtThreshold_;
  float candidateEnergyThreshold_;

  TTree* output_tree_;
  unsigned int eventNb_;
  std::unique_ptr<AbstractDNNInput> dnnInput_;
  std::vector<std::vector<float>> features_; // Outer index : feature number (split into branches), inner index : inference pair index
  std::vector<unsigned int> seedTracksterIdx_; // ID of seed trackster used for inference pair
  std::vector<unsigned int> candidateTracksterIdx_; // ID of candidate trackster used for inference pair
  std::vector<float> seedTracksterBestAssociationScore_; // Best association score of seed trackster (seedTracksterIdx) with CaloParticle
  std::vector<float> candidateTracksterAssociationScoreWithSeed_; // Association score of candidate trackster with the CaloParticle used for seedTracksterBestAssociationScore
};

SuperclusteringSampleDumper::SuperclusteringSampleDumper(const edm::ParameterSet& ps)
    : tracksters_clue3d_token_(consumes<std::vector<Trackster>>(ps.getParameter<edm::InputTag>("tracksters"))),
      tsRecoToSimCP_token_(
          consumes<hgcal::RecoToSimCollectionSimTracksters>(ps.getParameter<edm::InputTag>("recoToSimAssociatorCP"))),
      deltaEtaWindow_(ps.getParameter<double>("deltaEtaWindow")),
      deltaPhiWindow_(ps.getParameter<double>("deltaPhiWindow")),
      seedPtThreshold_(ps.getParameter<double>("seedPtThreshold")),
      candidateEnergyThreshold_(ps.getParameter<double>("candidateEnergyThreshold_")),
      eventNb_(0),
      dnnInput_(makeDNNInputFromString(ps.getParameter<std::string>("dnnVersion"))),
      features_(dnnInput_->featureCount()) {
  usesResource("TFileService");
}

void SuperclusteringSampleDumper::beginJob() {
  edm::Service<TFileService> fs;
  output_tree_ = fs->make<TTree>("superclusteringTraining", "Superclustering training samples");
  output_tree_->Branch("Event", &eventNb_);
  output_tree_->Branch("seedTracksterIdx", &seedTracksterIdx_);
  output_tree_->Branch("candidateTracksterIdx", &candidateTracksterIdx_);
  output_tree_->Branch("seedTracksterBestAssociationScore", &seedTracksterBestAssociationScore_);
  output_tree_->Branch("candidateTracksterAssociationScoreWithSeed", &candidateTracksterAssociationScoreWithSeed_);
  std::vector<std::string> featureNames = dnnInput_->featureNames();
  assert(featureNames.size() == dnnInput_->featureCount());
  for (unsigned int i = 0; i < dnnInput_->featureCount(); i++) {
    output_tree_->Branch(("feature_" + featureNames[i]).c_str(), &features_[i]);
  }
}


void SuperclusteringSampleDumper::analyze(const edm::Event& evt, const edm::EventSetup& iSetup) {
  eventNb_++;
  edm::Handle<std::vector<Trackster>> inputTracksters;
  evt.getByToken(tracksters_clue3d_token_, inputTracksters);

  edm::Handle<hgcal::RecoToSimCollectionSimTracksters> assoc_CP_recoToSim;
  evt.getByToken(tsRecoToSimCP_token_, assoc_CP_recoToSim);
  //edm::Handle<hgcal::RecoToSimCollectionSimTracksters> assoc_CP_simToReco;
  //evt.getByToken(tsSimToRecoCP_token_, assoc_CP_simToReco);

  //Sorting tracksters by decreasing order of pT (out-of-place sort)
  std::vector<unsigned int> trackstersIndicesPt(inputTracksters->size()); // Vector of indices into inputTracksters, to be sorted
  std::iota(trackstersIndicesPt.begin(), trackstersIndicesPt.end(), 0); // Fill trackstersIndicesPt with 0...tracksterCount-1
  std::stable_sort(trackstersIndicesPt.begin(), trackstersIndicesPt.end(), [&inputTracksters](unsigned int i1, unsigned int i2) {
    return (*inputTracksters)[i1].raw_pt() > (*inputTracksters)[i2].raw_pt();
  });

  // Order of loops are reversed compared to SuperclusteringProducer (here outer is seed, inner is candidate), for performance reasons. 
  // The same pairs seed-candidate should be present, just in a different order
  // First loop on seed tracksters
  for (unsigned int ts_seed_idx = 0; ts_seed_idx < inputTracksters->size(); ts_seed_idx++) {
    const unsigned int ts_seed_idx_input = trackstersIndicesPt[ts_seed_idx]; // Index of seed trackster in input collection (not in pT sorted collection)
    Trackster const& ts_seed = (*inputTracksters)[ts_seed_idx_input];

    if (ts_seed.raw_pt() < seedPtThreshold_)
        break; // All further seeds will have lower pT than threshold (due to pT sorting)

    // Find best associated CaloParticle
    auto seed_assocs = assoc_CP_recoToSim->find({inputTracksters, ts_seed_idx_input});
    if (seed_assocs == assoc_CP_recoToSim->end())
      continue; // No CaloParticle associations for the current trackster (should not happen in theory)
    
    // Best score is smallest score
    hgcal::RecoToSimCollectionSimTracksters::data_type const& assocWithBestScore = *std::min_element(seed_assocs->val.begin(), seed_assocs->val.end(), 
      [](hgcal::RecoToSimCollectionSimTracksters::data_type const& assoc_1, hgcal::RecoToSimCollectionSimTracksters::data_type const& assoc_2) {
        // assoc_* is of type : std::pair<edmRefIntoSimTracksterCollection, std::pair<sharedEnergy, associationScore>>
        return assoc_1.second.second < assoc_2.second.second; 
      });
    

    // Second loop on superclustering candidates tracksters
    // Look only at candidate tracksters with lower pT than the seed (so all pairs are only looked at once)
    for (unsigned int ts_cand_idx = ts_seed_idx+1; ts_cand_idx < inputTracksters->size(); ts_cand_idx++) {
      Trackster const& ts_cand = (*inputTracksters)[trackstersIndicesPt[ts_cand_idx]];
      // Check that the two tracksters are geometrically compatible for superclustering (using deltaEta, deltaPhi window)
      // There is no need to run inference for tracksters very far apart
      if (std::abs(ts_seed.barycenter().Eta() - ts_cand.barycenter().Eta()) < deltaEtaWindow_
          && deltaPhi(ts_seed.barycenter().Phi(), ts_cand.barycenter().Phi()) < deltaPhiWindow_
          &&  ts_cand.raw_energy() < candidateEnergyThreshold_) { 
        // Add to output
        std::vector<float> features = dnnInput_->computeVector(ts_seed, ts_cand);
        assert(features.size() == features_.size());
        for (unsigned int feature_idx = 0; feature_idx < features_.size(); feature_idx++) {
          features_[feature_idx].push_back(features[feature_idx]);
        }
        seedTracksterIdx_.push_back(trackstersIndicesPt[ts_seed_idx]);
        candidateTracksterIdx_.push_back(trackstersIndicesPt[ts_cand_idx]);

        // Compute the association score of the candidate trackster with the CaloParticle that has the best association with the seed trackster
        float cand_assocToSeedCP_score = 1.; // default value : 1 is the worst score (0 is best)
        // First find associated CaloParticles with candidate
        auto cand_assocCP = assoc_CP_recoToSim->find(edm::Ref<ticl::TracksterCollection>(inputTracksters, trackstersIndicesPt[ts_cand_idx]));
        
        if (cand_assocCP != assoc_CP_recoToSim->end()) {
          // find the association with 
          auto cand_assocWithSeedCP = std::find_if(cand_assocCP->val.begin(), cand_assocCP->val.end(), 
            [&assocWithBestScore](hgcal::RecoToSimCollectionSimTracksters::data_type const& assoc) {
              // assoc is of type : std::pair<edmRefIntoSimTracksterCollection, std::pair<sharedEnergy, associationScore>>
              return assoc.first == assocWithBestScore.first; 
            });
          if (cand_assocWithSeedCP != cand_assocCP->val.end()) {
            cand_assocToSeedCP_score = cand_assocWithSeedCP->second.second;
          }
        }
        seedTracksterBestAssociationScore_.push_back(assocWithBestScore.second.second);
        candidateTracksterAssociationScoreWithSeed_.push_back(cand_assocToSeedCP_score);
      }
    }
  }

  output_tree_->Fill();
  for (auto& feats : features_)
    feats.clear();
  seedTracksterIdx_.clear();
  candidateTracksterIdx_.clear();
  seedTracksterBestAssociationScore_.clear();
  candidateTracksterAssociationScoreWithSeed_.clear();
}



void SuperclusteringSampleDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracksters", edm::InputTag("ticlTrackstersCLUE3DHigh"))
    ->setComment("Input trackster collection. Should be CLUE3D tracksters.");
  desc.add<edm::InputTag>("recoToSimAssociatorCP",
                          edm::InputTag("tracksterSimTracksterAssociationLinkingbyCLUE3D", "recoToSim"));
  //desc.add<edm::InputTag>("simToRecoAssociatorCP",
  //                        edm::InputTag("tracksterSimTracksterAssociationLinkingbyCLUE3D", "simToReco"));
  desc.add<std::string>("dnnVersion", "alessandro-v2")
    ->setComment("DNN version tag. Can be alessandro-v1 or alessandro-v2");
  desc.add<double>("deltaEtaWindow", 0.1)
     ->setComment("Size of delta eta window to consider for superclustering. Seed-candidate pairs outside this window are not considered for DNN inference.");
  desc.add<double>("deltaPhiWindow", 0.5)
     ->setComment("Size of delta phi window to consider for superclustering. Seed-candidate pairs outside this window are not considered for DNN inference.");
  desc.add<double>("seedPtThreshold", 1.)
     ->setComment("Minimum transverse momentum of trackster to be considered as seed of a supercluster");
  desc.add<double>("candidateEnergyThreshold", 2.)
     ->setComment("Minimum energy of trackster to be considered as candidate for superclustering");
  descriptions.add("superclusteringSampleDumper", desc);
}

DEFINE_FWK_MODULE(SuperclusteringSampleDumper);