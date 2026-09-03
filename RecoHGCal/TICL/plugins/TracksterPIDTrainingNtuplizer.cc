#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

#include <TTree.h>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/HGCalObjects/interface/TICLGeomHost.h"
#include "CondFormats/HGCalObjects/interface/TICLGeomLayersHost.h"
#include "CondFormats/HGCalObjects/interface/TICLGeomLookupHost.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoHGCal/TICL/interface/TracksterPIDFeatures.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/TICLGeomTools.h"

class TracksterPIDTrainingNtuplizer
    : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns> {
public:
  explicit TracksterPIDTrainingNtuplizer(edm::ParameterSet const&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void beginJob() override;
  void beginRun(edm::Run const&, edm::EventSetup const&) override;
  void endRun(edm::Run const&, edm::EventSetup const&) override {}
  void analyze(edm::Event const&, edm::EventSetup const&) override;

  edm::EDGetTokenT<ticl::TracksterCollection> const trackstersToken_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> const layerClustersToken_;
  edm::EDGetTokenT<std::vector<int>> const signalMaskToken_;
  edm::EDGetTokenT<std::vector<int>> const fakeMaskToken_;
  edm::ESGetToken<TICLGeomHost, CaloGeometryRecord> const ticlGeomToken_;
  edm::ESGetToken<TICLGeomLookupHost, CaloGeometryRecord> const ticlGeomLookupToken_;
  edm::ESGetToken<TICLGeomLayersHost, CaloGeometryRecord> const ticlGeomLayersToken_;

  double const minTracksterEnergy_;
  unsigned int const backgroundPrescale_;
  bool saveSummary_ = false;
  bool saveLayer_ = false;
  bool saveCluster_ = false;
  bool saveImage_ = false;
  ticl::TracksterPIDFeatures const clusterFeaturesBuilder_;
  ticl::TracksterPIDFeatures const imageFeaturesBuilder_;
  ticlgeom::Tools rhtools_;

  TTree* tree_ = nullptr;
  TTree* eventTree_ = nullptr;
  uint32_t run_ = 0;
  uint32_t luminosityBlock_ = 0;
  uint64_t event_ = 0;
  uint32_t tracksterIndex_ = 0;
  int32_t label_ = -1;
  int32_t signalMask_ = -1;
  int32_t fakeMask_ = -1;
  float rawEnergy_ = 0.f;
  float rawEmEnergy_ = 0.f;
  float eta_ = 0.f;
  float phi_ = 0.f;
  uint32_t nVertices_ = 0;
  uint32_t nTracksters_ = 0;
  uint32_t nSignal_ = 0;
  uint32_t nBackground_ = 0;
  uint32_t nAmbiguous_ = 0;
  uint32_t nBelowEnergy_ = 0;
  std::vector<float> summaryFeatures_;
  std::vector<float> layerFeatures_;
  std::vector<float> clusterFeatures_;
  std::vector<float> imageFeatures_;
};

TracksterPIDTrainingNtuplizer::TracksterPIDTrainingNtuplizer(edm::ParameterSet const& config)
    : trackstersToken_(consumes(config.getParameter<edm::InputTag>("tracksters"))),
      layerClustersToken_(consumes(config.getParameter<edm::InputTag>("layerClusters"))),
      signalMaskToken_(consumes(config.getParameter<edm::InputTag>("signalMask"))),
      fakeMaskToken_(consumes(config.getParameter<edm::InputTag>("fakeMask"))),
      ticlGeomToken_(
          esConsumes<TICLGeomHost, CaloGeometryRecord, edm::Transition::BeginRun>(edm::ESInputTag("", ""))),
      ticlGeomLookupToken_(
          esConsumes<TICLGeomLookupHost, CaloGeometryRecord, edm::Transition::BeginRun>(edm::ESInputTag("", ""))),
      ticlGeomLayersToken_(
          esConsumes<TICLGeomLayersHost, CaloGeometryRecord, edm::Transition::BeginRun>(edm::ESInputTag("", ""))),
      minTracksterEnergy_(config.getParameter<double>("minTracksterEnergy")),
      backgroundPrescale_(config.getParameter<unsigned int>("backgroundPrescale")),
      clusterFeaturesBuilder_(config.getParameter<int>("nLayers"), config.getParameter<int>("maxClusters")),
      imageFeaturesBuilder_(config.getParameter<int>("nLayers"), config.getParameter<int>("maxClustersPerLayer")) {
  usesResource(TFileService::kSharedResource);
  if (minTracksterEnergy_ < 0. || backgroundPrescale_ == 0 || clusterFeaturesBuilder_.nLayers() <= 0 ||
      clusterFeaturesBuilder_.maxClusters() <= 0 || imageFeaturesBuilder_.maxClusters() <= 0) {
    throw cms::Exception("Configuration")
        << "minTracksterEnergy must be non-negative and all feature dimensions must be positive";
  }

  for (auto const& featureSet : config.getParameter<std::vector<std::string>>("featureSets")) {
    if (featureSet == "summary")
      saveSummary_ = true;
    else if (featureSet == "layer")
      saveLayer_ = true;
    else if (featureSet == "cluster")
      saveCluster_ = true;
    else if (featureSet == "image")
      saveImage_ = true;
    else
      throw cms::Exception("Configuration") << "Unknown featureSets entry '" << featureSet
                                            << "'; expected summary, layer, cluster, or image";
  }
  if (!(saveSummary_ || saveLayer_ || saveCluster_ || saveImage_))
    throw cms::Exception("Configuration") << "featureSets must contain at least one feature family";
}

void TracksterPIDTrainingNtuplizer::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("pidTraining", "Binary electromagnetic-versus-fake CLUE3D PID training samples");
  tree_->Branch("run", &run_);
  tree_->Branch("luminosityBlock", &luminosityBlock_);
  tree_->Branch("event", &event_);
  tree_->Branch("tracksterIndex", &tracksterIndex_);
  tree_->Branch("label", &label_);
  tree_->Branch("signalMask", &signalMask_);
  tree_->Branch("fakeMask", &fakeMask_);
  tree_->Branch("rawEnergy", &rawEnergy_);
  tree_->Branch("rawEmEnergy", &rawEmEnergy_);
  tree_->Branch("eta", &eta_);
  tree_->Branch("phi", &phi_);
  tree_->Branch("nVertices", &nVertices_);

  eventTree_ = fs->make<TTree>("pidTrainingEvents", "Per-event PID training sample bookkeeping");
  eventTree_->Branch("run", &run_);
  eventTree_->Branch("luminosityBlock", &luminosityBlock_);
  eventTree_->Branch("event", &event_);
  eventTree_->Branch("nTracksters", &nTracksters_);
  eventTree_->Branch("nSignal", &nSignal_);
  eventTree_->Branch("nBackground", &nBackground_);
  eventTree_->Branch("nAmbiguous", &nAmbiguous_);
  eventTree_->Branch("nBelowEnergy", &nBelowEnergy_);

  if (saveSummary_) {
    summaryFeatures_.resize(clusterFeaturesBuilder_.summarySize());
    tree_->Branch("summary", &summaryFeatures_);
  }
  if (saveLayer_) {
    layerFeatures_.resize(clusterFeaturesBuilder_.layerSize());
    tree_->Branch("layer", &layerFeatures_);
  }
  if (saveCluster_) {
    clusterFeatures_.resize(clusterFeaturesBuilder_.clusterSize());
    tree_->Branch("cluster", &clusterFeatures_);
  }
  if (saveImage_) {
    imageFeatures_.resize(imageFeaturesBuilder_.imageSize());
    tree_->Branch("image", &imageFeatures_);
  }
}

void TracksterPIDTrainingNtuplizer::beginRun(edm::Run const&, edm::EventSetup const& setup) {
  rhtools_.setGeometry(setup.getData(ticlGeomToken_),
                       setup.getData(ticlGeomLookupToken_),
                       setup.getData(ticlGeomLayersToken_));
}

void TracksterPIDTrainingNtuplizer::analyze(edm::Event const& event, edm::EventSetup const&) {
  auto const& tracksters = event.get(trackstersToken_);
  auto const& layerClusters = event.get(layerClustersToken_);
  auto const& signalMask = event.get(signalMaskToken_);
  auto const& fakeMask = event.get(fakeMaskToken_);
  if (signalMask.size() != tracksters.size() || fakeMask.size() != tracksters.size()) {
    throw cms::Exception("ProductMismatch") << "Trackster collection has " << tracksters.size()
                                            << " entries, signal mask has " << signalMask.size()
                                            << ", and fake mask has " << fakeMask.size();
  }

  run_ = event.id().run();
  luminosityBlock_ = event.id().luminosityBlock();
  event_ = event.id().event();
  nTracksters_ = tracksters.size();
  nSignal_ = 0;
  nBackground_ = 0;
  nAmbiguous_ = 0;
  nBelowEnergy_ = 0;
  std::vector<int> order;
  std::vector<int> seenClusters(imageFeaturesBuilder_.nLayers());

  for (std::size_t index = 0; index < tracksters.size(); ++index) {
    const bool isSignal = signalMask[index] == 0;
    const bool isFake = fakeMask[index] == 0;
    if (isSignal == isFake) {
      ++nAmbiguous_;
      continue;  // discard both the ambiguous score gap and any contradictory label
    }

    auto const& trackster = tracksters[index];
    float unsharedClusterEnergy = 0.f;
    for (auto vertex : trackster.vertices()) {
      unsharedClusterEnergy += layerClusters[vertex].energy();
      if (unsharedClusterEnergy >= minTracksterEnergy_)
        break;
    }
    if (unsharedClusterEnergy < minTracksterEnergy_) {
      ++nBelowEnergy_;
      continue;
    }
    if (isFake && backgroundPrescale_ > 1) {
      uint64_t key = event_ ^ (static_cast<uint64_t>(run_) << 32) ^
                     (static_cast<uint64_t>(luminosityBlock_) << 20) ^
                     (static_cast<uint64_t>(index) * 0x9E3779B185EBCA87ULL);
      key ^= key >> 30;
      key *= 0xBF58476D1CE4E5B9ULL;
      key ^= key >> 27;
      if (key % backgroundPrescale_ != 0)
        continue;
    }

    if (isSignal)
      ++nSignal_;
    else
      ++nBackground_;

    tracksterIndex_ = index;
    label_ = isSignal ? 1 : 0;
    signalMask_ = signalMask[index];
    fakeMask_ = fakeMask[index];
    rawEnergy_ = trackster.raw_energy();
    rawEmEnergy_ = trackster.raw_em_energy();
    eta_ = trackster.barycenter().eta();
    phi_ = trackster.barycenter().phi();
    nVertices_ = trackster.vertices().size();

    if (saveSummary_) {
      std::fill(summaryFeatures_.begin(), summaryFeatures_.end(), 0.f);
      clusterFeaturesBuilder_.fillSummary(summaryFeatures_.data(), trackster, layerClusters, rhtools_);
    }
    if (saveLayer_) {
      std::fill(layerFeatures_.begin(), layerFeatures_.end(), 0.f);
      clusterFeaturesBuilder_.fillLayer(layerFeatures_.data(), trackster, layerClusters, rhtools_);
    }
    if (saveCluster_) {
      std::fill(clusterFeatures_.begin(), clusterFeatures_.end(), 0.f);
      clusterFeaturesBuilder_.fillCluster(clusterFeatures_.data(), trackster, layerClusters, rhtools_, order);
    }
    if (saveImage_) {
      std::fill(imageFeatures_.begin(), imageFeatures_.end(), 0.f);
      imageFeaturesBuilder_.fillImage(
          imageFeatures_.data(), trackster, layerClusters, rhtools_, order, seenClusters);
    }
    tree_->Fill();
  }
  eventTree_->Fill();
}

void TracksterPIDTrainingNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracksters", edm::InputTag("ticlTrackstersCLUE3DHigh"));
  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalMergeLayerClusters"));
  desc.add<edm::InputTag>("signalMask", edm::InputTag("ticlValidSuperclusteringSeedMask"));
  desc.add<edm::InputTag>("fakeMask", edm::InputTag("ticlValidSuperclusteringSeedMask", "fakes"));
  desc.add<std::vector<std::string>>("featureSets", {"layer"})
      ->setComment("Branches to store: summary, layer, cluster, and/or image");
  desc.add<double>("minTracksterEnergy", 1.0)
      ->setComment("Same unshared layer-cluster energy eligibility cut used by PID inference");
  desc.add<unsigned int>("backgroundPrescale", 1)
      ->setComment("Deterministically retain one in N fake/background tracksters; signal is never prescaled");
  desc.add<int>("nLayers", 50);
  desc.add<int>("maxClusters", 128)->setComment("Top clusters retained for the cluster feature family");
  desc.add<int>("maxClustersPerLayer", 10)->setComment("Top clusters per layer retained for the image family");
  descriptions.add("tracksterPIDTrainingNtuplizer", desc);
}

DEFINE_FWK_MODULE(TracksterPIDTrainingNtuplizer);
