#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "PhysicsTools/TensorFlow/interface/TfGraphRecord.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "PhysicsTools/TensorFlow/interface/TfGraphDefWrapper.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/Common.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"

#include "RecoHGCal/TICL/interface/GlobalCache.h"

using namespace ticl;

class SuperclusteringProducer : public edm::stream::EDProducer<> {
public:
  explicit SuperclusteringProducer(const edm::ParameterSet &ps);
  ~SuperclusteringProducer() override{};
  void produce(edm::Event &, const edm::EventSetup &) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  // static methods for handling the global cache
  static std::unique_ptr<TrackstersCache> initializeGlobalCache(const edm::ParameterSet &);
  static void globalEndJob(TrackstersCache *);

  void beginJob();
  void endJob();

  void beginRun(edm::Run const &iEvent, edm::EventSetup const &es) override;

private:
  const edm::EDGetTokenT<std::vector<Trackster>> tracksters_clue3d_token_;
  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_token_;

  const edm::ESGetToken<TfGraphDefWrapper, TfGraphRecord> tfDnnToken_;
};



SuperclusteringProducer::SuperclusteringProducer(const edm::ParameterSet &ps)
    : tracksters_clue3d_token_(consumes<std::vector<Trackster>>(ps.getParameter<edm::InputTag>("trackstersclue3d"))),
      clusters_token_(consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layer_clusters"))),
      tfDnnToken_(esConsumes(edm::ESInputTag("", ps.getParameter<std::string>("tfDnnLabel")))) {
  produces<std::vector<Trackster>>("superclusteredTracksters");
}

void SuperclusteringProducer::beginJob() {}

void SuperclusteringProducer::endJob(){};

void SuperclusteringProducer::beginRun(edm::Run const &iEvent, edm::EventSetup const &es) {
};


void SuperclusteringProducer::produce(edm::Event &evt, const edm::EventSetup &es) {
  edm::Handle<std::vector<Trackster>> inputTracksters;
  evt.getByToken(tracksters_clue3d_token_, inputTracksters);

  auto outputTracksters = std::make_unique<std::vector<Trackster>>();

  tensorflow::Session const* tfSession = es.getData(tfDnnToken_).getSession();
  
  const int tracksterCount = inputTracksters->size();
  const int batch_size = tracksterCount * (tracksterCount - 1);
  const int feature_count = 1;
  tensorflow::Tensor tensor(tensorflow::DT_FLOAT, {batch_size, feature_count});
  Eigen::TensorMap<Eigen::Tensor<float, 2, 1, Eigen::DenseIndex>, 16> eigenTensor = tensor.tensor<float, 2>();
  //for (Trackster const& ts_base : *inputTracksters) {
  //  for (Trackster const& ts_toCluster : *inputTracksters) {
  for (int ts_base_idx = 0; ts_base_idx < tracksterCount; ts_base_idx++) {
    Trackster const& ts_base = (*inputTracksters)[ts_base_idx];
    int ts_toCluster_idx_tensor = 0; // index of trackster in tensor, possibly shifted by one from ts_toCluster_idx
    for (int ts_toCluster_idx = 0; ts_toCluster_idx < tracksterCount; ts_toCluster_idx++) {
        
        if (ts_base_idx != ts_toCluster_idx) { // Don't supercluster trackster with itself
            Trackster const& ts_toCluster = (*inputTracksters)[ts_toCluster_idx];
            eigenTensor(ts_base_idx*tracksterCount + ts_toCluster_idx_tensor, 0) = 1.; //dummy feature for now

            ts_toCluster_idx_tensor++; 
        }
    }
  }

  LogDebug("HGCalTICLSuperclustering") << tensor.SummarizeValue(100);

  evt.put(std::make_unique<std::vector<Trackster>>(2), "superclusteredTracksters"); // placeholder output
}


void SuperclusteringProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("tfDnnLabel", "tracksterSelectionTf"); // TODO maybe have a custom tensorflow session for superclustering ?
  descriptions.add("superclusteringProducer", desc);
}

DEFINE_FWK_MODULE(SuperclusteringProducer);
