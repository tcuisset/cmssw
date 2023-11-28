#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

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

  void beginJob();
  void endJob();

  void beginRun(edm::Run const &iEvent, edm::EventSetup const &es) override;

private:
  const edm::EDGetTokenT<std::vector<Trackster>> tracksters_clue3d_token_;
  //const edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_token_;

  const edm::ESGetToken<TfGraphDefWrapper, TfGraphRecord> tfDnnToken_; // ES Token to obtain tensorflow session and graph
  float nnWorkingPoint_; // Working point for neural network (above this score, consider the trackster candidate for superclustering)
};


SuperclusteringProducer::SuperclusteringProducer(const edm::ParameterSet &ps)
    : tracksters_clue3d_token_(consumes<std::vector<Trackster>>(ps.getParameter<edm::InputTag>("trackstersclue3d"))),
      //clusters_token_(consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layer_clusters"))),
      tfDnnToken_(esConsumes(edm::ESInputTag("", ps.getParameter<std::string>("tfDnnLabel")))),
      nnWorkingPoint_(ps.getParameter<double>("nnWorkingPoint")) {
  produces<std::vector<std::vector<std::size_t>>>("superclusteredTracksters"); // list of sets of indices of tracksters corresponding to a multicluster
}

void SuperclusteringProducer::beginJob() {}

void SuperclusteringProducer::endJob(){};

void SuperclusteringProducer::beginRun(edm::Run const &iEvent, edm::EventSetup const &es) {
};


/* class NNInput {
public:
  const int featureCount;
  void fillFeatures(Eigen::TensorMap<Eigen::Tensor<float, 2, 1, Eigen::DenseIndex>, 16> eigenTensor, int batchIndex, Trackster const& ts_base, Trackster const& ts_toCluster) = 0;
}; */

class PlaceholderNNInput{
public:
  PlaceholderNNInput() : featureCount(1) {}

  const int featureCount;

  void fillFeatures(Eigen::TensorMap<Eigen::Tensor<float, 2, 1, Eigen::DenseIndex>, 16> eigenTensor, int batchIndex, Trackster const& ts_base, Trackster const& ts_toCluster) {
    float deltaEta = (ts_toCluster.eigenvectors(0) - ts_base.eigenvectors(0)).Eta();
    eigenTensor(batchIndex, 0) = deltaEta; //dummy feature
  }
};

void SuperclusteringProducer::produce(edm::Event &evt, const edm::EventSetup &es) {
  edm::Handle<std::vector<Trackster>> inputTracksters;
  evt.getByToken(tracksters_clue3d_token_, inputTracksters);

  auto outputTracksters = std::make_unique<std::vector<Trackster>>();

  tensorflow::Session const* tfSession = es.getData(tfDnnToken_).getSession();
  
  const std::size_t tracksterCount = inputTracksters->size();
  const int batch_size = tracksterCount * (tracksterCount - 1);

  PlaceholderNNInput nnInput;
  const int feature_count = nnInput.featureCount;
  tensorflow::Tensor inputTensor(tensorflow::DT_FLOAT, {batch_size, feature_count});
  Eigen::TensorMap<Eigen::Tensor<float, 2, 1, Eigen::DenseIndex>, 16> eigenTensor = inputTensor.tensor<float, 2>();
  //for (Trackster const& ts_base : *inputTracksters) {
  //  for (Trackster const& ts_toCluster : *inputTracksters) {
  

  //Sorting tracksters by decreasing order of pT
  std::vector<std::size_t> trackstersIndicesPt(inputTracksters->size()); // Vector of indices into inputTracksters, sorted by decreasing order of pt
  std::iota(trackstersIndicesPt.begin(), trackstersIndicesPt.end(), 0);
  std::stable_sort(trackstersIndicesPt.begin(), trackstersIndicesPt.end(), [&inputTracksters](std::size_t i1, std::size_t i2) {
    return (*inputTracksters)[i1].raw_pt() > (*inputTracksters)[i2].raw_pt();
  });

  for (std::size_t ts_base_idx = 0; ts_base_idx < tracksterCount; ts_base_idx++) {
    Trackster const& ts_base = (*inputTracksters)[trackstersIndicesPt[ts_base_idx]];
    int ts_toCluster_idx_tensor = 0; // index of trackster in tensor, possibly shifted by one from ts_toCluster_idx
    for (std::size_t ts_toCluster_idx = 0; ts_toCluster_idx < tracksterCount; ts_toCluster_idx++) {
        
        if (ts_base_idx != ts_toCluster_idx) { // Don't supercluster trackster with itself
            Trackster const& ts_toCluster = (*inputTracksters)[ts_toCluster_idx]; // no need to sort by pt here

            nnInput.fillFeatures(eigenTensor, ts_base_idx*(tracksterCount-1) + ts_toCluster_idx_tensor, ts_base, ts_toCluster);
            ts_toCluster_idx_tensor++; 
        }
    }
  }
  LogDebug("HGCalTICLSuperclustering") << "Input tensor " << inputTensor.SummarizeValue(100);

  std::vector<tensorflow::Tensor> outputs;
  tensorflow::run(tfSession, {{"input_1", inputTensor}}, {"sequential/dense/Sigmoid"}, &outputs);
  assert(outputs.size() == 1);
  tensorflow::Tensor& outputTensor = outputs[0];
  const auto eigenOutputTensor = outputTensor.flat<float>();
  LogDebug("HGCalTICLSuperclustering") << "Output tensor " << outputTensor.SummarizeValue(100);

  // Build mask of tracksters already superclustered
  std::vector<bool> tracksterMask(tracksterCount, false); // Mask of tracksters, indexed with same indices as inputTracksters

  auto outputSuperclusters = std::make_unique<std::vector<std::vector<std::size_t>>>();
  // Supercluster tracksters
  for (std::size_t ts_base_idx = 0; ts_base_idx < tracksterCount; ts_base_idx++) {
    if (tracksterMask[trackstersIndicesPt[ts_base_idx]])
      continue;
    //Trackster const& ts_base = (*inputTracksters)[trackstersIndicesPt[ts_base_idx]];
    int ts_toCluster_idx_tensor = 0; // index of trackster in tensor, possibly shifted by one from ts_toCluster_idx

    // Build indices of tracksters in supercluster, starting with current trackster
    std::vector<std::size_t> superclusteredTracksterIndices{{trackstersIndicesPt[ts_base_idx]}};

    for (std::size_t ts_toCluster_idx = 0; ts_toCluster_idx < tracksterCount; ts_toCluster_idx++) {
        
        if (ts_base_idx != ts_toCluster_idx) { // Don't supercluster trackster with itself
            if (!tracksterMask[ts_toCluster_idx] // candidate trackster is not already part of a supercluster
                 // Candidate trackster passes the neural network working point
                 && eigenOutputTensor(ts_base_idx*(tracksterCount-1) + ts_toCluster_idx_tensor) > nnWorkingPoint_) {
                //float const& dnnScore = eigenOutputTensor(ts_base_idx*(tracksterCount-1) + ts_toCluster_idx_tensor);

                superclusteredTracksterIndices.push_back(ts_toCluster_idx);
                tracksterMask[ts_toCluster_idx] = true;

                LogDebug("HGCalTICLSuperclustering") << "Added trackster nb " << ts_toCluster_idx << " to supercluster (seed trackster nb " << trackstersIndicesPt[ts_base_idx] << ")";
            }
            ts_toCluster_idx_tensor++; 
        }
    }

    outputSuperclusters->push_back(std::move(superclusteredTracksterIndices));
  }

  evt.put(std::move(outputSuperclusters), "superclusteredTracksters"); // placeholder output
}


void SuperclusteringProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("tfDnnLabel", "superclusteringTf");
  desc.add<edm::InputTag>("trackstersclue3d", edm::InputTag("ticlTrackstersCLUE3DHigh"));
  desc.add<double>("nnWorkingPoint", 0.51);
  descriptions.add("superclusteringProducer", desc);
}

DEFINE_FWK_MODULE(SuperclusteringProducer);
