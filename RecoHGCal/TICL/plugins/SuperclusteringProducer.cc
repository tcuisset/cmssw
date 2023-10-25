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
  produces<SuperclusteringResult>("superclusteredTracksters"); // list of sets of indices of tracksters corresponding to a multicluster
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
  PlaceholderNNInput() {};

  static constexpr int featureCount{9};
  /*
  // SCaling values for network input, taken from Alessandro's notebook
  static constexpr std::array<float, featureCount>  scalingFactorTensor{{5.009924761018587,
  1.0018188597242355,
  0.010292174692060528,
  0.16288265085698453,
  0.04559109321673852,
  0.1682473240927026,
  0.15967887456557378,
  0.0015586544595969428,
  0.010722259319073196}};
  static constexpr std::array<float, featureCount> scalingMinTensor{{0.49977842782898685,
  0.4994941706919238,
  -0.020606080641338675,
  0.5002919753099604,
  -0.009400341457049127,
  0.49997065712986616,
  0.4996155459607123,
  -0.033067499716928135,
  -0.05139944510782373}};
  */
//Eigen::TensorFixedSize<float, Eigen::Sizes<featureCount>>
  void fillFeatures(tensorflow::Tensor& inputTensor, int batchIndex, Trackster const& ts_base, Trackster const& ts_toCluster) {
    //float deltaEta = (ts_toCluster.eigenvectors(0) - ts_base.eigenvectors(0)).Eta();
    /*  We use the barycenter for most of the variables below as that is what seems to have been used by Alessandro Tarabini, 
      but using PCA might be better. 
     (It would need retraining of the DNN)
    */
    assert(inputTensor.dims() == 2 && inputTensor.dim_size(1) == 9);
    assert(batchIndex < inputTensor.dim_size(0));
    auto eigenTensor = inputTensor.tensor<float, 2>();
    eigenTensor(batchIndex, 0) = std::abs(ts_toCluster.barycenter().Eta() - ts_base.barycenter().Eta());; //DeltaEtaBaryc
    eigenTensor(batchIndex, 1) = std::abs(ts_toCluster.barycenter().Phi() - ts_base.barycenter().phi());; //DeltaPhiBaryc
    eigenTensor(batchIndex, 2) = ts_toCluster.raw_energy(); //multi_en
    eigenTensor(batchIndex, 3) = ts_toCluster.barycenter().Eta(); //multi_eta
    eigenTensor(batchIndex, 4) = (ts_toCluster.raw_energy() * std::sin(ts_toCluster.barycenter().Theta())); //multi_pt
    eigenTensor(batchIndex, 5) = ts_base.barycenter().Eta(); //seedEta
    eigenTensor(batchIndex, 6) = ts_base.barycenter().Phi(); //seedPhi
    eigenTensor(batchIndex, 7) = ts_base.raw_energy(); //seedEn
    eigenTensor(batchIndex, 8) = (ts_base.raw_energy() * std::sin(ts_toCluster.barycenter().Theta())); //seedPt

    // TODO replace this with proper matrix multiplication (probably outside trackster loop)
    // or better : put inside tensorflow graph
    /*for (int i = 0; i < featureCount; i++) {
      eigenTensor(batchIndex, i) *= scalingFactorTensor[i];
      eigenTensor(batchIndex, i) += scalingMinTensor[i];
    }*/
  }
};

void SuperclusteringProducer::produce(edm::Event &evt, const edm::EventSetup &es) {
  edm::Handle<std::vector<Trackster>> inputTracksters;
  evt.getByToken(tracksters_clue3d_token_, inputTracksters);

  tensorflow::Session const* tfSession = es.getData(tfDnnToken_).getSession();
  
  const std::size_t tracksterCount = inputTracksters->size();
  const int mainBatchSize = tracksterCount * (tracksterCount - 1);

  PlaceholderNNInput nnInput;
  const int feature_count = nnInput.featureCount;
  tensorflow::Tensor inputTensor(tensorflow::DT_FLOAT, {mainBatchSize, feature_count});

  //Sorting tracksters by decreasing order of pT
  std::vector<std::size_t> trackstersIndicesPt(inputTracksters->size()); // Vector of indices into inputTracksters, sorted by decreasing order of pt
  std::iota(trackstersIndicesPt.begin(), trackstersIndicesPt.end(), 0);
  std::stable_sort(trackstersIndicesPt.begin(), trackstersIndicesPt.end(), [&inputTracksters](std::size_t i1, std::size_t i2) {
    return (*inputTracksters)[i1].raw_pt() > (*inputTracksters)[i2].raw_pt();
  });

  for (std::size_t ts_base_idx = 0; ts_base_idx < tracksterCount; ts_base_idx++) {
    Trackster const& ts_base = (*inputTracksters)[trackstersIndicesPt[ts_base_idx]];
    std::size_t ts_toCluster_idx_tensor = 0; // index of trackster in tensor, possibly shifted by one from ts_toCluster_idx
    for (std::size_t ts_toCluster_idx = 0; ts_toCluster_idx < tracksterCount; ts_toCluster_idx++) {
        
        if (trackstersIndicesPt[ts_base_idx] != ts_toCluster_idx) { // Don't supercluster trackster with itself
            Trackster const& ts_toCluster = (*inputTracksters)[ts_toCluster_idx]; // no need to sort by pt here

            nnInput.fillFeatures(inputTensor, ts_base_idx*(tracksterCount-1) + ts_toCluster_idx_tensor, ts_base, ts_toCluster);
            ts_toCluster_idx_tensor++; 
        }
    }
  }
  LogDebug("HGCalTICLSuperclustering") << "Input tensor " << inputTensor.SummarizeValue(100);

  /* Evaluate in minibatches since running with trackster count = 3000 leads to a short-lived ~15GB memory allocation
  As 3000^2 = 9e6 using a 
  */
  const int miniBatchSize = 1e6;
  std::vector<tensorflow::Tensor> miniBatchOutputs;
  for (int miniBatchStart = 0; miniBatchStart < mainBatchSize; miniBatchStart += miniBatchSize) {
    tensorflow::Tensor miniBatchInput = inputTensor.Slice(miniBatchStart, std::min(miniBatchStart+miniBatchSize, mainBatchSize));

    std::vector<tensorflow::Tensor> outputs;
    tensorflow::run(tfSession, {{"input", miniBatchInput}}, {"output_squeeze"}, &outputs);
    assert(outputs.size() == 1);
    miniBatchOutputs.push_back(std::move(outputs[0]));
  }

  // Helper function to access DNN results abstracting away the minibatches
  auto accessOutputTensor = [&miniBatchOutputs](int batchNumber) -> float {
    assert(static_cast<std::size_t>(batchNumber / miniBatchSize) < miniBatchOutputs.size());
    assert(miniBatchOutputs.at(batchNumber / miniBatchSize).dims() == 1);
    assert(static_cast<int>(batchNumber % miniBatchSize) < miniBatchOutputs.at(batchNumber / miniBatchSize).dim_size(0));

    return miniBatchOutputs.at(batchNumber / miniBatchSize).tensor<float, 1>()(batchNumber % miniBatchSize);
  };
  LogDebug("HGCalTICLSuperclustering") << "First output tensor " << miniBatchOutputs.at(0).SummarizeValue(100);

  // Build mask of tracksters already superclustered
  std::vector<bool> tracksterMask(tracksterCount, false); // Mask of tracksters, indexed with same indices as inputTracksters

  auto outputSuperclusters = std::make_unique<SuperclusteringResult>();
  // Supercluster tracksters
  for (std::size_t ts_base_idx = 0; ts_base_idx < tracksterCount; ts_base_idx++) {
    if (tracksterMask[trackstersIndicesPt[ts_base_idx]])
      continue;
    //Trackster const& ts_base = (*inputTracksters)[trackstersIndicesPt[ts_base_idx]];
    int ts_toCluster_idx_tensor = 0; // index of trackster in tensor, possibly shifted by one from ts_toCluster_idx

    // Build indices of tracksters in supercluster, starting with current trackster
    std::vector<std::size_t> superclusteredTracksterIndices{{trackstersIndicesPt[ts_base_idx]}};

    for (std::size_t ts_toCluster_idx = 0; ts_toCluster_idx < tracksterCount; ts_toCluster_idx++) {
        
        if (trackstersIndicesPt[ts_base_idx] != ts_toCluster_idx) { // Don't supercluster trackster with itself
            if (!tracksterMask[ts_toCluster_idx] // candidate trackster is not already part of a supercluster
                 // Candidate trackster passes the neural network working point
                 && accessOutputTensor(ts_base_idx*(tracksterCount-1) + ts_toCluster_idx_tensor) > nnWorkingPoint_) {
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

  evt.put(std::move(outputSuperclusters), "superclusteredTracksters");
}


void SuperclusteringProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("tfDnnLabel", "superclusteringTf");
  desc.add<edm::InputTag>("trackstersclue3d", edm::InputTag("ticlTrackstersCLUE3DHigh"));
  desc.add<double>("nnWorkingPoint", 0.51);
  descriptions.add("superclusteringProducer", desc);
}

DEFINE_FWK_MODULE(SuperclusteringProducer);
