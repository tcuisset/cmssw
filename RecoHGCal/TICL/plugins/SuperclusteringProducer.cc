#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

#include <Math/Vector2D.h>
#include <Math/Vector3D.h>
#include <Math/Rotation3D.h>
#include <Math/VectorUtil.h>

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
  const std::string dnnVersion_; 
  float nnWorkingPoint_; // Working point for neural network (above this score, consider the trackster candidate for superclustering)
};


SuperclusteringProducer::SuperclusteringProducer(const edm::ParameterSet &ps)
    : tracksters_clue3d_token_(consumes<std::vector<Trackster>>(ps.getParameter<edm::InputTag>("trackstersclue3d"))),
      //clusters_token_(consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layer_clusters"))),
      tfDnnToken_(esConsumes(edm::ESInputTag("", ps.getParameter<std::string>("tfDnnLabel")))),
      dnnVersion_(ps.getParameter<std::string>("dnnVersion")),
      nnWorkingPoint_(ps.getParameter<double>("nnWorkingPoint")) {
  produces<SuperclusteringResult>("superclusteredTracksters"); // list of sets of indices of tracksters corresponding to a multicluster
}

void SuperclusteringProducer::beginJob() {}

void SuperclusteringProducer::endJob(){};

void SuperclusteringProducer::beginRun(edm::Run const &iEvent, edm::EventSetup const &es) {
};

class AbstractDNNInput {
public:
  virtual ~AbstractDNNInput() {};
  virtual int featureCount() const = 0;
  virtual void fillFeatures(tensorflow::Tensor& inputTensor, int batchIndex, Trackster const& ts_base, Trackster const& ts_toCluster) = 0;
};

/*
DNN by Alessandro Tarabini comes in 2 versions with the following observables : 
 'v1_dEtadPhi': ['DeltaEta', 'DeltaPhi', 'multi_en', 'multi_eta', 'multi_pt', 'seedEta','seedPhi','seedEn', 'seedPt'],
  'v2_dEtadPhi': ['DeltaEta', 'DeltaPhi', 'multi_en', 'multi_eta', 'multi_pt', 'seedEta','seedPhi','seedEn', 'seedPt', 'theta', 'theta_xz_seedFrame', 'theta_yz_seedFrame', 'theta_xy_cmsFrame', 'theta_yz_cmsFrame', 'theta_xz_cmsFrame', 'explVar', 'explVarRatio'],
*/

class DNNInputAlessandroV1 : public AbstractDNNInput {
public:
  int featureCount() const override { return 9; }

//Eigen::TensorFixedSize<float, Eigen::Sizes<featureCount>>
  void fillFeatures(tensorflow::Tensor& inputTensor, int batchIndex, Trackster const& ts_base, Trackster const& ts_toCluster) override {
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
  }
};

/*
//class CMSCoordinateTag {};
class SeedCoordinateTag {};

using Vector3Seed = ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<float>, SeedCoordinateTag>;
using Vector3CMS = ROOT::Math::XYZVectorF;
using Vector2Seed = ROOT::Math::DisplacementVector2D<ROOT::Math::Cartesian3D<float>, SeedCoordinateTag>;
using Vector2CMS = ROOT::Math::XYVectorF;
*/

template <class Vector1, class Vector2>
double CosTheta2D( const Vector1 &  v1, const Vector2  & v2) {
  double arg;
  double v1_r2 = v1.X()*v1.X() + v1.Y()*v1.Y();
  double v2_r2 = v2.X()*v2.X() + v2.Y()*v2.Y();
  double ptot2 = v1_r2*v2_r2;
  if(ptot2 <= 0) {
      arg = 0.0;
  }else{
      double pdot = v1.X()*v2.X() + v1.Y()*v2.Y();
      using std::sqrt;
      arg = pdot/sqrt(ptot2);
      if(arg >  1.0) arg =  1.0;
      if(arg < -1.0) arg = -1.0;
  }
  return arg;
}
template <class Vector1, class Vector2>
inline double Angle2D( const  Vector1 & v1, const Vector2 & v2) {
  return std::acos( CosTheta2D(v1, v2) );
}

class DNNInputAlessandroV2 : public AbstractDNNInput{
public:
  int featureCount() const override { return 17; }

//Eigen::TensorFixedSize<float, Eigen::Sizes<featureCount>>
  void fillFeatures(tensorflow::Tensor& inputTensor, int batchIndex, Trackster const& ts_base, Trackster const& ts_toCluster) override {
    /*  We use the barycenter for most of the variables below as that is what seems to have been used by Alessandro Tarabini, 
      but using PCA might be better. 
     (It would need retraining of the DNN)
    */
    assert(inputTensor.dims() == 2 && inputTensor.dim_size(1) == featureCount());
    assert(batchIndex < inputTensor.dim_size(0));
    auto eigenTensor = inputTensor.tensor<float, 2>();

    using ROOT::Math::VectorUtil::Angle;
    using ROOT::Math::XYZVectorF;
    using ROOT::Math::XYVectorF;
    XYZVectorF const& pca_base_cmsFrame(ts_base.eigenvectors(0));
    XYZVectorF const& pca_cand_cmsFrame(ts_toCluster.eigenvectors(0));
    XYZVectorF xs(pca_base_cmsFrame.Cross(XYZVectorF(0, 0, 1)).Unit());
    ROOT::Math::Rotation3D rot(
        xs, 
        xs.Cross(pca_base_cmsFrame).Unit(),
        pca_base_cmsFrame);
      
    XYZVectorF pca_cand_seedFrame = rot(pca_cand_cmsFrame); // seed coordinates

    eigenTensor(batchIndex, 0) = std::abs(ts_toCluster.barycenter().Eta() - ts_base.barycenter().Eta());; //DeltaEtaBaryc
    eigenTensor(batchIndex, 1) = std::abs(ts_toCluster.barycenter().Phi() - ts_base.barycenter().phi());; //DeltaPhiBaryc
    eigenTensor(batchIndex, 2) = ts_toCluster.raw_energy(); //multi_en
    eigenTensor(batchIndex, 3) = ts_toCluster.barycenter().Eta(); //multi_eta
    eigenTensor(batchIndex, 4) = (ts_toCluster.raw_energy() * std::sin(ts_toCluster.barycenter().Theta())); //multi_pt
    eigenTensor(batchIndex, 5) = ts_base.barycenter().Eta(); //seedEta
    eigenTensor(batchIndex, 6) = ts_base.barycenter().Phi(); //seedPhi
    eigenTensor(batchIndex, 7) = ts_base.raw_energy(); //seedEn
    eigenTensor(batchIndex, 8) = (ts_base.raw_energy() * std::sin(ts_toCluster.barycenter().Theta())); //seedPt
    eigenTensor(batchIndex, 9) = Angle(pca_cand_cmsFrame, pca_base_cmsFrame); // theta : angle between seed and candidate
    eigenTensor(batchIndex, 10) = Angle2D(XYVectorF(pca_cand_seedFrame.y(), pca_cand_seedFrame.z()), XYVectorF(0, 1)); // theta_xz_seedFrame
    eigenTensor(batchIndex, 11) = Angle2D(XYVectorF(pca_cand_seedFrame.y(), pca_cand_seedFrame.z()), XYVectorF(0, 1)); // theta_yz_seedFrame
    eigenTensor(batchIndex, 12) = Angle2D(XYVectorF(pca_cand_cmsFrame.x(), pca_cand_cmsFrame.y()), XYVectorF(pca_base_cmsFrame.x(), pca_base_cmsFrame.y())); // theta_xy_cmsFrame
    eigenTensor(batchIndex, 13) = Angle2D(XYVectorF(pca_cand_cmsFrame.y(), pca_cand_cmsFrame.z()), XYVectorF(pca_base_cmsFrame.y(), pca_base_cmsFrame.z())); // theta_yz_cmsFrame
    eigenTensor(batchIndex, 14) = Angle2D(XYVectorF(pca_cand_cmsFrame.x(), pca_cand_cmsFrame.z()), XYVectorF(pca_base_cmsFrame.x(), pca_base_cmsFrame.z())); // theta_xz_cmsFrame
    eigenTensor(batchIndex, 15) = ts_toCluster.eigenvalues()[0]; // explVar
    eigenTensor(batchIndex, 16) = ts_toCluster.eigenvalues()[0] / std::accumulate(std::begin(ts_toCluster.eigenvalues()), std::end(ts_toCluster.eigenvalues()), 0, std::plus<float>()); // explVarRatio
  }
};

void SuperclusteringProducer::produce(edm::Event &evt, const edm::EventSetup &es) {
  edm::Handle<std::vector<Trackster>> inputTracksters;
  evt.getByToken(tracksters_clue3d_token_, inputTracksters);

  tensorflow::Session const* tfSession = es.getData(tfDnnToken_).getSession();
  
  const std::size_t tracksterCount = inputTracksters->size();
  const int mainBatchSize = tracksterCount * (tracksterCount - 1);

  std::unique_ptr<AbstractDNNInput> nnInput;
  if (dnnVersion_ == "alessandro-v1")
    nnInput = std::make_unique<DNNInputAlessandroV1>();
  else if (dnnVersion_ == "alessandro-v2")
    nnInput = std::make_unique<DNNInputAlessandroV2>();

  const int feature_count = nnInput->featureCount();
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

            nnInput->fillFeatures(inputTensor, ts_base_idx*(tracksterCount-1) + ts_toCluster_idx_tensor, ts_base, ts_toCluster);
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
  desc.add<std::string>("dnnVersion", "alessandro-v1");
  desc.add<edm::InputTag>("trackstersclue3d", edm::InputTag("ticlTrackstersCLUE3DHigh"));
  desc.add<double>("nnWorkingPoint", 0.51);
  descriptions.add("superclusteringProducer", desc);
}

DEFINE_FWK_MODULE(SuperclusteringProducer);
