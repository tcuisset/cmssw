// Author : Theo Cuisset - theo.cuisset@polytechnique.edu
// Date : 11/2023
/*
Electron superclustering in HGCAL using a DNN. 
DNN designed and trained by Alessandro Tarabini.

Consumes TICL tracksters, outputs superclusters (as vectors of IDs of tracksters)
Tracksters are ordered by decreasing pT. Then, we iterate over seed tracksters by decreasing pT. 
For each seed, we look at candidate tracksters with lower pT than the seed, that are in a compatible eta-phi window. 
We compute DNN inputs for every seed-candidate pairs thus found. Then the DNN inference is run (in batches to reduce memory usage).
DNN inference results are then iterated over in the same way as tracksters were. If the current seed-candidate pair has a DNN score above the working point,
the candidate is added to the supercluster of the seed (creating it if it does not exist yet), and the candidate is masked.
Masked candidates are not considered further (thus trackster are preferentially superclustered to higher pT seeds).

If macro EDM_ML_DEBUG is set, will also produce DNN score in the event (as ticl::SuperclusteringDNNScore)
*/
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort, std::copy
#include <sstream>      // std::stringstream

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

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/Supercluster.h"

using namespace ticl;

class SuperclusteringProducer : public edm::stream::EDProducer<> {
public:
  explicit SuperclusteringProducer(const edm::ParameterSet &ps);
  ~SuperclusteringProducer() override{};
  void produce(edm::Event &, const edm::EventSetup &) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  const edm::EDGetTokenT<std::vector<Trackster>> tracksters_clue3d_token_;
  const edm::ESGetToken<TfGraphDefWrapper, TfGraphRecord> tfDnnToken_; // ES Token to obtain tensorflow session and graph
  const std::string dnnVersion_; // Version identifier of the DNN (to choose which inputs to use)
  double nnWorkingPoint_; // Working point for neural network (above this score, consider the trackster candidate for superclustering)
  float deltaEtaWindow_; // Delta eta window to consider trackster seed-candidate pairs for inference
  float deltaPhiWindow_; // Delta phi window
};


SuperclusteringProducer::SuperclusteringProducer(const edm::ParameterSet &ps)
    : tracksters_clue3d_token_(consumes<std::vector<Trackster>>(ps.getParameter<edm::InputTag>("tracksters"))),
      tfDnnToken_(esConsumes(edm::ESInputTag("", ps.getParameter<std::string>("tfDnnLabel")))),
      dnnVersion_(ps.getParameter<std::string>("dnnVersion")),
      nnWorkingPoint_(ps.getParameter<double>("nnWorkingPoint")),
      deltaEtaWindow_(ps.getParameter<double>("deltaEtaWindow")),
      deltaPhiWindow_(ps.getParameter<double>("deltaPhiWindow")) {
  produces<SuperclusteringResult>("superclusteredTracksters"); // Produces std::vector<edm::RefVector<std::vector<Trackster>>>
#ifdef EDM_ML_DEBUG
  produces<SuperclusteringDNNScore>("superclusteringTracksterDNNScore"); // DNN scores of trackster -> trackster
#endif
}

// Abstract base class for DNN input preparation.
class AbstractDNNInput {
public:
  virtual ~AbstractDNNInput() {};
  virtual int featureCount() const = 0; // How many features the DNN uses as inputs
  /* Fill features in input tensor from pair of seed and candidate tracksters.
  inputTensor is a 2D tensor. This function should fill inputTensor(batchIndex, 0...featureCount) with features 
  */
  virtual void fillFeatures(tensorflow::Tensor& inputTensor, int batchIndex, Trackster const& ts_base, Trackster const& ts_toCluster) = 0;
};

/*
DNN by Alessandro Tarabini comes in 2 versions with the following observables : 
 'v1_dEtadPhi': ['DeltaEta', 'DeltaPhi', 'multi_en', 'multi_eta', 'multi_pt', 'seedEta','seedPhi','seedEn', 'seedPt'],
  'v2_dEtadPhi': ['DeltaEta', 'DeltaPhi', 'multi_en', 'multi_eta', 'multi_pt', 'seedEta','seedPhi','seedEn', 'seedPt', 'theta', 'theta_xz_seedFrame', 'theta_yz_seedFrame', 'theta_xy_cmsFrame', 'theta_yz_cmsFrame', 'theta_xz_cmsFrame', 'explVar', 'explVarRatio'],
*/

/* First version of DNN by Alessandro Tarabini. Meant as a DNN equivalent of Moustache algorithm (superclustering algo in ECAL)
Uses features : ['DeltaEta', 'DeltaPhi', 'multi_en', 'multi_eta', 'multi_pt', 'seedEta','seedPhi','seedEn', 'seedPt']
*/
class DNNInputAlessandroV1 : public AbstractDNNInput {
public:
  int featureCount() const override { return 9; }

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


// Helper functions for angles. Adapted from ROOT (3D vectors -> 2D vectors)
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

/* Second version of DNN by Alessandro Tarabini. 
Uses features : ['DeltaEta', 'DeltaPhi', 'multi_en', 'multi_eta', 'multi_pt', 'seedEta','seedPhi','seedEn', 'seedPt', 'theta', 'theta_xz_seedFrame', 'theta_yz_seedFrame', 'theta_xy_cmsFrame', 'theta_yz_cmsFrame', 'theta_xz_cmsFrame', 'explVar', 'explVarRatio']
*/
class DNNInputAlessandroV2 : public AbstractDNNInput{
public:
  int featureCount() const override { return 17; }

  void fillFeatures(tensorflow::Tensor& inputTensor, int batchIndex, Trackster const& ts_base, Trackster const& ts_toCluster) override {
    assert(inputTensor.dims() == 2 && inputTensor.dim_size(1) == featureCount());
    assert(batchIndex < inputTensor.dim_size(0));
    auto eigenTensor = inputTensor.tensor<float, 2>();

    using ROOT::Math::VectorUtil::Angle;
    using ROOT::Math::XYZVectorF;
    using ROOT::Math::XYVectorF;
    XYZVectorF const& pca_seed_cmsFrame(ts_base.eigenvectors(0));
    XYZVectorF const& pca_cand_cmsFrame(ts_toCluster.eigenvectors(0));
    XYZVectorF xs(pca_seed_cmsFrame.Cross(XYZVectorF(0, 0, 1)).Unit());
    ROOT::Math::Rotation3D rot(
        xs, 
        xs.Cross(pca_seed_cmsFrame).Unit(),
        pca_seed_cmsFrame);
      
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
    eigenTensor(batchIndex, 9) = Angle(pca_cand_cmsFrame, pca_seed_cmsFrame); // theta : angle between seed and candidate
    eigenTensor(batchIndex, 10) = Angle2D(XYVectorF(pca_cand_seedFrame.y(), pca_cand_seedFrame.z()), XYVectorF(0, 1)); // theta_xz_seedFrame
    eigenTensor(batchIndex, 11) = Angle2D(XYVectorF(pca_cand_seedFrame.y(), pca_cand_seedFrame.z()), XYVectorF(0, 1)); // theta_yz_seedFrame
    eigenTensor(batchIndex, 12) = Angle2D(XYVectorF(pca_cand_cmsFrame.x(), pca_cand_cmsFrame.y()), XYVectorF(pca_seed_cmsFrame.x(), pca_seed_cmsFrame.y())); // theta_xy_cmsFrame
    eigenTensor(batchIndex, 13) = Angle2D(XYVectorF(pca_cand_cmsFrame.y(), pca_cand_cmsFrame.z()), XYVectorF(pca_seed_cmsFrame.y(), pca_seed_cmsFrame.z())); // theta_yz_cmsFrame
    eigenTensor(batchIndex, 14) = Angle2D(XYVectorF(pca_cand_cmsFrame.x(), pca_cand_cmsFrame.z()), XYVectorF(pca_seed_cmsFrame.x(), pca_seed_cmsFrame.z())); // theta_xz_cmsFrame
    eigenTensor(batchIndex, 15) = ts_toCluster.eigenvalues()[0]; // explVar
    float explVar_denominator = std::accumulate(std::begin(ts_toCluster.eigenvalues()), std::end(ts_toCluster.eigenvalues()), 0, std::plus<float>());
    if (explVar_denominator != 0.)
      eigenTensor(batchIndex, 16) = ts_toCluster.eigenvalues()[0] / explVar_denominator; // explVarRatio
    else {
      eigenTensor(batchIndex, 16) = 0.;  // TODO Study what would be best in this case (or if it could ever happen)
      LogDebug("HGCalTICLSuperclustering") << "Sum of eigenvalues was zero for trackster. Could not compute explained variance ratio.";
    }
  }
};

void SuperclusteringProducer::produce(edm::Event &evt, const edm::EventSetup &es) {
  edm::Handle<std::vector<Trackster>> inputTracksters;
  evt.getByToken(tracksters_clue3d_token_, inputTracksters);
  const unsigned int tracksterCount = inputTracksters->size();

  tensorflow::Session const* tfSession = es.getData(tfDnnToken_).getSession();

  std::unique_ptr<AbstractDNNInput> nnInput;
  if (dnnVersion_ == "alessandro-v1")
    nnInput = std::make_unique<DNNInputAlessandroV1>();
  else if (dnnVersion_ == "alessandro-v2")
    nnInput = std::make_unique<DNNInputAlessandroV2>();
  const int feature_count = nnInput->featureCount();


  //Sorting tracksters by decreasing order of pT (out-of-place sort)
  std::vector<unsigned int> trackstersIndicesPt(inputTracksters->size()); // Vector of indices into inputTracksters, to be sorted
  std::iota(trackstersIndicesPt.begin(), trackstersIndicesPt.end(), 0); // Fill trackstersIndicesPt with 0...tracksterCount-1
  std::stable_sort(trackstersIndicesPt.begin(), trackstersIndicesPt.end(), [&inputTracksters](unsigned int i1, unsigned int i2) {
    return (*inputTracksters)[i1].raw_pt() > (*inputTracksters)[i2].raw_pt();
  });


  /* Evaluate in minibatches since running with trackster count = 3000 leads to a short-lived ~15GB memory allocation
  Also we do not know in advance how many superclustering candidate pairs there are going to be
  */
  const unsigned int miniBatchSize = 1e6;

  std::vector<tensorflow::Tensor> inputTensorBatches; // DNN input features tensors, in minibatches
  // How far along in the latest tensor of inputTensorBatches are we. Set to miniBatchSize to trigger the creation of the tensor batch on first run
  unsigned int candidateIndexInCurrentBatch = miniBatchSize;
  // List of all (ts_seed_id; ts_cand_id) selected for DNN inference (same layout as inputTensorBatches) 
  // Index is in global trackster collection (not pt ordered collection)
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> tracksterIndicesUsedInDNN; 

  // First loop on seed tracksters
  for (unsigned int ts_seed_idx = 0; ts_seed_idx < tracksterCount; ts_seed_idx++) {
    Trackster const& ts_seed = (*inputTracksters)[trackstersIndicesPt[ts_seed_idx]];

    // Second loop on superclustering candidates tracksters
    // Look only at candidate tracksters with lower pT than the seed (so all pairs are only looked at once)
    for (unsigned int ts_cand_idx = ts_seed_idx+1; ts_cand_idx < tracksterCount; ts_cand_idx++) {
      Trackster const& ts_cand = (*inputTracksters)[trackstersIndicesPt[ts_cand_idx]];
      // Check that the two tracksters are geometrically compatible for superclustering (using deltaEta, deltaPhi window)
      // There is no need to run inference for tracksters very far apart
      if (std::abs(ts_seed.barycenter().Eta() - ts_cand.barycenter().Eta()) < deltaEtaWindow_
          && deltaPhi(ts_seed.barycenter().Phi(), ts_cand.barycenter().Phi()) < deltaPhiWindow_) { 
        // First check if we need to add an additional minibatch
        if (candidateIndexInCurrentBatch >= miniBatchSize) {
          // Create new minibatch
          assert(candidateIndexInCurrentBatch == miniBatchSize);
          //Estimate how many seed-candidate pairs are remaining and don't allocate a full batch in this case
          // static_cast<long> is required to avoid narrowing error
          inputTensorBatches.emplace_back(tensorflow::DT_FLOAT, 
              tensorflow::TensorShape({static_cast<long>(std::min(miniBatchSize, (tracksterCount-ts_seed_idx) * (tracksterCount - ts_seed_idx - 1u)/2)), feature_count}));
          candidateIndexInCurrentBatch = 0;
          tracksterIndicesUsedInDNN.emplace_back();
        }

        nnInput->fillFeatures(inputTensorBatches.back(), candidateIndexInCurrentBatch, ts_seed, ts_cand);
        candidateIndexInCurrentBatch++;
        tracksterIndicesUsedInDNN.back().emplace_back(trackstersIndicesPt[ts_seed_idx], trackstersIndicesPt[ts_cand_idx]);
      }
    }
  }

  /* The last tensor of inputTensorBatches will be larger than necessary. The end of it will be uninitialized, then passed to the DNN
  We do not look at the output of the DNN on this so it should not matter.
  TODO consider slicing the last batch to avoid unneeded inference */

  if (inputTensorBatches.size() == 0) {
    LogDebug("HGCalTICLSuperclustering") << "No superclustering candidate pairs passed preselection before DNN. There are " << tracksterCount << " tracksters in this event.";
    evt.put(std::make_unique<SuperclusteringResult>(), "superclusteredTracksters");
    #ifdef EDM_ML_DEBUG
    evt.put(std::make_unique<SuperclusteringDNNScore>(), "superclusteringTracksterDNNScore");
    #endif
    return;
  }

  LogDebug("HGCalTICLSuperclustering") << "Input tensor " << inputTensorBatches[0].SummarizeValue(100);

  // Run the DNN inference
  std::vector<tensorflow::Tensor> batchOutputs;
  for (tensorflow::Tensor& singleBatch : inputTensorBatches) {
    std::vector<tensorflow::Tensor> outputs;
    tensorflow::run(tfSession, {{"input", singleBatch}}, {"output_squeeze"}, &outputs);
    assert(outputs.size() == 1);
    batchOutputs.push_back(std::move(outputs[0]));
  }

  auto outputSuperclusters = std::make_unique<SuperclusteringResult>();
#ifdef EDM_ML_DEBUG
  auto outputTracksterDNNScore = std::make_unique<SuperclusteringDNNScore>();
#endif

  // Build mask of tracksters already superclustered as candidates (seeds are not added). Uses global trackster ids
  std::vector<bool> tracksterMask(tracksterCount, false);
  // note that tracksterMask for the last seed trackster is never filled, as it is not needed.
  unsigned int previousSeedTrackster_idx; // Index of the seed trackster of the previous iteration
  ticl::Supercluster currentSupercluster; 
  //Iterate over minibatches
  for (unsigned int batchIndex = 0; batchIndex < batchOutputs.size(); batchIndex++) {
    auto outputEigenTensor = batchOutputs[batchIndex].tensor<float, 1>(); // 1D-tensor of DNN score outputs
    // Iterate over seed-candidate pairs inside current minibatch
    for (unsigned int indexInBatch = 0; indexInBatch < tracksterIndicesUsedInDNN[batchIndex].size(); indexInBatch++) {
      assert(indexInBatch < static_cast<unsigned int>(batchOutputs[batchIndex].dim_size(0)));

      unsigned int currentSeedTrackster_idx = tracksterIndicesUsedInDNN[batchIndex][indexInBatch].first;
      if (currentSeedTrackster_idx != previousSeedTrackster_idx) {
        // There is a transition from one seed to the next
        if (!currentSupercluster.empty()) {
          // Register the seed as part of a supercluster if at least one candidate got superclustered
          assert(!tracksterMask[previousSeedTrackster_idx]);

          #ifdef EDM_ML_DEBUG
          std::ostringstream s;
          for (ticl::Supercluster::const_iterator it = currentSupercluster.begin(); it != currentSupercluster.end(); it++)
            s << it->key() << " ";
          LogDebug("HGCalTICLSuperclustering") << "Created supercluster of size " << currentSupercluster.size() << " holding tracksters (first one is seed) " << s.str();
          #endif
          
          // tracksterMask[previousSeedTrackster_idx] = true; // no need to mask seeds as they will not be considered again in the loop
          outputSuperclusters->push_back(std::move(currentSupercluster));
          currentSupercluster.clear(); // bring the vector from moved-from state to empty state
        }
        previousSeedTrackster_idx = currentSeedTrackster_idx;
      }
      #ifdef EDM_ML_DEBUG
      // Map the DNN score from float in [0, 1] to unsigned short
      outputTracksterDNNScore->emplace_back(currentSeedTrackster_idx, tracksterIndicesUsedInDNN[batchIndex][indexInBatch].second, static_cast<ticl::SuperclusteringDNNScoreValuePacked>(outputEigenTensor(indexInBatch)*std::numeric_limits<SuperclusteringDNNScoreValuePacked>::max()));
      #endif

      if (outputEigenTensor(indexInBatch) > nnWorkingPoint_) {
        unsigned int ts_cand_idx = tracksterIndicesUsedInDNN[batchIndex][indexInBatch].second;
        if (!tracksterMask[currentSeedTrackster_idx] && !tracksterMask[ts_cand_idx]) { // Check that the seed and candidates are not in a supercluster
          // Note that the seed is not added to the trackster mask
          if (currentSupercluster.size() == 0) {
             // only add the seed once to the supercluster
            currentSupercluster.push_back({inputTracksters, currentSeedTrackster_idx});
          }
          currentSupercluster.push_back({inputTracksters, ts_cand_idx});
          tracksterMask[ts_cand_idx] = true; // Mask the candidate
        }
      }
    }
  }

  evt.put(std::move(outputSuperclusters), "superclusteredTracksters");
  #ifdef EDM_ML_DEBUG
  evt.put(std::move(outputTracksterDNNScore), "superclusteringTracksterDNNScore");
  #endif
}


void SuperclusteringProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("tfDnnLabel", "superclusteringTf")
    ->setComment("Name of an TfGraphDefWrapper (PhysicsTools/TensorFlow) in the event holding the DNN. It is ususally created by TfGraphDefProducer in RecoHGCal/TICL/python/superclusteringTf_cff.py");
  desc.add<std::string>("dnnVersion", "alessandro-v2")
    ->setComment("DNN version tag. Can be alessandro-v1 or alessandro-v2");
  desc.add<edm::InputTag>("tracksters", edm::InputTag("ticlTrackstersCLUE3DHigh"))
    ->setComment("Input trackster collection. Should be CLUE3D tracksters.");
  desc.add<double>("nnWorkingPoint", 0.51)
     ->setComment("Working point of DNN (in [0, 1]). DNN score above WP will attempt to supercluster.");
  desc.add<double>("deltaEtaWindow", 0.1)
     ->setComment("Size of delta eta window to consider for superclustering. Seed-candidate pairs outside this window are not considered for DNN inference.");
  desc.add<double>("deltaPhiWindow", 0.5)
     ->setComment("Size of delta phi window to consider for superclustering. Seed-candidate pairs outside this window are not considered for DNN inference.");
  descriptions.add("superclusteringProducer", desc);
}

DEFINE_FWK_MODULE(SuperclusteringProducer);
