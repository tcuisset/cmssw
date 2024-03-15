// Author : Theo Cuisset - theo.cuisset@polytechnique.edu
// Date : 11/2023
/*
TICL plugin for electron superclustering in HGCAL using a DNN. 
DNN designed and trained by Alessandro Tarabini.

Inputs are CLUE3D EM tracksters. Outputs are superclusters (as vectors of IDs of trackster)
"Seed trackster" : seed of supercluster, always highest pT trackster of supercluster, normally should be an electron
"Candidate trackster" : trackster that is considered for superclustering with a seed

Algorithm description :
1) Tracksters are ordered by decreasing pT.
2) We iterate over candidate tracksters, then over seed tracksters with higher pT than the candidate.
   If the pair seed-candidate is in a compatible eta-phi window and passes some selections (seed pT, energy, etc), then we add the DNN features of the pair to a tensor for later inference.
3) We run the inference with the DNN on the pairs (in minibatches to reduce memory usage)
4) We iterate over candidate and seed pairs inference results. For each candidate, we take the seed for which the DNN score for the seed-candidate score is best.
   If the score is also above a working point, then we add the candidate to the supercluster of the seed, and mask the candidate so it cannot be considered as a seed further

The loop is first on candidate, then on seeds as it is more efficient for step 4 to find the best seed for each candidate.
If macro SUPERCLUSTERING_DNN_SAVESCORE is set, will also produce DNN score in the event (as ticl::SuperclusteringDNNScore)
*/

#include <string>
#include <memory>

#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"
#include "RecoHGCal/TICL/plugins/TracksterLinkingbySuperClustering.h"
#include "RecoHGCal/TICL/interface/TracksterLinkingAlgoBase.h"
#include "RecoHGCal/TICL/plugins/SuperclusteringDNNInputs.h"

using namespace ticl;

std::unique_ptr<AbstractDNNInput> makeDNNInputFromString(std::string dnnVersion) {
  if (dnnVersion == "alessandro-v1")
    return std::make_unique<DNNInputAlessandroV1>();
  else if (dnnVersion == "alessandro-v2")
    return std::make_unique<DNNInputAlessandroV2>();
  assert(false);
}

TracksterLinkingbySuperClustering::TracksterLinkingbySuperClustering(const edm::ParameterSet& ps, edm::ConsumesCollector iC, cms::Ort::ONNXRuntime const* onnxRuntime)
      : TracksterLinkingAlgoBase(ps, iC, onnxRuntime),
      dnnVersion_(ps.getParameter<std::string>("dnnVersion")),
      nnWorkingPoint_(ps.getParameter<double>("nnWorkingPoint")),
      deltaEtaWindow_(ps.getParameter<double>("deltaEtaWindow")),
      deltaPhiWindow_(ps.getParameter<double>("deltaPhiWindow")),
      seedPtThreshold_(ps.getParameter<double>("seedPtThreshold")),
      candidateEnergyThreshold_(ps.getParameter<double>("candidateEnergyThreshold"))
{
  assert(onnxRuntime_ && "TracksterLinkingbySuperClustering : ONNXRuntime was not provided, the model should have been set in onnxModelPath in the plugin config");
} 

void TracksterLinkingbySuperClustering::initialize(const HGCalDDDConstants *hgcons,
                                             const hgcal::RecHitTools rhtools,
                                             const edm::ESHandle<MagneticField> bfieldH,
                                             const edm::ESHandle<Propagator> propH) {
}

/**
 * resultTracksters : should be all EM tracksters (including those not in Superclusters)
 * outputSuperclusters : indices into resultsTracksters. Should include ts not in SC as one-element vectors
 * linkedTracksterIdToInputTracksterId : map indices from output to input
*/
void TracksterLinkingbySuperClustering::linkTracksters(const Inputs& input, std::vector<Trackster>& resultTracksters,
                    std::vector<std::vector<unsigned int>>& outputSuperclusters,
                    std::vector<std::vector<unsigned int>>& linkedTracksterIdToInputTracksterId) {
  // For now we use all input tracksters for superclustering. At some point there might be a filter here for EM tracksters (electromagnetic identification with DNN ?)
  //resultTracksters = std::vector<Trackster>(input.tracksters.begin(), input.tracksters.end());
  auto const& inputTracksters = input.tracksters;
  const unsigned int tracksterCount = inputTracksters.size();

  std::unique_ptr<AbstractDNNInput> nnInput = makeDNNInputFromString(dnnVersion_);

  //Sorting tracksters by decreasing order of pT (out-of-place sort). 
  std::vector<unsigned int> trackstersIndicesPt(inputTracksters.size()); // Vector of indices into inputTracksters, to be sorted
  std::iota(trackstersIndicesPt.begin(), trackstersIndicesPt.end(), 0); // Fill trackstersIndicesPt with 0...tracksterCount-1
  std::stable_sort(trackstersIndicesPt.begin(), trackstersIndicesPt.end(), [&inputTracksters](unsigned int i1, unsigned int i2) {
    return inputTracksters[i1].raw_pt() > inputTracksters[i2].raw_pt();
  });

  /* Evaluate in minibatches since running with trackster count = 3000 leads to a short-lived ~15GB memory allocation
  Also we do not know in advance how many superclustering candidate pairs there are going to be
  So we pick (arbitrarily) a batch size of 10^5. It needs to be rounded to featureCount
  */
  const unsigned int miniBatchSize = static_cast<unsigned int>(1e5) / nnInput->featureCount() * nnInput->featureCount();

  std::vector<std::vector<float>> inputTensorBatches; // DNN input features tensors, in minibatches. Outer array : minibatches, inner array : 2D (flattened) array of features (indexed by batchIndex, featureId)
  // How far along in the latest tensor of inputTensorBatches are we. Set to miniBatchSize to trigger the creation of the tensor batch on first run
  unsigned int candidateIndexInCurrentBatch = miniBatchSize;
  // List of all (ts_seed_id; ts_cand_id) selected for DNN inference (same layout as inputTensorBatches) 
  // Index is in global trackster collection (not pt ordered collection)
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> tracksterIndicesUsedInDNN; 

  // Use TracksterTiles to speed up search of tracksters in eta-phi window. One per endcap
  std::array<TICLLayerTile, 2> tracksterTilesBothEndcaps; // one per endcap
  for (unsigned int i = 0; i < trackstersIndicesPt.size(); ++i)
  {
    Trackster const& ts = inputTracksters[trackstersIndicesPt[i]];
    tracksterTilesBothEndcaps[ts.barycenter().eta() > 0.].fill(ts.barycenter().eta(), ts.barycenter().phi(), i);
  }

  // First loop on candidate tracksters (start at 1 since the highest pt trackster can only be a seed, not a candidate)
  for (unsigned int ts_cand_idx = 1; ts_cand_idx < tracksterCount; ts_cand_idx++) {
    Trackster const& ts_cand = inputTracksters[trackstersIndicesPt[ts_cand_idx]];

    // Cut on candidate trackster energy, to match what was used for training by Alessandro
    if (ts_cand.raw_energy() < candidateEnergyThreshold_)
      continue;
    
    /* Cut on explained variance ratio. The DNN was trained by Alessandro Tarabini using this cut on the explained variance ratio.
    Therefore we reproduce it here. It is expected that this cut will be removed when the network for EM/hadronic differentiation is in place
    (would need retraining of the superclustering DNN) */
    float explVar_denominator = std::accumulate(std::begin(ts_cand.eigenvalues()), std::end(ts_cand.eigenvalues()), 0.f, std::plus<float>());
    if (explVar_denominator != 0.) {
      float explVarRatio = ts_cand.eigenvalues()[0] / explVar_denominator;
      if ((ts_cand.raw_energy() > 50 && explVarRatio <= 0.95) || (ts_cand.raw_energy() <= 50 && explVarRatio <= 0.92))
        continue;
    } else
      continue;

    auto& tracksterTiles = tracksterTilesBothEndcaps[ts_cand.barycenter().eta()>0];
    std::array<int, 4> search_box = tracksterTiles.searchBoxEtaPhi(
      ts_cand.barycenter().Eta() - deltaEtaWindow_,
      ts_cand.barycenter().Eta() + deltaEtaWindow_,
      ts_cand.barycenter().Phi() - deltaPhiWindow_,
      ts_cand.barycenter().Phi() + deltaPhiWindow_
    );
    // Second loop on superclustering seed tracksters
    // Use the search box for performance
    for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
      for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
        for (unsigned int ts_seed_idx : tracksterTiles[tracksterTiles.globalBin(eta_i, (phi_i % TileConstants::nPhiBins))]) {
          if (ts_seed_idx >= ts_cand_idx)
            continue; // Look only at seed tracksters with higher pT than the candidate (so all pairs are only looked at once)
          
          Trackster const& ts_seed = inputTracksters[trackstersIndicesPt[ts_seed_idx]];

          if (ts_seed.raw_pt() < seedPtThreshold_)
            break; // All further seeds will have lower pT than threshold (due to pT sorting)

          // Check that the two tracksters are geometrically compatible for superclustering
          if (std::abs(ts_seed.barycenter().Eta() - ts_cand.barycenter().Eta()) < deltaEtaWindow_
              && deltaPhi(ts_seed.barycenter().Phi(), ts_cand.barycenter().Phi()) < deltaPhiWindow_) { 
            if (candidateIndexInCurrentBatch >= miniBatchSize) {
              // Create new minibatch
              assert(candidateIndexInCurrentBatch == miniBatchSize);

              //Estimate how many seed-candidate pairs are remaining and don't allocate a full batch in this case. Use worst-case scenario of all pairs passing geometrical window
              // Also assume ts_seed_idx=0 (worst case) 
              // TODO we could probably use std::vector::reserve with a basic estimate for less memory usage on last batch
              inputTensorBatches.emplace_back(std::min(miniBatchSize, (tracksterCount*(tracksterCount-1) - ts_cand_idx*(ts_cand_idx-1))/2) * nnInput->featureCount());

              candidateIndexInCurrentBatch = 0;
              tracksterIndicesUsedInDNN.emplace_back();
            }

            std::vector<float> features = nnInput->computeVector(ts_seed, ts_cand); // Compute DNN features
            assert(features.size() == nnInput->featureCount());
            assert((candidateIndexInCurrentBatch+1)*nnInput->featureCount() <= inputTensorBatches.back().size());
            // Copy the features into the batch (TODO : could probably avoid the copy and fill straight in the batch vector)
            std::copy(features.begin(), features.end(), inputTensorBatches.back().begin() + candidateIndexInCurrentBatch*nnInput->featureCount());
            candidateIndexInCurrentBatch++;
            tracksterIndicesUsedInDNN.back().emplace_back(trackstersIndicesPt[ts_seed_idx], trackstersIndicesPt[ts_cand_idx]);
          }
        }
      }
    }
  }

  /* The last tensor of inputTensorBatches will be larger than necessary. The end of it will be uninitialized, then passed to the DNN
  We do not look at the output of the DNN on this so it should not matter.
  TODO consider slicing the last batch to avoid unneeded inference */

  if (inputTensorBatches.empty()) {
    LogDebug("HGCalTICLSuperclustering") << "No superclustering candidate pairs passed preselection before DNN. There are " << tracksterCount << " tracksters in this event.";
  }

#ifdef EDM_ML_DEBUG
  std::ostringstream s;
  for (unsigned int i = 0; i < std::min(nnInput->featureCount()*20, static_cast<unsigned int>(inputTensorBatches[0].size())); i++) {
    s << inputTensorBatches[0][i] << " ";
    if (i != 0 && i % nnInput->featureCount() == 0)
      s << "],\t[";
  }
  LogDebug("HGCalTICLSuperclustering") << inputTensorBatches.size() <<  " batches were created. First batch starts as follows : [" << s.str() << "]";
#endif

  // Run the DNN inference
  std::vector<std::vector<float>> batchOutputs; // Outer index : minibatch, inner index : inference index in minibatch, value : DNN score
  for (std::vector<float>& singleBatch : inputTensorBatches) {
    // ONNXRuntime takes std::vector<std::vector<float>>& as input (non-const reference) so we have to make a new vector
    std::vector<std::vector<float>> inputs_for_onnx{{std::move(singleBatch)}}; // Don't use singleBatch after this as it is in moved-from state
    std::vector<float> outputs = onnxRuntime_->run({"input"}, inputs_for_onnx, {}, {}, inputs_for_onnx[0].size()/nnInput->featureCount())[0];
    batchOutputs.push_back(std::move(outputs));
  }

#ifdef SUPERCLUSTERING_DNN_SAVESCORE
  auto outputTracksterDNNScore = std::make_unique<SuperclusteringDNNScore>();
#endif

  /* Build mask of tracksters already superclustered as candidates, as well as seeds (only needed to add tracksters not in a supercluster to the output).
  Uses global trackster ids
  */
  std::vector<bool> tracksterMask(tracksterCount, false);
  
  /* Index of the seed trackster of the previous iteration 
  Initialized with an id that cannot be obtained in input */
  unsigned int previousCandTrackster_idx = std::numeric_limits<unsigned int>::max(); 
  unsigned int bestSeedForCurrentCandidate_idx = std::numeric_limits<unsigned int>::max(); 
  float bestSeedForCurrentCandidate_dnnScore = nnWorkingPoint_;

  // Lambda to be called when there is a transition from one candidate to the next (as well as after the last iteration)
  // Does the actual supercluster creation
  auto onCandidateTransition = [&](unsigned ts_cand_idx) {
    if (bestSeedForCurrentCandidate_idx < std::numeric_limits<unsigned int>::max()) {
      // At least one seed can be superclustered with the candidate
      tracksterMask[ts_cand_idx] = true; // Mask the candidate so it is not considered as seed in later iterations

      // Look for a supercluster of the seed
      std::vector<std::vector<unsigned int>>::iterator seed_supercluster_it = std::find_if(outputSuperclusters.begin(), outputSuperclusters.end(), [bestSeedForCurrentCandidate_idx](std::vector<unsigned int> const& sc){
        return sc[0] == bestSeedForCurrentCandidate_idx;
      });

      if (seed_supercluster_it == outputSuperclusters.end()) {
        // No supercluster exists yet for the seed. Create one.
        outputSuperclusters.emplace_back();
        outputSuperclusters.back().push_back(bestSeedForCurrentCandidate_idx);
        resultTracksters.emplace_back(inputTracksters[bestSeedForCurrentCandidate_idx]);
        linkedTracksterIdToInputTracksterId.emplace_back(std::initializer_list<unsigned int>{bestSeedForCurrentCandidate_idx});
        seed_supercluster_it = outputSuperclusters.end()-1;
        tracksterMask[bestSeedForCurrentCandidate_idx] = true; // mask the seed as well (needed to find tracksters not in any supercluster)
      }
      // Index of the supercluster into resultTracksters, outputSuperclusters and linkedTracksterIdToInputTracksterId collections (the indices are the same)
      unsigned int indexIntoOutputTracksters = seed_supercluster_it - outputSuperclusters.begin();
      seed_supercluster_it->push_back(ts_cand_idx);
      resultTracksters[indexIntoOutputTracksters].mergeTracksters({inputTracksters[ts_cand_idx]});
      linkedTracksterIdToInputTracksterId[indexIntoOutputTracksters].push_back(ts_cand_idx);
      
      assert(outputSuperclusters.size() == resultTracksters.size() && outputSuperclusters.size() == linkedTracksterIdToInputTracksterId.size());
      assert(seed_supercluster_it->size() == linkedTracksterIdToInputTracksterId[indexIntoOutputTracksters].size());
      
      bestSeedForCurrentCandidate_idx = std::numeric_limits<unsigned int>::max(); 
      bestSeedForCurrentCandidate_dnnScore = nnWorkingPoint_;
    }
  };

  //Iterate over minibatches
  for (unsigned int batchIndex = 0; batchIndex < batchOutputs.size(); batchIndex++) {
    std::vector<float> const& currentBatchOutputs = batchOutputs[batchIndex]; // DNN score outputs
    // Iterate over seed-candidate pairs inside current minibatch
    for (unsigned int indexInBatch = 0; indexInBatch < tracksterIndicesUsedInDNN[batchIndex].size(); indexInBatch++) {
      assert(indexInBatch < static_cast<unsigned int>(batchOutputs[batchIndex].size()));

      const unsigned int ts_seed_idx = tracksterIndicesUsedInDNN[batchIndex][indexInBatch].first;
      const unsigned int ts_cand_idx = tracksterIndicesUsedInDNN[batchIndex][indexInBatch].second;
      const float currentDnnScore = currentBatchOutputs[indexInBatch];

      if (previousCandTrackster_idx != std::numeric_limits<unsigned int>::max() && ts_cand_idx != previousCandTrackster_idx) {
        // There is a transition from one seed to the next (don't make a transition for the first iteration)
        onCandidateTransition(previousCandTrackster_idx);
      }
      #ifdef SUPERCLUSTERING_DNN_SAVESCORE
      // Map the DNN score from float in [0, 1] to unsigned short
      outputTracksterDNNScore->emplace_back(ts_seed_idx, ts_cand_idx, static_cast<ticl::SuperclusteringDNNScoreValuePacked>(currentDnnScore*std::numeric_limits<SuperclusteringDNNScoreValuePacked>::max()));
      #endif

      if (currentDnnScore > bestSeedForCurrentCandidate_dnnScore && !tracksterMask[ts_seed_idx]) {
        // Check that the DNN suggests superclustering, that this seed-candidate assoc is better than previous ones, and that the seed is not already in a supercluster as candidate
        bestSeedForCurrentCandidate_idx = ts_seed_idx;
        bestSeedForCurrentCandidate_dnnScore = currentDnnScore;
      }
      previousCandTrackster_idx = ts_cand_idx;
    }
  }
  onCandidateTransition(previousCandTrackster_idx);

  // Adding one-trackster superclusters for all tracksters not in a supercluster already
  for (unsigned int ts_id = 0; ts_id < tracksterCount; ts_id++) {
    if (!tracksterMask[ts_id] // && inputTracksters[ts_id].raw_energy() >= candidateEnergyThreshold_
          && inputTracksters[ts_id].raw_pt() >= seedPtThreshold_
    ) {
      outputSuperclusters.emplace_back(std::initializer_list<unsigned int>{ts_id});
      resultTracksters.emplace_back(inputTracksters[ts_id]);
      linkedTracksterIdToInputTracksterId.emplace_back(std::initializer_list<unsigned int>{ts_id});
    } 
  }

  #ifdef EDM_ML_DEBUG
  for (std::vector<unsigned int> const& sc : outputSuperclusters) {
    std::ostringstream s;
    for (unsigned int trackster_id : sc)
      s << trackster_id << " ";
    LogDebug("HGCalTICLSuperclustering") << "Created supercluster of size " << sc.size() << " holding tracksters (first one is seed) " << s.str();
  }
  #endif
}

void TracksterLinkingbySuperClustering::fillPSetDescription(edm::ParameterSetDescription &desc) {
  TracksterLinkingAlgoBase::fillPSetDescription(desc); // adds algo_verbosity
  desc.add<edm::FileInPath>("onnxModelPath", edm::FileInPath("RecoHGCal/TICL/data/tf_models/supercls_v2.onnx"))
    ->setComment("Path to DNN (as ONNX model)");
  desc.add<std::string>("dnnVersion", "alessandro-v2")
    ->setComment("DNN version tag. Can be alessandro-v1 or alessandro-v2");
  desc.add<double>("nnWorkingPoint", 0.51)
     ->setComment("Working point of DNN (in [0, 1]). DNN score above WP will attempt to supercluster.");
  desc.add<double>("deltaEtaWindow", 0.1)
     ->setComment("Size of delta eta window to consider for superclustering. Seed-candidate pairs outside this window are not considered for DNN inference.");
  desc.add<double>("deltaPhiWindow", 0.5)
     ->setComment("Size of delta phi window to consider for superclustering. Seed-candidate pairs outside this window are not considered for DNN inference.");
  desc.add<double>("seedPtThreshold", 5.)
     ->setComment("Minimum transverse momentum of trackster to be considered as seed of a supercluster");
  desc.add<double>("candidateEnergyThreshold", 2.) // set the same as Alessandro
     ->setComment("Minimum energy of trackster to be considered as candidate for superclustering");
}
