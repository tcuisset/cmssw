#include <string>
#include <memory>

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

void TracksterLinkingbySuperClustering::linkTracksters(const Inputs& input, std::vector<Trackster>& resultTracksters,
                    std::vector<std::vector<unsigned int>>& outputSuperclusters,
                    std::vector<std::vector<unsigned int>>& linkedTracksterIdToInputTracksterId) {
  auto const& inputTracksters = input.tracksters;
  const unsigned int tracksterCount = inputTracksters.size();

  std::unique_ptr<AbstractDNNInput> nnInput = makeDNNInputFromString(dnnVersion_);

  //Sorting tracksters by decreasing order of pT (out-of-place sort)
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

  // First loop on candidate tracksters (start at 1 since the highest pt trackster can only be a seed, not a candidate)
  for (unsigned int ts_cand_idx = 1; ts_cand_idx < tracksterCount; ts_cand_idx++) {
    Trackster const& ts_cand = inputTracksters[trackstersIndicesPt[ts_cand_idx]];

    // Cut on candidate trackster energy, to match what was used for training by Alessandro
    if (ts_cand.raw_energy() < candidateEnergyThreshold_)
      continue;

    // Second loop on superclustering seed tracksters
    // Look only at seed tracksters with higher pT than the candidate (so all pairs are only looked at once)
    for (unsigned int ts_seed_idx = 0; ts_seed_idx < ts_cand_idx; ts_seed_idx++) {
      Trackster const& ts_seed = inputTracksters[trackstersIndicesPt[ts_seed_idx]];

      if (ts_seed.raw_pt() < seedPtThreshold_)
        break; // All further seeds will have lower pT than threshold (due to pT sorting)

      // Check that the two tracksters are geometrically compatible for superclustering (using deltaEta, deltaPhi window)
      // There is no need to run inference for tracksters very far apart
      if (std::abs(ts_seed.barycenter().Eta() - ts_cand.barycenter().Eta()) < deltaEtaWindow_
          && deltaPhi(ts_seed.barycenter().Phi(), ts_cand.barycenter().Phi()) < deltaPhiWindow_) { 
        /* Cut on explained variance ratio. The DNN was trained by Alessandro Tarabini using this cut on the explained variance ratio.
        Therefore we reproduce it here. It is expected that this cut will be removed when the network for EM/hadronic differentitation is in place
        (would need retraining of the superclustering DNN) */
        float explVar_denominator = std::accumulate(std::begin(ts_cand.eigenvalues()), std::end(ts_cand.eigenvalues()), 0.f, std::plus<float>());
        if (explVar_denominator != 0.) {
          float explVarRatio = ts_cand.eigenvalues()[0] / explVar_denominator; // explVarRatio
          if ((ts_cand.raw_energy() > 50 && explVarRatio <= 0.95) || (ts_cand.raw_energy() <= 50 && explVarRatio <= 0.92))
            continue;
        } else
          continue;

        // First check if we need to add an additional minibatch
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

  /* The last tensor of inputTensorBatches will be larger than necessary. The end of it will be uninitialized, then passed to the DNN
  We do not look at the output of the DNN on this so it should not matter.
  TODO consider slicing the last batch to avoid unneeded inference */

  if (inputTensorBatches.size() == 0) {
    LogDebug("HGCalTICLSuperclustering") << "No superclustering candidate pairs passed preselection before DNN. There are " << tracksterCount << " tracksters in this event.";
    //evt.put(std::make_unique<SuperclusteringResult>(), "superclusteredTracksters");
    //#ifdef SUPERCLUSTERING_DNN_SAVESCORE
    //evt.put(std::make_unique<SuperclusteringDNNScore>(), "superclusteringTracksterDNNScore");
    //#endif
    return;
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
    std::vector<float> outputs = onnxRuntime_.run({"input"}, inputs_for_onnx, {}, {}, inputs_for_onnx[0].size()/nnInput->featureCount())[0];
    batchOutputs.push_back(std::move(outputs));
  }

//  auto outputSuperclusters = std::make_unique<SuperclusteringResult>();
#ifdef SUPERCLUSTERING_DNN_SAVESCORE
  auto outputTracksterDNNScore = std::make_unique<SuperclusteringDNNScore>();
#endif

  // Build mask of tracksters already superclustered as candidates (seeds are not added). Uses global trackster ids
  std::vector<bool> tracksterMask(tracksterCount, false);
  // note that tracksterMask for the last seed trackster is never filled, as it is not needed.
  /* Index of the seed trackster of the previous iteration 
  Initialized with an id that cannot be obtained in input */
  unsigned int previousCandTrackster_idx = std::numeric_limits<unsigned int>::max(); 
  unsigned int bestSeedForCurrentCandidate_idx = std::numeric_limits<unsigned int>::max(); 
  float bestSeedForCurrentCandidate_dnnScore = nnWorkingPoint_;

  // Lambda to be called when there is a transition from one candidate to the next (as well as after the last iteration)
  auto onCandidateTransition = [&](unsigned ts_cand_idx) {
    if (bestSeedForCurrentCandidate_idx < std::numeric_limits<unsigned int>::max()) {
      // At least one seed can be superclustered with the candidate
      tracksterMask[ts_cand_idx] = true; // Mask the candidate so it is not considered as seed in later iterations

      // Look for a supercluster of the seed
      std::vector<std::vector<unsigned int>>::iterator seed_supercluster_it = std::find_if(outputSuperclusters.begin(), outputSuperclusters.end(), [bestSeedForCurrentCandidate_idx](std::vector<unsigned int> const& sc){
        return sc[0] == bestSeedForCurrentCandidate_idx;
      });

      //SuperclusteringResult::iterator seed_supercluster_it = std::find_if(outputSuperclusters->begin(), outputSuperclusters->end(), [bestSeedForCurrentCandidate_idx](Supercluster const& sc){
      //  return sc[0].key() == bestSeedForCurrentCandidate_idx;
      //});
      if (seed_supercluster_it == outputSuperclusters.end()) {
        // No supercluster exists yet for the seed. Create one.
        outputSuperclusters.emplace_back();
        outputSuperclusters.back().push_back(bestSeedForCurrentCandidate_idx);
        seed_supercluster_it = outputSuperclusters.end()-1;
      }
      seed_supercluster_it->push_back(ts_cand_idx);
      // Reset variables
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

  #ifdef EDM_ML_DEBUG
  /*
  for (std::vector<unsigned int> const& sc : *outputSuperclusters) {
    std::ostringstream s;
    for (auto trackster_it = sc.begin(); trackster_it != sc.end(); trackster_it++)
      s << trackster_it->key() << " ";
    LogDebug("HGCalTICLSuperclustering") << "Created supercluster of size " << sc.size() << " holding tracksters (first one is seed) " << s.str();
  }
  */
  #endif

  //evt.put(std::move(outputSuperclusters), "superclusteredTracksters");
  //#ifdef SUPERCLUSTERING_DNN_SAVESCORE
  //evt.put(std::move(outputTracksterDNNScore), "superclusteringTracksterDNNScore");
  //#endif
}

void TracksterLinkingbySuperClustering::fillPSetDescription(edm::ParameterSetDescription &desc) {
  TracksterLinkingAlgoBase::fillPSetDescription(desc); // adds algo_verbosity
  desc.add<edm::FileInPath>("dnn_model_path", edm::FileInPath("RecoHGCal/TICL/data/tf_models/supercls_v2.onnx"))
    ->setComment("Path to DNN (as ONNX model)");
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
  desc.add<double>("seedPtThreshold", 1.)
     ->setComment("Minimum transverse momentum of trackster to be considered as seed of a supercluster");
  desc.add<double>("candidateEnergyThreshold", 2.) // set the same as Alessandro
     ->setComment("Minimum energy of trackster to be considered as candidate for superclustering");
}
