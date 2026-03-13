#include "RecoHGCal/TICL/interface/TracksterInferenceByTransformer.h"
#include "RecoHGCal/TICL/interface/TracksterInferenceAlgoFactory.h"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

namespace ticl {

  TracksterInferenceByTransformer::TracksterInferenceByTransformer(const edm::ParameterSet& conf, TICLONNXGlobalCache const* cache)
      : TracksterInferenceAlgoBase(conf, cache),
        inputNames_(conf.getParameter<std::vector<std::string>>("inputNames")),
        // output_en_(conf.getParameter<std::vector<std::string>>("output_en")),
        output_id_(conf.getParameter<std::vector<std::string>>("output_id")),
        eidMinClusterEnergy_(conf.getParameter<double>("eid_min_cluster_energy")),
        // eidNLayers_(conf.getParameter<int>("eid_n_layers")),
        eidNClusters_(conf.getParameter<int>("eid_n_clusters")),
        doPID_(conf.getParameter<int>("doPID")),
        // doRegression_(conf.getParameter<int>("doRegression")),
        miniBatchSize_(conf.getUntrackedParameter<int>("miniBatchSize", 256)) {
    const std::string pidModel = conf.getParameter<std::string>("onnxPIDModelPath");
    const std::string energyModel = conf.getParameter<std::string>("onnxEnergyModelPath");

    if (cache_ != nullptr) {
      if (!pidModel.empty()) {
        onnxPIDSession_ = cache_->getByModelPathString(pidModel);
      }
    }

    enabled_ = ((doPID_ != 0 && onnxPIDSession_ != nullptr));

    ortScratch_.inputs.resize(3);
    ortScratch_.input_shapes.resize(3);
  }

  void TracksterInferenceByTransformer::runInference(const std::vector<reco::CaloCluster>& layerClusters,
                                             std::vector<Trackster>& tracksters,
                                             const hgcal::RecHitTools& rhtools) const {
    if (!enabled_ || tracksters.empty()) {
      return;
    }

    // ---- select tracksters (same physics logic), reset outputs once
    std::vector<int> indices;
    indices.reserve(tracksters.size());

    for (int i = 0; i < static_cast<int>(tracksters.size()); ++i) {
      float sumClusterEnergy = 0.f;

      // Note: keep the same semantics you had (skip barrel clusters, sum endcap energy)
      for (const unsigned int& v : tracksters[i].vertices()) {
        if (rhtools.isBarrel(layerClusters[v].seed())) {
          continue;
        }
        sumClusterEnergy += static_cast<float>(layerClusters[v].energy());
        if (sumClusterEnergy >= eidMinClusterEnergy_) {
          tracksters[i].setRegressedEnergy(0.f);
          tracksters[i].zeroProbabilities();
          indices.push_back(i);
          break;
        }
      }
    }

    const int total = static_cast<int>(indices.size());
    if (total == 0) {
      return;
    }

    const int mb = std::max(1, miniBatchSize_);

    // Reuse buffers across events
    ortScratch_.clearPerEvent();

    // // Per-minibatch reusable temporaries to avoid churn
    // std::vector<int> seenClusters;
    // seenClusters.resize(eidNLayers_);

    std::vector<int> clusterIndices;

    // Alias for input tensor
    std::vector<float> layerClusterFeatures;// = ortScratch_.inputs[0];
    std::vector<uint8_t> layerClusterMask; // = ortScratch_.inputs[1];
    std::vector<float> tracksterFeatures; // = ortScratch_.inputs[2];

    for (int start = 0; start < total; start += mb) {
      const int nTrackstersInBatch = std::min(mb, total - start); ///< nTracksters in batch

      // shape: layer cluster features = [trackster, maxLayerClusterCount, nFeats]
      ortScratch_.input_shapes[0] = {nTrackstersInBatch, eidNClusters_, eidNFeatures_};
      // shape: layer cluster mask = [trackster, maxLayerClusterCount]
      ortScratch_.input_shapes[1] = {nTrackstersInBatch, eidNClusters_};
      // shape: trackster features : [trackster, nTracksterFeatures]
      ortScratch_.input_shapes[2] = {nTrackstersInBatch, 5};

      const size_t nFloats = static_cast<size_t>(nTrackstersInBatch) * eidNClusters_ * eidNFeatures_;
      layerClusterFeatures.assign(nFloats, 0.f);  // sparse fill -> must zero

      layerClusterMask.assign(static_cast<size_t>(nTrackstersInBatch)*eidNClusters_, 0);
      tracksterFeatures.assign(static_cast<size_t>(nTrackstersInBatch)*5, 0.f);

      // ---- build sparse tensor for this minibatch
      for (int tracksterIdxInBatch = 0; tracksterIdxInBatch < nTrackstersInBatch; ++tracksterIdxInBatch) { // loop on trackster
        const int tsIdx = indices[start + tracksterIdxInBatch];
        Trackster const& ts = tracksters[tsIdx];

        const int vtxCount = static_cast<int>(ts.vertices().size());
        clusterIndices.resize(vtxCount);
        std::iota(clusterIndices.begin(), clusterIndices.end(), 0);

        // the layer cluster mask is 0 for LCs in use and padded with ones
        std::fill(layerClusterMask.begin()+vtxCount,layerClusterMask.end(), 1); 

        // std::sort(clusterIndices.begin(), clusterIndices.end(), [&layerClusters, &ts](int a, int b) {
        //   return layerClusters[ts.vertices(a)].energy() > layerClusters[ts.vertices(b)].energy();
        // });

        int minLayer = -100;
        int maxLayer = 100;

        int layerClusterCount = 0;
        for (int k : clusterIndices) {
          const unsigned int v = ts.vertices(k);
          auto const& cl = layerClusters[v];

          const int layer = rhtools.getLayerWithOffset(cl.hitsAndFractions()[0].first);
          minLayer = std::min(layer, minLayer);
          maxLayer = std::max(layer, maxLayer);

          if (layerClusterCount >= eidNClusters_)
            break;
          
          // trackster x LCs x features
          const size_t base =
              static_cast<size_t>(tracksterIdxInBatch) * (eidNClusters_ * eidNFeatures_) +
              static_cast<size_t>(layerClusterCount) * eidNFeatures_;
          
          // feats are clusX,clusY,clusZ,clusE,clusT,clusL
          // taking [0, 1, 5, 3] -> X, Y, abs(L), E
          layerClusterFeatures[base + 0] = static_cast<float>(cl.x());
          layerClusterFeatures[base + 1] = static_cast<float>(cl.y());
          layerClusterFeatures[base + 2] = static_cast<float>(std::abs(layer));
          layerClusterFeatures[base + 3] = static_cast<float>(cl.energy() / static_cast<float>(ts.vertex_multiplicity(k)));

          ++layerClusterCount;
        }

        // clus3d_feat[[0,1,2,4,5]]
        // clus3d_feat = (trkcluseta,trkclusphi,trkclusen,trkclustime, min(clusL),max(clusL))
        // feats : eta, phi, en, minLayer, maxLayer
        tracksterFeatures[static_cast<size_t>(tracksterIdxInBatch)*5 + 0] = ts.barycenter().eta();
        tracksterFeatures[static_cast<size_t>(tracksterIdxInBatch)*5 + 1] = ts.barycenter().phi();
        tracksterFeatures[static_cast<size_t>(tracksterIdxInBatch)*5 + 0] = ts.raw_energy();
        tracksterFeatures[static_cast<size_t>(tracksterIdxInBatch)*5 + 0] = minLayer;
        tracksterFeatures[static_cast<size_t>(tracksterIdxInBatch)*5 + 0] = maxLayer;
      }

      // ---- PID
      ortScratch_.outputs.clear();

      onnxPIDSession_->runIntoVariable(std::make_tuple(
          cms::Ort::ONNXRuntime::TensorArgs{inputNames_[0], layerClusterFeatures, ortScratch_.input_shapes[0]},
          cms::Ort::ONNXRuntime::TensorArgsBool{inputNames_[1], layerClusterMask, ortScratch_.input_shapes[1]},
          cms::Ort::ONNXRuntime::TensorArgs{inputNames_[2], tracksterFeatures, ortScratch_.input_shapes[2]}
        ),
          output_id_, ortScratch_.outputs, {}, nTrackstersInBatch);

      if (!ortScratch_.outputs.empty() && !output_id_.empty()) {
        for (int bi = 0; bi < nTrackstersInBatch; ++bi) {
          auto& ts = tracksters[indices[start + bi]];
          ts.setIdProbability(Trackster::ParticleType::charged_hadron, ortScratch_.outputs[0][bi*2+0]);
          ts.setIdProbability(Trackster::ParticleType::photon, ortScratch_.outputs[0][bi*2+1]);
        }
      }
    }
  }

  void TracksterInferenceByTransformer::fillPSetDescription(edm::ParameterSetDescription& iDesc) {
    TracksterInferenceAlgoBase::fillPSetDescription(iDesc);

    iDesc.add<std::string>("onnxPIDModelPath", "")
        ->setComment("Path to ONNX PID model. If empty, PID inference is skipped.");
    iDesc.add<std::string>("onnxEnergyModelPath", "")
        ->setComment("Path to ONNX energy model. If empty, energy regression is skipped.");

    iDesc.add<std::vector<std::string>>("inputNames", {"layerCluster_features", "layerCluster_mask", "trackster_features"});
    // iDesc.add<std::vector<std::string>>("output_en", {"enreg_output"});
    iDesc.add<std::vector<std::string>>("output_id", {"pid_output"});

    iDesc.add<double>("eid_min_cluster_energy", 1.0);
    // iDesc.add<int>("eid_n_layers", 50);
    iDesc.add<int>("eid_n_clusters", 150);
    iDesc.add<int>("doPID", 1);
    // iDesc.add<int>("doRegression", 0);

    iDesc.addUntracked<int>("miniBatchSize", 256)
        ->setComment("Mini-batch size for inference to limit peak memory usage.");
  }

}  // namespace ticl
