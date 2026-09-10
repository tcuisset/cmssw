#include "RecoHGCal/TICL/interface/TracksterInferenceByONNX.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/TICLGeomTools.h"

namespace {
  ticl::TracksterInferenceByONNX::InputFormat parseInputFormat(std::string const& value) {
    using InputFormat = ticl::TracksterInferenceByONNX::InputFormat;
    if (value == "legacyCNN")
      return InputFormat::legacyCNN;
    if (value == "legacyPFN")
      return InputFormat::legacyPFN;
    if (value == "summary")
      return InputFormat::summary;
    if (value == "cluster")
      return InputFormat::cluster;
    throw cms::Exception("Configuration")
        << "Unknown Trackster PID inputFormat '" << value << "'; expected legacyCNN, legacyPFN, summary, or cluster";
  }
}  // namespace

namespace ticl {

  TracksterInferenceByONNX::TracksterInferenceByONNX(edm::ParameterSet const& conf, TICLONNXGlobalCache const* cache)
      : TracksterInferenceAlgoBase(conf, cache),
        inputNames_(conf.getParameter<std::vector<std::string>>("inputNames")),
        outputNamesPID_(conf.getParameter<std::vector<std::string>>("outputNamesPID")),
        outputNamesEnergy_(conf.getParameter<std::vector<std::string>>("outputNamesEnergy")),
        outputProbabilityIndices_(conf.getParameter<std::vector<unsigned int>>("outputProbabilityIndices")),
        minTracksterEnergy_(conf.getParameter<double>("minTracksterEnergy")),
        nLayers_(conf.getParameter<int>("nLayers")),
        maxClusters_(conf.getParameter<int>("maxClusters")),
        miniBatchSize_(conf.getUntrackedParameter<int>("miniBatchSize", 64)),
        doPID_(conf.getParameter<bool>("doPID")),
        doRegression_(conf.getParameter<bool>("doRegression")),
        reportTiming_(conf.getUntrackedParameter<bool>("reportTiming", false)),
        inputFormat_(parseInputFormat(conf.getParameter<std::string>("inputFormat"))),
        features_(nLayers_, maxClusters_) {
    if (minTracksterEnergy_ < 0. || nLayers_ <= 0 || maxClusters_ <= 0 || miniBatchSize_ <= 0) {
      throw cms::Exception("Configuration")
          << "minTracksterEnergy must be non-negative and nLayers, maxClusters, and miniBatchSize must be positive";
    }
    if (doPID_ && (outputNamesPID_.empty() || outputProbabilityIndices_.empty())) {
      throw cms::Exception("Configuration")
          << "PID inference requires non-empty outputNamesPID and outputProbabilityIndices";
    }
    if (doRegression_ && outputNamesEnergy_.empty()) {
      throw cms::Exception("Configuration") << "Energy regression requires a non-empty outputNamesEnergy";
    }
    if (inputFormat_ == InputFormat::legacyPFN && inputNames_.size() != 2) {
      throw cms::Exception("Configuration") << "legacyPFN inputFormat requires exactly two inputNames";
    }
    if (inputFormat_ != InputFormat::legacyPFN && inputNames_.size() != 1) {
      throw cms::Exception("Configuration") << "The selected inputFormat requires exactly one inputNames entry";
    }

    auto const pidModel = conf.getParameter<std::string>("onnxPIDModelPath");
    auto const energyModel = conf.getParameter<std::string>("onnxEnergyModelPath");
    if (cache_ != nullptr) {
      if (doPID_ && !pidModel.empty())
        onnxPIDSession_ = cache_->getByModelPathString(pidModel);
      if (doRegression_ && !energyModel.empty())
        onnxEnergySession_ = cache_->getByModelPathString(energyModel);
    }
  }

  bool TracksterInferenceByONNX::select(Trackster& trackster,
                                        std::vector<reco::CaloCluster> const& layerClusters,
                                        ticlgeom::Tools const& rhtools) const {
    if (inputFormat_ == InputFormat::legacyPFN) {
      for (auto vertex : trackster.vertices()) {
        if (rhtools.isBarrel(layerClusters[vertex].seed()))
          return false;
      }
      trackster.setRegressedEnergy(0.f);
      trackster.zeroProbabilities();
      return true;
    }

    float energy = 0.f;
    for (auto vertex : trackster.vertices()) {
      energy += static_cast<float>(layerClusters[vertex].energy());
      if (energy >= minTracksterEnergy_) {
        trackster.zeroProbabilities();
        return true;
      }
    }
    return false;
  }

  void TracksterInferenceByONNX::fillInputs(OrtScratch& scratch,
                                            int batchSize,
                                            int start,
                                            std::vector<int> const& indices,
                                            std::vector<Trackster> const& tracksters,
                                            std::vector<reco::CaloCluster> const& layerClusters,
                                            ticlgeom::Tools const& rhtools,
                                            std::vector<int>& order,
                                            std::vector<int>& seenClusters) const {
    if (inputFormat_ == InputFormat::summary || inputFormat_ == InputFormat::cluster) {
      const std::size_t stride =
          inputFormat_ == InputFormat::summary ? features_.summarySize() : features_.clusterSize();
      scratch.input_shapes[0] =
          inputFormat_ == InputFormat::summary
              ? std::vector<int64_t>{batchSize, TracksterPIDFeatures::summaryFeatureCount}
              : std::vector<int64_t>{batchSize, maxClusters_, TracksterPIDFeatures::clusterFeatureCount};
      scratch.inputs[0].assign(static_cast<std::size_t>(batchSize) * stride, 0.f);
      for (int batchIndex = 0; batchIndex < batchSize; ++batchIndex) {
        auto const& trackster = tracksters[indices[start + batchIndex]];
        float* output = scratch.inputs[0].data() + static_cast<std::size_t>(batchIndex) * stride;
        if (inputFormat_ == InputFormat::summary)
          features_.fillSummary(output, trackster, layerClusters, rhtools);
        else
          features_.fillCluster(output, trackster, layerClusters, rhtools, order);
      }
      return;
    }

    const int featureCount = inputFormat_ == InputFormat::legacyCNN ? 3 : 7;
    scratch.input_shapes[0] = {batchSize, nLayers_, maxClusters_, featureCount};
    scratch.inputs[0].assign(static_cast<std::size_t>(batchSize) * nLayers_ * maxClusters_ * featureCount, 0.f);
    if (inputFormat_ == InputFormat::legacyPFN) {
      scratch.input_shapes[1] = {batchSize, featureCount};
      scratch.inputs[1].resize(static_cast<std::size_t>(batchSize) * featureCount);
    }

    for (int batchIndex = 0; batchIndex < batchSize; ++batchIndex) {
      auto const& trackster = tracksters[indices[start + batchIndex]];
      if (inputFormat_ == InputFormat::legacyPFN) {
        float* tracksterFeatures = scratch.inputs[1].data() + static_cast<std::size_t>(batchIndex) * featureCount;
        tracksterFeatures[0] = trackster.raw_energy();
        tracksterFeatures[1] = trackster.raw_em_energy();
        tracksterFeatures[2] = trackster.barycenter().x();
        tracksterFeatures[3] = trackster.barycenter().y();
        tracksterFeatures[4] = std::abs(trackster.barycenter().z());
        tracksterFeatures[5] = std::abs(trackster.barycenter().eta());
        tracksterFeatures[6] = trackster.barycenter().phi();
      }

      order.resize(trackster.vertices().size());
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), [&layerClusters, &trackster](int a, int b) {
        return layerClusters[trackster.vertices(a)].energy() > layerClusters[trackster.vertices(b)].energy();
      });
      std::fill(seenClusters.begin(), seenClusters.end(), 0);
      for (int vertexPosition : order) {
        auto const& cluster = layerClusters[trackster.vertices(vertexPosition)];
        const int layer = rhtools.getLayerWithOffset(cluster.hitsAndFractions()[0].first) - 1;
        if (layer < 0 || layer >= nLayers_ || seenClusters[layer] >= maxClusters_)
          continue;
        const std::size_t offset =
            ((static_cast<std::size_t>(batchIndex) * nLayers_ + layer) * maxClusters_ + seenClusters[layer]++) *
            featureCount;
        float* clusterFeatures = scratch.inputs[0].data() + offset;
        clusterFeatures[0] = cluster.energy() / trackster.vertex_multiplicity(vertexPosition);
        clusterFeatures[1] = std::abs(cluster.eta());
        clusterFeatures[2] = cluster.phi();
        if (inputFormat_ == InputFormat::legacyPFN) {
          clusterFeatures[3] = cluster.x();
          clusterFeatures[4] = cluster.y();
          clusterFeatures[5] = std::abs(cluster.z());
          clusterFeatures[6] = cluster.hitsAndFractions().size();
        }
      }
    }
  }

  void TracksterInferenceByONNX::assignPID(std::vector<float> const& output,
                                           int batchSize,
                                           int start,
                                           std::vector<int> const& indices,
                                           std::vector<Trackster>& tracksters) const {
    const std::size_t outputWidth = outputProbabilityIndices_.size();
    if (output.size() < static_cast<std::size_t>(batchSize) * outputWidth) {
      throw cms::Exception("InvalidONNXOutput") << "PID model returned " << output.size() << " values for " << batchSize
                                                << " tracksters and " << outputWidth << " configured probability slots";
    }
    for (int batchIndex = 0; batchIndex < batchSize; ++batchIndex) {
      auto& trackster = tracksters[indices[start + batchIndex]];
      trackster.zeroProbabilities();
      for (std::size_t outputIndex = 0; outputIndex < outputWidth; ++outputIndex) {
        const auto probabilityIndex = outputProbabilityIndices_[outputIndex];
        if (probabilityIndex >= trackster.id_probabilities().size()) {
          throw cms::Exception("Configuration")
              << "outputProbabilityIndices contains out-of-range index " << probabilityIndex;
        }
        trackster.setIdProbability(static_cast<Trackster::ParticleType>(probabilityIndex),
                                   output[static_cast<std::size_t>(batchIndex) * outputWidth + outputIndex]);
      }
    }
  }

  void TracksterInferenceByONNX::runInference(std::vector<reco::CaloCluster> const& layerClusters,
                                              std::vector<Trackster>& tracksters,
                                              ticlgeom::Tools const& rhtools) const {
    if ((!onnxPIDSession_ && !onnxEnergySession_) || tracksters.empty())
      return;
    auto const timingStart = std::chrono::steady_clock::now();

    std::vector<int> indices;
    indices.reserve(tracksters.size());
    for (int index = 0, size = tracksters.size(); index < size; ++index) {
      if (select(tracksters[index], layerClusters, rhtools))
        indices.push_back(index);
    }
    if (indices.empty())
      return;

    OrtScratch scratch;
    const int inputCount = inputFormat_ == InputFormat::legacyPFN ? 2 : 1;
    scratch.inputs.resize(inputCount);
    scratch.input_shapes.resize(inputCount);
    std::vector<int> order;
    std::vector<int> seenClusters(nLayers_);

    for (int start = 0, total = indices.size(); start < total; start += miniBatchSize_) {
      const int batchSize = std::min(miniBatchSize_, total - start);
      fillInputs(scratch, batchSize, start, indices, tracksters, layerClusters, rhtools, order, seenClusters);

      if (doRegression_ && onnxEnergySession_) {
        scratch.outputs.clear();
        onnxEnergySession_->runInto(
            inputNames_, scratch.inputs, scratch.input_shapes, outputNamesEnergy_, scratch.outputs, {}, batchSize);
        if (scratch.outputs.empty() || scratch.outputs[0].size() < static_cast<std::size_t>(batchSize))
          throw cms::Exception("InvalidONNXOutput") << "Energy model returned too few values";
        for (int batchIndex = 0; batchIndex < batchSize; ++batchIndex) {
          auto& trackster = tracksters[indices[start + batchIndex]];
          const float regressedEnergy = scratch.outputs[0][batchIndex];
          trackster.setRegressedEnergy(trackster.raw_energy() > minTracksterEnergy_ ? regressedEnergy
                                                                                    : trackster.raw_energy());
        }
      }
      if (doPID_ && onnxPIDSession_) {
        scratch.outputs.clear();
        onnxPIDSession_->runInto(
            inputNames_, scratch.inputs, scratch.input_shapes, outputNamesPID_, scratch.outputs, {}, batchSize);
        if (scratch.outputs.empty())
          throw cms::Exception("InvalidONNXOutput") << "PID model returned no output tensors";
        assignPID(scratch.outputs[0], batchSize, start, indices, tracksters);
      }
    }

    if (reportTiming_) {
      auto const elapsed = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - timingStart);
      edm::LogVerbatim("TracksterInferenceTiming")
          << "milliseconds=" << elapsed.count() << " eligible=" << indices.size();
    }
  }

  void TracksterInferenceByONNX::fillPSetDescription(edm::ParameterSetDescription& desc) {
    TracksterInferenceAlgoBase::fillPSetDescription(desc);
    desc.add<std::string>("onnxPIDModelPath", "")
        ->setComment("FileInPath-relative ONNX PID model path; empty disables PID inference");
    desc.add<std::string>("onnxEnergyModelPath", "")
        ->setComment("FileInPath-relative ONNX regression model path; empty disables energy regression");
    desc.add<std::vector<std::string>>("inputNames", {"input"});
    desc.add<std::vector<std::string>>("outputNamesPID", {"pid_output"});
    desc.add<std::vector<std::string>>("outputNamesEnergy", {"enreg_output"});
    desc.add<std::vector<unsigned int>>("outputProbabilityIndices", {0, 1, 2, 3, 4, 5, 6, 7})
        ->setComment("Trackster ParticleType slots corresponding to consecutive PID model outputs");
    desc.add<std::string>("inputFormat", "legacyCNN")
        ->setComment("Input tensor layout: legacyCNN, legacyPFN, summary, or cluster");
    desc.add<double>("minTracksterEnergy", 1.0);
    desc.add<int>("nLayers", 50);
    desc.add<int>("maxClusters", 10);
    desc.add<bool>("doPID", true);
    desc.add<bool>("doRegression", false);
    desc.addUntracked<int>("miniBatchSize", 64)
        ->setComment("Maximum inference batch size, bounding temporary tensor memory");
    desc.addUntracked<bool>("reportTiming", false);
  }

}  // namespace ticl
