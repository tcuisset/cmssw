#include "RecoHGCal/TICL/interface/TracksterInferenceByBenchmark.h"

#include <algorithm>
#include <chrono>

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

namespace {
  ticl::TracksterInferenceByBenchmark::InputFormat parseInputFormat(std::string const& value) {
    using InputFormat = ticl::TracksterInferenceByBenchmark::InputFormat;
    if (value == "summary")
      return InputFormat::summary;
    if (value == "layer")
      return InputFormat::layer;
    if (value == "cluster")
      return InputFormat::cluster;
    if (value == "image")
      return InputFormat::image;
    throw cms::Exception("Configuration")
        << "Unknown PID benchmark inputFormat '" << value << "'; expected summary, layer, cluster, or image";
  }
}  // namespace

namespace ticl {

  TracksterInferenceByBenchmark::TracksterInferenceByBenchmark(const edm::ParameterSet& conf,
                                                               TICLONNXGlobalCache const* cache)
      : TracksterInferenceAlgoBase(conf, cache),
        inputNames_(conf.getParameter<std::vector<std::string>>("inputNames")),
        outputNames_(conf.getParameter<std::vector<std::string>>("outputNames")),
        minTracksterEnergy_(conf.getParameter<double>("minTracksterEnergy")),
        nLayers_(conf.getParameter<int>("nLayers")),
        maxClusters_(conf.getParameter<int>("maxClusters")),
        miniBatchSize_(conf.getUntrackedParameter<int>("miniBatchSize", 64)),
        reportTiming_(conf.getUntrackedParameter<bool>("reportTiming", false)),
        inputFormat_(parseInputFormat(conf.getParameter<std::string>("inputFormat"))),
        features_(nLayers_, maxClusters_) {
    if (nLayers_ <= 0 || maxClusters_ <= 0 || miniBatchSize_ <= 0) {
      throw cms::Exception("Configuration") << "nLayers, maxClusters, and miniBatchSize must all be positive";
    }

    const auto modelPath = conf.getParameter<std::string>("onnxModelPath");
    if (cache_ != nullptr && !modelPath.empty()) {
      onnxSession_ = cache_->getByModelPathString(modelPath);
    }
    enabled_ = (onnxSession_ != nullptr);
  }

  std::size_t TracksterInferenceByBenchmark::floatsPerTrackster() const {
    switch (inputFormat_) {
      case InputFormat::summary:
        return features_.summarySize();
      case InputFormat::layer:
        return features_.layerSize();
      case InputFormat::cluster:
        return features_.clusterSize();
      case InputFormat::image:
        return features_.imageSize();
    }
    return 0;
  }

  std::vector<int64_t> TracksterInferenceByBenchmark::inputShape(int batchSize) const {
    switch (inputFormat_) {
      case InputFormat::summary:
        return {batchSize, TracksterPIDFeatures::summaryFeatureCount};
      case InputFormat::layer:
        return {batchSize, nLayers_, TracksterPIDFeatures::layerFeatureCount};
      case InputFormat::cluster:
        return {batchSize, maxClusters_, TracksterPIDFeatures::clusterFeatureCount};
      case InputFormat::image:
        return {batchSize, TracksterPIDFeatures::imageChannelCount, nLayers_, maxClusters_};
    }
    return {};
  }

  void TracksterInferenceByBenchmark::runInference(const std::vector<reco::CaloCluster>& layerClusters,
                                                   std::vector<Trackster>& tracksters,
                                                   const ticlgeom::Tools& rhtools) const {
    if (!enabled_ || tracksters.empty())
      return;

    const auto timingStart = std::chrono::steady_clock::now();

    std::vector<int> indices;
    indices.reserve(tracksters.size());
    for (int index = 0, size = tracksters.size(); index < size; ++index) {
      float energy = 0.f;
      for (unsigned int vertex : tracksters[index].vertices()) {
        energy += layerClusters[vertex].energy();
        if (energy >= minTracksterEnergy_) {
          indices.push_back(index);
          break;
        }
      }
    }
    if (indices.empty())
      return;

    OrtScratch scratch;
    scratch.inputs.resize(1);
    scratch.input_shapes.resize(1);
    std::vector<int> order;
    std::vector<int> seenClusters(nLayers_);
    const std::size_t stride = floatsPerTrackster();

    for (int start = 0, total = indices.size(); start < total; start += miniBatchSize_) {
      const int batchSize = std::min(miniBatchSize_, total - start);
      scratch.input_shapes[0] = inputShape(batchSize);
      scratch.inputs[0].assign(static_cast<std::size_t>(batchSize) * stride, 0.f);

      for (int batchIndex = 0; batchIndex < batchSize; ++batchIndex) {
        auto const& ts = tracksters[indices[start + batchIndex]];
        float* output = scratch.inputs[0].data() + static_cast<std::size_t>(batchIndex) * stride;
        switch (inputFormat_) {
          case InputFormat::summary:
            features_.fillSummary(output, ts, layerClusters, rhtools);
            break;
          case InputFormat::layer:
            features_.fillLayer(output, ts, layerClusters, rhtools);
            break;
          case InputFormat::cluster:
            features_.fillCluster(output, ts, layerClusters, rhtools, order);
            break;
          case InputFormat::image:
            features_.fillImage(output, ts, layerClusters, rhtools, order, seenClusters);
            break;
        }
      }

      scratch.outputs.clear();
      onnxSession_->runInto(
          inputNames_, scratch.inputs, scratch.input_shapes, outputNames_, scratch.outputs, {}, batchSize);
      if (scratch.outputs.empty() || scratch.outputs[0].size() < static_cast<std::size_t>(2 * batchSize)) {
        throw cms::Exception("InvalidONNXOutput") << "Binary PID model must return at least two floats per trackster";
      }

      auto const& probabilities = scratch.outputs[0];
      for (int batchIndex = 0; batchIndex < batchSize; ++batchIndex) {
        auto& ts = tracksters[indices[start + batchIndex]];
        ts.zeroProbabilities();
        ts.setIdProbability(Trackster::ParticleType::electron, probabilities[2 * batchIndex]);
        ts.setIdProbability(Trackster::ParticleType::charged_hadron, probabilities[2 * batchIndex + 1]);
      }
    }

    if (reportTiming_) {
      const auto elapsed = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - timingStart);
      edm::LogVerbatim("PIDBenchmarkTiming") << "milliseconds=" << elapsed.count() << " eligible=" << indices.size();
    }
  }

  void TracksterInferenceByBenchmark::fillPSetDescription(edm::ParameterSetDescription& desc) {
    TracksterInferenceAlgoBase::fillPSetDescription(desc);
    desc.add<std::string>("onnxModelPath", "");
    desc.add<std::vector<std::string>>("inputNames", {"input"});
    desc.add<std::vector<std::string>>("outputNames", {"pid_output"});
    desc.add<std::string>("inputFormat", "layer");
    desc.add<double>("minTracksterEnergy", 1.0);
    desc.add<int>("nLayers", 50);
    desc.add<int>("maxClusters", 64);
    desc.addUntracked<int>("miniBatchSize", 64);
    desc.addUntracked<bool>("reportTiming", false);
  }

}  // namespace ticl
