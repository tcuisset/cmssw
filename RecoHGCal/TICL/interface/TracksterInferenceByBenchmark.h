#ifndef RecoHGCal_TICL_TracksterInferenceByBenchmark_h
#define RecoHGCal_TICL_TracksterInferenceByBenchmark_h

#include <string>
#include <vector>

#include "RecoHGCal/TICL/interface/TICLONNXGlobalCache.h"
#include "RecoHGCal/TICL/interface/TracksterInferenceAlgoBase.h"
#include "RecoHGCal/TICL/interface/TracksterPIDFeatures.h"

namespace ticl {

  class TracksterInferenceByBenchmark final : public TracksterInferenceAlgoBase {
  public:
    enum class InputFormat { summary, layer, cluster, image };

    TracksterInferenceByBenchmark(const edm::ParameterSet& conf, TICLONNXGlobalCache const* cache);

    void runInference(const std::vector<reco::CaloCluster>& layerClusters,
                      std::vector<Trackster>& tracksters,
                      const ticlgeom::Tools& rhtools) const override;

    static void fillPSetDescription(edm::ParameterSetDescription& iDesc);

  private:
    std::vector<int64_t> inputShape(int batchSize) const;
    std::size_t floatsPerTrackster() const;

    const std::vector<std::string> inputNames_;
    const std::vector<std::string> outputNames_;
    const double minTracksterEnergy_;
    const int nLayers_;
    const int maxClusters_;
    const int miniBatchSize_;
    const bool reportTiming_;
    const InputFormat inputFormat_;
    const TracksterPIDFeatures features_;

    cms::Ort::ONNXRuntime const* onnxSession_ = nullptr;
    bool enabled_ = false;

  };

}  // namespace ticl

#endif
