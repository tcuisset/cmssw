#ifndef RecoHGCal_TICL_TracksterInferenceByONNX_h
#define RecoHGCal_TICL_TracksterInferenceByONNX_h

#include <cstddef>
#include <string>
#include <vector>

#include "RecoHGCal/TICL/interface/TracksterInferenceAlgoBase.h"
#include "RecoHGCal/TICL/interface/TracksterPIDFeatures.h"

namespace ticl {

  // Configurable ONNX inference for both the established TICL models and the
  // compact feature families used by newer PID models.
  class TracksterInferenceByONNX final : public TracksterInferenceAlgoBase {
  public:
    enum class InputFormat { legacyCNN, legacyPFN, summary, cluster };

    TracksterInferenceByONNX(edm::ParameterSet const& conf, TICLONNXGlobalCache const* cache);

    void runInference(std::vector<reco::CaloCluster> const& layerClusters,
                      std::vector<Trackster>& tracksters,
                      ticlgeom::Tools const& rhtools) const override;

    static void fillPSetDescription(edm::ParameterSetDescription& desc);

  private:
    bool select(Trackster& trackster,
                std::vector<reco::CaloCluster> const& layerClusters,
                ticlgeom::Tools const& rhtools) const;
    void fillInputs(OrtScratch& scratch,
                    int batchSize,
                    int start,
                    std::vector<int> const& indices,
                    std::vector<Trackster> const& tracksters,
                    std::vector<reco::CaloCluster> const& layerClusters,
                    ticlgeom::Tools const& rhtools,
                    std::vector<int>& order,
                    std::vector<int>& seenClusters) const;
    void assignPID(std::vector<float> const& output,
                   int batchSize,
                   int start,
                   std::vector<int> const& indices,
                   std::vector<Trackster>& tracksters) const;

    std::vector<std::string> const inputNames_;
    std::vector<std::string> const outputNamesPID_;
    std::vector<std::string> const outputNamesEnergy_;
    std::vector<unsigned int> const outputProbabilityIndices_;
    double const minTracksterEnergy_;
    int const nLayers_;
    int const maxClusters_;
    int const miniBatchSize_;
    bool const doPID_;
    bool const doRegression_;
    bool const reportTiming_;
    InputFormat const inputFormat_;
    TracksterPIDFeatures const features_;

    cms::Ort::ONNXRuntime const* onnxPIDSession_ = nullptr;
    cms::Ort::ONNXRuntime const* onnxEnergySession_ = nullptr;
  };

}  // namespace ticl

#endif
