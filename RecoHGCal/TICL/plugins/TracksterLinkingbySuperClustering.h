#ifndef RecoHGCal_TICL_TracksterLinkingSuperClustering_H
#define RecoHGCal_TICL_TracksterLinkingSuperClustering_H

#include "RecoHGCal/TICL/interface/TracksterLinkingAlgoBase.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

namespace ticl {

  class TracksterLinkingbySuperClustering : public TracksterLinkingAlgoBase {
  public:
    TracksterLinkingbySuperClustering(const edm::ParameterSet& ps,
                                      edm::ConsumesCollector iC,
                                      cms::Ort::ONNXRuntime const* onnxRuntime = nullptr);
    /* virtual */ ~TracksterLinkingbySuperClustering() override {}
    static void fillPSetDescription(edm::ParameterSetDescription& iDesc);

    void linkTracksters(const Inputs& input,
                        std::vector<Trackster>& resultTracksters,
                        std::vector<std::vector<unsigned int>>& linkedResultTracksters,
                        std::vector<std::vector<unsigned int>>& linkedTracksterIdToInputTracksterId) override;
    void initialize(const HGCalDDDConstants* hgcons,
                    const hgcal::RecHitTools rhtools,
                    const edm::ESHandle<MagneticField> bfieldH,
                    const edm::ESHandle<Propagator> propH) override;
    
  private:
    const std::string dnnVersion_;  // Version identifier of the DNN (to choose which inputs to use)
    double
        nnWorkingPoint_;  // Working point for neural network (above this score, consider the trackster candidate for superclustering)
    float deltaEtaWindow_;            // Delta eta window to consider trackster seed-candidate pairs for inference
    float deltaPhiWindow_;            // Delta phi window
    float seedPtThreshold_;           // Min pT for a trackster to be considered as supercluster seed
    float candidateEnergyThreshold_;  // Min energy for a trackster to be superclustered as candidate
  };

}  // namespace ticl

#endif
