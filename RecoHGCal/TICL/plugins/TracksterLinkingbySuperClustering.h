#ifndef RecoHGCal_TICL_TracksterLinkingSuperClustering_H
#define RecoHGCal_TICL_TracksterLinkingSuperClustering_H

#include "RecoHGCal/TICL/interface/TracksterLinkingAlgoBase.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"


namespace ticl {

class TracksterLinkingbySuperClustering : public TracksterLinkingAlgoBase {
public:
  TracksterLinkingbySuperClustering(const edm::ParameterSet& ps, edm::ConsumesCollector iC)
      : TracksterLinkingAlgoBase(ps, iC),
      onnxRuntime_(ps.getParameter<edm::FileInPath>("dnn_model_path").fullPath()),
      dnnVersion_(ps.getParameter<std::string>("dnnVersion")),
      nnWorkingPoint_(ps.getParameter<double>("nnWorkingPoint")),
      deltaEtaWindow_(ps.getParameter<double>("deltaEtaWindow")),
      deltaPhiWindow_(ps.getParameter<double>("deltaPhiWindow")),
      seedPtThreshold_(ps.getParameter<double>("seedPtThreshold")),
      candidateEnergyThreshold_(ps.getParameter<double>("candidateEnergyThreshold"))
  {}
  virtual ~TracksterLinkingbySuperClustering() {}
  static void fillPSetDescription(edm::ParameterSetDescription& iDesc);

  void linkTracksters(const Inputs& input, std::vector<Trackster>& resultTracksters,
                      std::vector<std::vector<unsigned int>>& linkedResultTracksters,
                      std::vector<std::vector<unsigned int>>& linkedTracksterIdToInputTracksterId) override;
  
private:
  cms::Ort::ONNXRuntime onnxRuntime_; // TODO this only needs to be instantiated globally, but currently it 

  const std::string dnnVersion_; // Version identifier of the DNN (to choose which inputs to use)
  double nnWorkingPoint_; // Working point for neural network (above this score, consider the trackster candidate for superclustering)
  float deltaEtaWindow_; // Delta eta window to consider trackster seed-candidate pairs for inference
  float deltaPhiWindow_; // Delta phi window
  float seedPtThreshold_; // Min pT for a trackster to be considered as supercluster seed
  float candidateEnergyThreshold_; // Min energy for a trackster to be superclustered as candidate
};

}  // namespace ticl

#endif
