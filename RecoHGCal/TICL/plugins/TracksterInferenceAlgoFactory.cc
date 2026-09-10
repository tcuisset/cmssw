#include "RecoHGCal/TICL/interface/TracksterInferenceAlgoFactory.h"
#include "RecoHGCal/TICL/interface/TracksterInferenceByONNX.h"

#include "FWCore/ParameterSet/interface/ValidatedPluginFactoryMacros.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginMacros.h"

EDM_REGISTER_VALIDATED_PLUGINFACTORY(TracksterInferenceAlgoFactory, "TracksterInferenceAlgoFactory");
DEFINE_EDM_VALIDATED_PLUGIN(TracksterInferenceAlgoFactory, ticl::TracksterInferenceByONNX, "TracksterInferenceByONNX");

std::unique_ptr<ticl::TracksterInferenceAlgoBase> ticl::makeTracksterInferenceAlgo(edm::ParameterSet const& modulePSet,
                                                                                   TICLONNXGlobalCache const* cache) {
  auto const plugin = modulePSet.getParameter<std::string>("inferenceAlgo");
  if (plugin.empty())
    return nullptr;

  auto const pluginPSet = modulePSet.getParameter<edm::ParameterSet>("pluginInferenceAlgo" + plugin);
  auto const hasModel = [&pluginPSet](char const* parameter) {
    return pluginPSet.existsAs<std::string>(parameter, true) &&
           !pluginPSet.getParameter<std::string>(parameter).empty();
  };
  if (!hasModel("onnxPIDModelPath") && !hasModel("onnxEnergyModelPath"))
    return nullptr;

  return std::unique_ptr<TracksterInferenceAlgoBase>(
      TracksterInferenceAlgoFactory::get()->create(plugin, pluginPSet, cache));
}
