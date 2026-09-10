// Author: Felice Pantaleo - felice.pantaleo@cern.ch
// Date: 07/2024

#ifndef RecoHGCal_TICL_TracksterInferenceAlgoFactory_h
#define RecoHGCal_TICL_TracksterInferenceAlgoFactory_h

#include <memory>
#include <string>

#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "RecoHGCal/TICL/interface/TracksterInferenceAlgoBase.h"

typedef edmplugin::PluginFactory<ticl::TracksterInferenceAlgoBase*(const edm::ParameterSet&,
                                                                   ticl::TICLONNXGlobalCache const*)>
    TracksterInferenceAlgoFactory;

namespace ticl {
  std::unique_ptr<TracksterInferenceAlgoBase> makeTracksterInferenceAlgo(edm::ParameterSet const& modulePSet,
                                                                         TICLONNXGlobalCache const* cache);
}  // namespace ticl

#endif
