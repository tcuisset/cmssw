#ifndef SimGeneral_CaloAnalysis_SimClusterBuilder_h
#define SimGeneral_CaloAnalysis_SimClusterBuilder_h

#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Framework/interface/ProducesCollector.h"

#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"

#include "SimGeneral/MixingModule/interface/DecayGraph.h"

#include <string>


namespace {
  struct SimClusterConfig {
    SimClusterCollection outputClusters;
    edm::EDPutTokenT<SimClusterCollection> outputClusters_token;

    // For the map back to "CaloParticle" it can be either to CaloParticle dataformat or "CaloParticle" SimCluster dataformat.
    // We choose SimCluster dataformat
    SimClusterRefVector clustersToCaloParticleMap;
    edm::EDPutTokenT<SimClusterRefVector> clustersToCaloParticleMap_token;

    SimClusterConfig(edm::ProducesCollector& c, std::string tag) 
      : outputClusters_token(c.produces<SimClusterCollection>(tag)), clustersToCaloParticleMap_token(c.produces<SimClusterRefVector>(tag))
    {}
    void clear() {
      outputClusters.clear();
      clustersToCaloParticleMap.clear();
    }
  };

    class SubClusterMergerBase {
  public:
    void start_subcluster(DecayChain::edge_descriptor e, const DecayChain &g, std::size_t currentSubClusterIndex);
    void end_parentCluster(std::span<SimCluster> subClustersBuilt, std::size_t currentCaloParticleIndex);
  };


}  // namespace

#endif