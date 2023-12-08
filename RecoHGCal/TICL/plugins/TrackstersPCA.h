#ifndef RECOHGCAL_TICL_TRACKSTERSPCA_H
#define RECOHGCAL_TICL_TRACKSTERSPCA_H

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include <vector>
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

typedef math::XYZVector Vector;

namespace ticl {
  void assignPCAtoTracksters(std::vector<Trackster> &,
			     const std::vector<reco::CaloCluster> &,
			     const edm::ValueMap<std::pair<float, float>> &,
			     double,
			     bool energyWeight = true,
			     const hgcal::RecHitTools = hgcal::RecHitTools(), 
			     int minLayer = 10,
			     int maxLayer = 10,
			     bool clean = false);
  
  inline unsigned getLayerFromLC(const reco::CaloCluster LC, hgcal::RecHitTools rhtools) {
    std::vector<std::pair<DetId, float>> thisclusterHits = LC.hitsAndFractions();
    auto layer = rhtools.getLayerWithOffset(thisclusterHits[0].first);
    return layer;
  }
  
  inline std::vector<std::vector<unsigned>> sortByLayer(const Trackster &ts,
							const std::vector<reco::CaloCluster> &layerClusters,
							hgcal::RecHitTools rhtools) {
    size_t N = ts.vertices().size();
    
    std::vector<std::vector<unsigned>> result;
    result.resize(rhtools.lastLayer() + 1);
    
    for (unsigned i = 0; i < N; ++i) {
      auto thisLC = layerClusters[ts.vertices(i)];
      auto layer = getLayerFromLC(thisLC, rhtools);
      result[layer].push_back(i);
    }
    return result;
  }
  
}
#endif



