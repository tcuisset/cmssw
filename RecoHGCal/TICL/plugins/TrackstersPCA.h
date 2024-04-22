#ifndef RECOHGCAL_TICL_TRACKSTERSPCA_H
#define RECOHGCAL_TICL_TRACKSTERSPCA_H

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include <vector>
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

namespace ticl {
  void assignPCAtoTracksters(std::vector<Trackster> &,
                             const std::vector<reco::CaloCluster> &,
                             const edm::ValueMap<std::pair<float, float>> &,
                             double,
                             bool computeLocalTime = false,
                             bool energyWeight = true,
                             const hgcal::RecHitTools = hgcal::RecHitTools(),
                             int minLayer = 10,
                             int maxLayer = 10,
                             bool clean = false);
  std::pair<float, float> computeLocalTracksterTime(const Trackster &trackster,
                                                    const std::vector<reco::CaloCluster> &layerClusters,
                                                    const edm::ValueMap<std::pair<float, float>> &layerClustersTime,
                                                    const Eigen::Vector3f &barycenter,
                                                    size_t N);
  std::pair<float, float> computeTracksterTime(const Trackster &trackster,
                                               const edm::ValueMap<std::pair<float, float>> &layerClustersTime,
                                               size_t N);

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
      const auto &thisLC = layerClusters[ts.vertices(i)];
      auto layer = getLayerFromLC(thisLC, rhtools);
      result[layer].push_back(i);
    }
    return result;
  }
}  // namespace ticl
#endif
