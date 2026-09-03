#ifndef RecoHGCal_TICL_TracksterPIDFeatures_h
#define RecoHGCal_TICL_TracksterPIDFeatures_h

#include <cstddef>
#include <vector>

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"

namespace ticlgeom {
  class Tools;
}

namespace ticl {

  // Shared feature definition used by both CMSSW inference and training ntupling.
  // Keeping the transformations here prevents train/inference feature drift.
  class TracksterPIDFeatures {
  public:
    static constexpr int summaryFeatureCount = 16;
    static constexpr int layerFeatureCount = 8;
    static constexpr int clusterFeatureCount = 8;
    static constexpr int imageChannelCount = 3;

    TracksterPIDFeatures(int nLayers, int maxClusters) : nLayers_(nLayers), maxClusters_(maxClusters) {}

    int nLayers() const { return nLayers_; }
    int maxClusters() const { return maxClusters_; }
    std::size_t summarySize() const { return summaryFeatureCount; }
    std::size_t layerSize() const { return static_cast<std::size_t>(nLayers_) * layerFeatureCount; }
    std::size_t clusterSize() const { return static_cast<std::size_t>(maxClusters_) * clusterFeatureCount; }
    std::size_t imageSize() const {
      return static_cast<std::size_t>(imageChannelCount) * nLayers_ * maxClusters_;
    }

    void fillSummary(float* output,
                     Trackster const& trackster,
                     std::vector<reco::CaloCluster> const& layerClusters,
                     ticlgeom::Tools const& rhtools) const;
    void fillLayer(float* output,
                   Trackster const& trackster,
                   std::vector<reco::CaloCluster> const& layerClusters,
                   ticlgeom::Tools const& rhtools) const;
    void fillCluster(float* output,
                     Trackster const& trackster,
                     std::vector<reco::CaloCluster> const& layerClusters,
                     ticlgeom::Tools const& rhtools,
                     std::vector<int>& order) const;
    void fillImage(float* output,
                   Trackster const& trackster,
                   std::vector<reco::CaloCluster> const& layerClusters,
                   ticlgeom::Tools const& rhtools,
                   std::vector<int>& order,
                   std::vector<int>& seenClusters) const;

  private:
    int nLayers_;
    int maxClusters_;
  };

}  // namespace ticl

#endif
