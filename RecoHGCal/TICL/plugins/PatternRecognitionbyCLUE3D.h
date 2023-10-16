// Author: Marco Rovere - marco.rovere@cern.ch
// Date: 04/2021

#ifndef __RecoHGCal_TICL_PRbyCLUE3D_H__
#define __RecoHGCal_TICL_PRbyCLUE3D_H__
#include <memory>  // unique_ptr
#include "RecoHGCal/TICL/interface/PatternRecognitionAlgoBase.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

namespace ticl {
  template <typename TILES>
  class PatternRecognitionbyCLUE3D final : public PatternRecognitionAlgoBaseT<TILES> {
  public:
    PatternRecognitionbyCLUE3D(const edm::ParameterSet& conf, edm::ConsumesCollector);
    ~PatternRecognitionbyCLUE3D() override = default;

    void makeTracksters(const typename PatternRecognitionAlgoBaseT<TILES>::Inputs& input,
                        std::vector<Trackster>& result,
                        std::unordered_map<int, std::vector<int>>& seedToTracksterAssociation) override;

    void energyRegressionAndID(const std::vector<reco::CaloCluster>& layerClusters,
                               const tensorflow::Session*,
                               std::vector<Trackster>& result);

    static void fillPSetDescription(edm::ParameterSetDescription& iDesc);

  private:
    /**
     * Information of all layer clusters on a given layer (as a Struct of Arrays)
    */
    struct ClustersOnLayer {
      std::vector<float> x;
      std::vector<float> y;
      std::vector<float> z;
      std::vector<float> r_over_absz;
      std::vector<float> radius;
      std::vector<float> eta;
      std::vector<float> phi;
      std::vector<int> cells;
      std::vector<uint8_t> isSilicon;

      std::vector<float> energy;
      std::vector<float> rho;
      std::vector<float> z_extension;///< depth in z where layer clusters are searched, taking into account boundaries of detector

      /**
       * For each layer cluster : 
       * pair.first : distance to nearest higher (how it is computed depends on parameters)
       * pair.second : number of layers of difference to the nearest higher
      */
      std::vector<std::pair<float, int>> delta;
      /**
       * For each layer cluster
       * pair.first : layer ID of nearest higher
       * pair.second : ID of the nearest higher (in local-to-layer ID)
      */
      std::vector<std::pair<int, int>> nearestHigher;
      std::vector<int> clusterIndex;
      std::vector<unsigned int> layerClusterOriginalIdx;
      std::vector<std::vector<std::pair<int, int>>> followers;
      std::vector<bool> isSeed;

      void clear() {
        x.clear();
        y.clear();
        z.clear();
        r_over_absz.clear();
        radius.clear();
        eta.clear();
        phi.clear();
        cells.clear();
        isSilicon.clear();
        energy.clear();
        rho.clear();
        z_extension.clear();
        delta.clear();
        nearestHigher.clear();
        clusterIndex.clear();
        layerClusterOriginalIdx.clear();
        followers.clear();
        isSeed.clear();
      }

      void shrink_to_fit() {
        x.shrink_to_fit();
        y.shrink_to_fit();
        z.shrink_to_fit();
        r_over_absz.shrink_to_fit();
        radius.shrink_to_fit();
        eta.shrink_to_fit();
        phi.shrink_to_fit();
        cells.shrink_to_fit();
        isSilicon.shrink_to_fit();
        energy.shrink_to_fit();
        rho.shrink_to_fit();
        z_extension.shrink_to_fit();
        delta.shrink_to_fit();
        nearestHigher.shrink_to_fit();
        clusterIndex.shrink_to_fit();
        layerClusterOriginalIdx.shrink_to_fit();
        followers.shrink_to_fit();
        isSeed.shrink_to_fit();
      }
    };

    void reset() {
      for (auto& c : clusters_) {
        c.clear();
        c.shrink_to_fit();
      }
    }

    /**
     * Compute local density of all layer clusters on the given layer
     * \param layerIdx2layerandSoa Vector holding (layer ID; ID of cluster on layer) indexed by input.layerClusters ID. If an entry is (-1, -1) then cluster is masked
    */
    void calculateLocalDensity(const TILES& tiles, const int layerId, const std::vector<std::pair<int, int>>& layerIdx2layerandSoa);
    
    /**
     * Compute distance to nearest higher of all layer clusters on the given layer
     * \param layerIdx2layerandSoa Vector holding (layer ID; ID of cluster on layer) indexed by input.layerClusters ID. If an entry is (-1, -1) then cluster is masked
    */
    void calculateDistanceToHigher(const TILES&, const int layerId, const std::vector<std::pair<int, int>>& layerIdx2layerandSoa);
    
    /**
     * For all layer clusters, compute whether it is a seed/outlier/follower, and expand tracksters from seeds
     * \param layerIdx2layerandSoa Vector holding (layer ID; ID of cluster on layer) indexed by input.layerClusters ID. If an entry is (-1, -1) then cluster is masked
    */
    int findAndAssignTracksters(const TILES&, const std::vector<std::pair<int, int>>& layerIdx2layerandSoa);
    void dumpClusters(const TILES& tiles,
                      const std::vector<std::pair<int, int>>& layerIdx2layerandSoa,
                      const int) const;
    void dumpTracksters(const std::vector<std::pair<int, int>>& layerIdx2layerandSoa,
                        const int,
                        const std::vector<Trackster>&) const;
    void dumpTiles(const TILES&) const;

    /**
     * Vector indexed by layer number holding information of all layer clusters on this layer
    */
    std::vector<ClustersOnLayer> clusters_;
    std::vector<float> layersPosZ_;

    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
    const double criticalDensity_; ///< Threshold of density for seeding (default 4 GeV)
    const double criticalSelfDensity_;
    const int densitySiblingLayers_; ///< Layer window for density calculation (default 3)
    const double densityEtaPhiDistanceSqr_; ///< Square of distance threshold in (eta; phi) plane for density calculation
    const double densityXYDistanceSqr_;
    const double kernelDensityFactor_;
    const bool densityOnSameLayer_; ///< Whether to consider 2D clusters on same layer to make a trackster
    const bool nearestHigherOnSameLayer_;
    const bool useAbsoluteProjectiveScale_;
    const bool useClusterDimensionXY_;
    const bool rescaleDensityByZ_;
    const double criticalEtaPhiDistance_; ///< if the distance to higher density cell is larger than this value, the cell will be promoted as a seed.
    const double criticalXYDistance_;
    const int criticalZDistanceLyr_;
    const double outlierMultiplier_; ///< Factor to multiply deltac (critical distance) to determine if 2D cluster is an outlier 
    const int minNumLayerCluster_; ///< Minimum number of layer clusters (default 5)
    const std::vector<int> filter_on_categories_;
    const std::string eidInputName_;
    const std::string eidOutputNameEnergy_;
    const std::string eidOutputNameId_;
    const float eidMinClusterEnergy_;
    const int eidNLayers_;
    const int eidNClusters_;

    hgcal::RecHitTools rhtools_;
    tensorflow::Session* eidSession_;

    static const int eidNFeatures_ = 3;
  };

}  // namespace ticl
#endif
