#ifndef RecoLocalCalo_HGCalRecProducers_HGCalCLUEAlgo_h
#define RecoLocalCalo_HGCalRecProducers_HGCalCLUEAlgo_h

#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalClusteringAlgoBase.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalLayerTiles.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

// C/C++ headers
#include <set>
#include <string>
#include <vector>

using Density = hgcal_clustering::Density; // std::map<DetId, float>

template <typename TILE>
class HGCalCLUEAlgoT : public HGCalClusteringAlgoBase {
public:
  HGCalCLUEAlgoT(const edm::ParameterSet& ps, edm::ConsumesCollector iC)
      : HGCalClusteringAlgoBase(
            (HGCalClusteringAlgoBase::VerbosityLevel)ps.getUntrackedParameter<unsigned int>("verbosity", 3),
            reco::CaloCluster::undefined,
            iC),
        thresholdW0_(ps.getParameter<std::vector<double>>("thresholdW0")),
        positionDeltaRho2_(ps.getParameter<double>("positionDeltaRho2")),
        vecDeltas_(ps.getParameter<std::vector<double>>("deltac")),
        kappa_(ps.getParameter<double>("kappa")),
        ecut_(ps.getParameter<double>("ecut")),
        dependSensor_(ps.getParameter<bool>("dependSensor")),
        dEdXweights_(ps.getParameter<std::vector<double>>("dEdXweights")),
        thicknessCorrection_(ps.getParameter<std::vector<double>>("thicknessCorrection")),
        sciThicknessCorrection_(ps.getParameter<double>("sciThicknessCorrection")),
        deltasi_index_regemfac_(ps.getParameter<int>("deltasi_index_regemfac")),
        maxNumberOfThickIndices_(ps.getParameter<unsigned>("maxNumberOfThickIndices")),
        fcPerMip_(ps.getParameter<std::vector<double>>("fcPerMip")),
        fcPerEle_(ps.getParameter<double>("fcPerEle")),
        nonAgedNoises_(ps.getParameter<std::vector<double>>("noises")),
        noiseMip_(ps.getParameter<edm::ParameterSet>("noiseMip").getParameter<double>("noise_MIP")),
        use2x2_(ps.getParameter<bool>("use2x2")),
        initialized_(false) {}

  ~HGCalCLUEAlgoT() override {}

  void getEventSetupPerAlgorithm(const edm::EventSetup& es) override;

  void populate(const HGCRecHitCollection& hits) override;

  // this is the method that will start the clusterisation (it is possible to invoke this method
  // more than once - but make sure it is with different hit collections (or else use reset)

  void makeClusters() override;

  ///< this is the method to get the cluster collection out (parameter is ignored)
  std::vector<reco::BasicCluster> getClusters(bool) override;

  void reset() override {
    clusters_v_.clear();
    clusters_v_.shrink_to_fit();
    for (auto& cl : numberOfClustersPerLayer_) {
      cl = 0;
    }

    for (auto& cells : cells_) {
      cells.clear();
      cells.shrink_to_fit();
    }
    density_.clear();
  }

  Density getDensity() override;

  void computeThreshold();

  static void fillPSetDescription(edm::ParameterSetDescription& iDesc) {
    iDesc.add<std::vector<double>>("thresholdW0", {2.9, 2.9, 2.9});
    iDesc.add<double>("positionDeltaRho2", 1.69);
    iDesc.add<std::vector<double>>("deltac",
                                   {
                                       1.3,
                                       1.3,
                                       1.3,
                                       0.0315,  // for scintillator
                                   });
    iDesc.add<bool>("dependSensor", true);
    iDesc.add<double>("ecut", 3.0);
    iDesc.add<double>("kappa", 9.0);
    iDesc.addUntracked<unsigned int>("verbosity", 3);
    iDesc.add<std::vector<double>>("dEdXweights", {});
    iDesc.add<std::vector<double>>("thicknessCorrection", {});
    iDesc.add<double>("sciThicknessCorrection", 0.9);
    iDesc.add<int>("deltasi_index_regemfac", 3);
    iDesc.add<unsigned>("maxNumberOfThickIndices", 6);
    iDesc.add<std::vector<double>>("fcPerMip", {});
    iDesc.add<double>("fcPerEle", 0.0);
    iDesc.add<std::vector<double>>("noises", {});
    edm::ParameterSetDescription descNestedNoiseMIP;
    descNestedNoiseMIP.add<bool>("scaleByDose", false);
    descNestedNoiseMIP.add<unsigned int>("scaleByDoseAlgo", 0);
    descNestedNoiseMIP.add<double>("scaleByDoseFactor", 1.);
    descNestedNoiseMIP.add<std::string>("doseMap", "");
    descNestedNoiseMIP.add<std::string>("sipmMap", "");
    descNestedNoiseMIP.add<double>("referenceIdark", -1);
    descNestedNoiseMIP.add<double>("referenceXtalk", -1);
    descNestedNoiseMIP.add<double>("noise_MIP", 1. / 100.);
    iDesc.add<edm::ParameterSetDescription>("noiseMip", descNestedNoiseMIP);
    iDesc.add<bool>("use2x2", true);  // use 2x2 or 3x3 scenario for scint density calculation
  }

  /// point in the space
  typedef math::XYZPoint Point;

private:
  // To compute the cluster position
  std::vector<double> thresholdW0_;
  /**
   * When computing layer cluster position,
   *  ignore cells that have the distance squared to the highest energy cell in layer cluster greater than this squared distance
  */
  const double positionDeltaRho2_;

  // The two parameters used to identify clusters
  std::vector<double> vecDeltas_;
  double kappa_; ///< Minimum value, in terms of sigma(noise), to promote a hit to become a seed candidate for cluster (default is 9)

  // The hit energy cutoff
  double ecut_;

  // For keeping the density per hit
  Density density_;

  // various parameters used for calculating the noise levels for a given sensor (and whether to use
  // them)
  bool dependSensor_;
  std::vector<double> dEdXweights_;
  std::vector<double> thicknessCorrection_;
  double sciThicknessCorrection_;
  int deltasi_index_regemfac_;
  unsigned maxNumberOfThickIndices_;
  std::vector<double> fcPerMip_;
  double fcPerEle_;
  std::vector<double> nonAgedNoises_;
  double noiseMip_;
  std::vector<std::vector<double>> thresholds_;
  std::vector<std::vector<double>> v_sigmaNoise_;

  bool use2x2_;

  // initialization bool
  bool initialized_;

  /**
   * Factor to multiply deltac to determine if hit is an outlier 
   * Has to be greater than 1 (calculateDistanceToHigher works properly only when this is the case)
  */
  float outlierDeltaFactor_ = 2.f; 

  struct CellsOnLayer {
    std::vector<DetId> detid;
    std::vector<bool> isSi;
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> eta;
    std::vector<float> phi;

    std::vector<float> weight;
    std::vector<float> rho;

    std::vector<float> delta;
    std::vector<int> nearestHigher;
    std::vector<int> clusterIndex;
    std::vector<float> sigmaNoise;
    std::vector<std::vector<int>> followers;
    std::vector<bool> isSeed;

    void clear() {
      detid.clear();
      isSi.clear();
      x.clear();
      y.clear();
      eta.clear();
      phi.clear();
      weight.clear();
      rho.clear();
      delta.clear();
      nearestHigher.clear();
      clusterIndex.clear();
      sigmaNoise.clear();
      followers.clear();
      isSeed.clear();
    }

    void shrink_to_fit() {
      detid.shrink_to_fit();
      isSi.shrink_to_fit();
      x.shrink_to_fit();
      y.shrink_to_fit();
      eta.shrink_to_fit();
      phi.shrink_to_fit();
      weight.shrink_to_fit();
      rho.shrink_to_fit();
      delta.shrink_to_fit();
      nearestHigher.shrink_to_fit();
      clusterIndex.shrink_to_fit();
      sigmaNoise.shrink_to_fit();
      followers.shrink_to_fit();
      isSeed.shrink_to_fit();
    }
  };

  std::vector<CellsOnLayer> cells_; ///< List of CellsOnLayer indexed by layer ID

  std::vector<int> numberOfClustersPerLayer_; ///< Computed by findAndAssignClusters

  /**
   * Compute distance squared between two cells on the same layer
   * \param isEtaPhi if true, compute distance in (eta, phi) ie eta^2 + phi^2. If false, compute in (x;y)
  */
  inline float distance2(int cell1, int cell2, int layerId, bool isEtaPhi) const {  // distance squared
    if (isEtaPhi) {
      const float dphi = reco::deltaPhi(cells_[layerId].phi[cell1], cells_[layerId].phi[cell2]);
      const float deta = cells_[layerId].eta[cell1] - cells_[layerId].eta[cell2];
      return (deta * deta + dphi * dphi);
    } else {
      const float dx = cells_[layerId].x[cell1] - cells_[layerId].x[cell2];
      const float dy = cells_[layerId].y[cell1] - cells_[layerId].y[cell2];
      return (dx * dx + dy * dy);
    }
  }

  ///< Compute 2D distance between two cells on the same layer. Square root of \a HCCalCLUEAlgoT::distance2
  inline float distance(int cell1, int cell2, int layerId, bool isEtaPhi) const {
    return std::sqrt(distance2(cell1, cell2, layerId, isEtaPhi));
  }

  void prepareDataStructures(const unsigned int layerId);

  /**
   * Compute rho, the local energy density, for all cells in a given layer (ie fill cells_[layerId].rho[all cells])
   * Computed differently for sillicon and scintillator cells 
   * \param delta_c distance parameter for sillicon
   * \param delta_r distance parameter for scintillator
  */
  void calculateLocalDensity(const TILE& lt,
                             const unsigned int layerId,
                             float delta_c,
                             float delta_r);
  
  /**
   * Compute delta (distance to nearest higher) and nearestHigher (ID of nearest higher) for all cells in given layer
   * (ie fills cells_[layerId].delta[all cells] and cells_[layerId].nearestHigher[all cells])
   * The range to search for is outlierDeltaFactor_ * delta_* (depends on sillicon/scintillator)
   * \param delta_c distance parameter for sillicon
   * \param delta_r distance parameter for scintillator
  */
  void calculateDistanceToHigher(const TILE& lt, const unsigned int layerId, float delta_c, float delta_r);

  /**
   * For all hits on a given layer, compute whether it is a seed, outlier, or follower (in this case register to the nearest higher)
   * Then expand clusters from seeds (setting point.clusterIndex for all points in each cluster)
   * \param delta_c distance parameter for sillicon
   * \param delta_r distance parameter for scintillator
   * Depends also on \a kappa_ , sigmaNoise of each cell, and outlierDeltaFactor_
  */
  int findAndAssignClusters(const unsigned int layerId, float delta_c, float delta_r);

  /**
   * COmpute layer cluster position
   * \param v list of cells id in a cluster
  */
  math::XYZPoint calculatePosition(const std::vector<int>& v, const unsigned int layerId) const;

  ///< Fill density_ from computed rho (computed in calculateLocalDensity) for a given layer
  void setDensity(const unsigned int layerId);
};

// explicit template instantiation
extern template class HGCalCLUEAlgoT<HGCalLayerTiles>;
extern template class HGCalCLUEAlgoT<HFNoseLayerTiles>;

using HGCalCLUEAlgo = HGCalCLUEAlgoT<HGCalLayerTiles>;
using HFNoseCLUEAlgo = HGCalCLUEAlgoT<HFNoseLayerTiles>;

#endif
