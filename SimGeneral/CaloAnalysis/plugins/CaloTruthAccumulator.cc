#define DEBUG false

#if DEBUG
#pragma GCC diagnostic pop
#endif

#include <iterator>
#include <algorithm>
#include <memory>
#include <vector>
#include <map>
#include <unordered_map>
#include <string>
#include <tuple>

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ProducesCollector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalTestNumbering.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloTest/interface/HGCalTestNumbering.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Track/interface/SimTrack.h"

#include "SimGeneral/MixingModule/interface/DecayGraph.h"
#include "SimGeneral/MixingModule/interface/DigiAccumulatorMixMod.h"
#include "SimGeneral/MixingModule/interface/DigiAccumulatorMixModFactory.h"
#include "SimGeneral/MixingModule/interface/PileUpEventPrincipal.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/HcalCommonData/interface/HcalHitRelabeller.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include <CLHEP/Units/SystemOfUnits.h>

namespace {
  using Index_t = unsigned;
  using Barcode_t = int;
  const std::string messageCategoryGraph_("CaloTruthAccumulatorGraphProducer");
}  // namespace

class CaloTruthAccumulator : public DigiAccumulatorMixMod {
public:
  explicit CaloTruthAccumulator(const edm::ParameterSet &config, edm::ProducesCollector, edm::ConsumesCollector &iC);

private:
  void initializeEvent(const edm::Event &event, const edm::EventSetup &setup) override;
  void accumulate(const edm::Event &event, const edm::EventSetup &setup) override;
  void accumulate(const PileUpEventPrincipal &event, const edm::EventSetup &setup, edm::StreamID const &) override;
  void finalizeEvent(edm::Event &event, const edm::EventSetup &setup) override;

  /** @brief Both forms of accumulate() delegate to this templated method. */
  template <class T>
  void accumulateEvent(const T &event,
                       const edm::EventSetup &setup,
                       const edm::Handle<edm::HepMCProduct> &hepMCproduct);

  /**
 * @brief Fills the supplied vector with pointers to the SimHits, checking or bad modules if required.
 * Iterate over all SimHits collections (vectors of PCaloHit) and build maps DetId->Hit
 * 
 * @param[out] returnValue : return parameter holding vector of (DetId, PCaloHit*) pairs
 * @param[out] simTrackDetIdEnergyMap : return parameter, nested map G4 Track ID -> DetId -> accumulated SimHit energy (over same track, in case of loops)
 * 
 * Also this->m_detIdToTotalSimEnergy is filled (map DetId-> accumulated sim energy), for normalization into fractions in finalizeEvent
 */
  template <class T>
  void fillSimHits(std::vector<std::pair<DetId, const PCaloHit *>> &returnValue,
                   std::unordered_map<int, std::map<int, float>> &simTrackDetIdEnergyMap,
                   const T &event,
                   const edm::EventSetup &setup);

  const std::string messageCategory_;

  std::unordered_map<Index_t, float> m_detIdToTotalSimEnergy;  // keep track of cell normalizations
  std::unordered_multimap<Barcode_t, Index_t> m_simHitBarcodeToIndex;

  /** The maximum bunch crossing BEFORE the signal crossing to create
      TrackinParticles for. Use positive values. If set to zero no
      previous bunches are added and only in-time, signal and after bunches
      (defined by maximumSubsequentBunchCrossing_) are used.
  */
  const unsigned int maximumPreviousBunchCrossing_;
  /** The maximum bunch crossing AFTER the signal crossing to create
      TrackinParticles for. E.g. if set to zero only
      uses the signal and in time pileup (and previous bunches defined by the
      maximumPreviousBunchCrossing_ parameter).
  */
  const unsigned int maximumSubsequentBunchCrossing_;

  const edm::InputTag simTrackLabel_;
  const edm::InputTag simVertexLabel_;

  std::vector<edm::InputTag> collectionTags_;
  edm::InputTag genParticleLabel_;
  /// Needed to add HepMC::GenVertex to SimVertex
  edm::InputTag hepMCproductLabel_;
  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geomToken_;
  edm::ESWatcher<CaloGeometryRecord> geomWatcher_;

  const double minEnergy_, maxPseudoRapidity_;
  const bool premixStage1_;
  
  CaloParticleCollection outputCaloParticles_;
  edm::EDPutTokenT<CaloParticleCollection> outputCaloParticles_token_;

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
  SimClusterConfig legacySimClusters_config_;   ///< Legacy SimCluster from every SimTrack with simhits
  SimClusterConfig boundarySimClusters_config_; ///< SimClusters from each SimTrack crossing boundary
  SimClusterConfig caloParticleSimClusters_config_; ///< SimCluster that are identical to CaloParticle (to make it easier on downstream code, only one dataformat for everything)

  edm::RefProd<SimClusterCollection> caloParticles_refProd_;  // To build the RefVector from SimCluster to "CaloParticle" SimCluster collection

  const HGCalTopology *hgtopo_[3] = {nullptr, nullptr, nullptr};
  const HGCalDDDConstants *hgddd_[3] = {nullptr, nullptr, nullptr};
  const HcalDDDRecConstants *hcddd_ = nullptr;
  // geometry type (0 pre-TDR; 1 TDR)
  int geometryType_;
  bool doHGCAL;
};

/* Graph utility functions */

namespace {
  class VisitorHelper {
  public:
    VisitorHelper(std::unordered_multimap<Barcode_t, Index_t> const &simHitBarcodeToIndex,
                  std::unordered_map<int, std::map<int, float>> const &simTrackDetIdEnergyMap,
                  std::unordered_map<uint32_t, float> const &vertex_time_map)
        : simHitBarcodeToIndex_(simHitBarcodeToIndex),
          simTrackDetIdEnergyMap_(simTrackDetIdEnergyMap),
          vertex_time_map_(vertex_time_map) {}

    // bool simTrackHasSimHits(unsigned int trackIdx) const { return simHitBarcodeToIndex_.count(trackIdx); } // Not needed as already saved in EdgeProperty
    /** Given a track G4 ID, give map DetId->accumulated SimHit energy  */
    std::map<int, float> const &hits_and_energies(unsigned int trackIdx) const {
      return simTrackDetIdEnergyMap_.at(trackIdx);
    }

    float getVertexTime(uint32_t simVertexId) { return vertex_time_map_.at(simVertexId); }

  private:
    std::unordered_multimap<Barcode_t, Index_t> const &
        simHitBarcodeToIndex_;  // Map G4 track ID -> index in simHitPointers collection holding (DetId, PCaloHit*) pair. Used to quickly determine if a track has simhits
    std::unordered_map<int, std::map<int, float>> const &
        simTrackDetIdEnergyMap_;  // nested map G4 Track ID -> DetId -> accumulated SimHit energy (over same track, in case of loops)
    std::unordered_map<uint32_t, float> const &vertex_time_map_;  // Map SimVertex index to sim time
  };

  /** 
   * Accumlates SimHits energies and DetId as each graph edge (=SimTrack) is visited
   * @tparam SimClassNameT either SimCluster or CaloParticle
   */
  template <typename SimClassNameT>
  class ClusterEnergyAccumulator {
  public:
    ClusterEnergyAccumulator(VisitorHelper const &helper) : helper_(helper) {}

    template <class Edge, class Graph>
    void accumulate_edge_in_cluster(Edge e, const Graph &g) {
      auto edge_property = get(edge_weight, g, e);
      const SimTrack *edge_simTrack = edge_property.simTrack;

      if (!edge_simTrack)
        return;  // Should not happen
      auto trackIdx = edge_simTrack->trackId();
      if (edge_property.simHits != 0) {
        for (auto const &hit_and_energy : helper_.hits_and_energies(trackIdx)) {
          acc_energy[hit_and_energy.first] += hit_and_energy.second;
          // cluster.addSimHit()
        }
      }
    }

    void doEndCluster(SimClassNameT &cluster) {
      for (auto const &hit_and_energy : acc_energy) {
        cluster.addRecHitAndFraction(hit_and_energy.first, hit_and_energy.second);
      }
      acc_energy.clear();
    }

  private:
    VisitorHelper const &helper_;
    std::unordered_map<uint32_t, float>
        acc_energy;  // Map DetId->simHit energies for the current cluster (only vaild if insideCluster_)
  };

  /**
   * Visitor class for depth_first_search.
   *  - on a root edge, check the associated GenParticle. Check if it passes selection
   *     if it does not, discard 
   *     if it does, start a new CaloParticle (record edge)
   *  - on any edge within a CaloParticle, check SimCluster conditions. 
   *    If they pass, create a SimCluster. 
   * @tparam Selector_t lambda for selections to create a CaloParticle
   * @tparam SubClusterVisitorTuple std::tuple of visitors for creating SimCluster collections nested inside the CaloParticle
   */
  template <typename Selector_t, typename SubClusterVisitorTuple>
  class PrimaryVisitor : public boost::default_dfs_visitor {
  public:
    /**
   * Use like PrimaryVisitor(helper, caloParticles, [](auto edgeProperty){return ...;}, 
   *  std::make_tuple(SubClusterVisitor<SimCluster>(legacySimClusters, ....), (SubClusterVisitor<SimCluster> boundarySimClusters, ....)
   */
    PrimaryVisitor(VisitorHelper &helper,
                   std::vector<CaloParticle> &caloParticles,
                   Selector_t caloParticleSelector,
                   SubClusterVisitorTuple subClusterVisitors)
        : helper_(helper),
          caloParticleSelector_(caloParticleSelector),
          caloParticleAccumulator_(helper),
          caloParticles_(caloParticles),
          subClusterVisitors_(subClusterVisitors) {}

    void examine_edge(DecayChain::edge_descriptor e, const DecayChain &g) {
      auto edge_property = get(edge_weight, g, e);
      if (!insideCluster_ && caloParticleSelector_(edge_property)) {
        insideCluster_ = true;
        caloParticleEdge_ = e;
        // Create a new CaloParticle
        caloParticles_.emplace_back(*edge_property.simTrack);
        caloParticles_.back().setSimTime(helper_.getVertexTime(
            edge_property.simTrack->vertIndex()));  // For a CaloParticle the simTime is set at the vertex !
      }

      if (insideCluster_) {
        caloParticleAccumulator_.accumulate_edge_in_cluster(e, g); // accumulate simhits energies

        // This loops over all elements of subClusterVisitors_ tuple and calls examine_edge on them
        // Another approach could be (compile-time) recursion with nested std::pair
        std::apply([&](auto &...x) { (..., x.examine_edge(e, g, caloParticles_.size() - 1)); }, subClusterVisitors_);
      }
    }

    void finish_edge(DecayChain::edge_descriptor e, const DecayChain &g) {
      if (insideCluster_) {
        // This loops over all elements of subClusterVisitors_ tuple and calls finish_edge on them
        std::apply([&](auto &...x) { (..., x.finish_edge(e, g)); }, subClusterVisitors_);
      }
      if (insideCluster_ && e == caloParticleEdge_) {
        insideCluster_ = false;
        caloParticleAccumulator_.doEndCluster(caloParticles_.back());
      }
    }

  private:
    VisitorHelper &helper_;

    Selector_t caloParticleSelector_;
    ClusterEnergyAccumulator<CaloParticle> caloParticleAccumulator_;
    std::vector<CaloParticle> &caloParticles_; // Output collection

    SubClusterVisitorTuple subClusterVisitors_;  // std::tuple of SubClusterVisitor

    bool insideCluster_{false};
    DecayChain::edge_descriptor caloParticleEdge_; // Keep track of CaloParticle root edge to finish cluster creation
  };

  /** 
   * Fills RefVector to keep mapping from SimCluster to another collection (typically CaloParticle).
   * @tparam ParentClusterCollectionT type of the collection to refer to (eg: std::vector<SimCluster> or std::vector<CaloParticle>)
   */
  template <typename ParentClusterCollectionT>
  class ClusterParentIndexRecorder {
  public:
    /** The RefProd is to fill the RefVector with edm::Ref (build it with getRefBeforePut) */
    ClusterParentIndexRecorder(edm::RefVector<ParentClusterCollectionT> &refVector,
                               edm::RefProd<ParentClusterCollectionT> const &refProd)
        : refVector_(refVector), refProd_(refProd) {}

    void recordParentClusterIndex(std::size_t parentClusterIndex) {
      refVector_.push_back(edm::Ref<ParentClusterCollectionT>(refProd_, parentClusterIndex));
    }

  private:
    edm::RefVector<ParentClusterCollectionT> &refVector_;
    edm::RefProd<ParentClusterCollectionT> const
        &refProd_;  // Need the RefProd to build the edm::Ref to insert into RefVector
  };

  /** Does nothing */
  class ClusterParentIndexRecorderVoid {
    public:
      void recordParentClusterIndex(std::size_t parentClusterIndex) {}
  };

  /** 
   * Visitor that creates cluster with reference to another "parent" cluster collection, like SimCluster inside CaloParticle
   * Does not make "nested sub-clusters"
   * @tparam SubClusterT typically SimCluster
   * @tparam ClusterParentIndexRecorderT type of the object in charge of building the RefVector to parent (@see ClusterParentIndexRecorder)
   * @tparam Selector_t lambda for criteria for creating a subcluster
   */
  template <typename SubClusterT, typename ClusterParentIndexRecorderT, typename Selector_t>
  class SubClusterVisitor {
  public:
    SubClusterVisitor(std::vector<SubClusterT> &clusters,
                      ClusterParentIndexRecorderT parentIndexRecorder,
                      VisitorHelper const &helper,
                      Selector_t selector)
        : clusters_(clusters), accumulator(helper), indexRecorder(parentIndexRecorder), selector_(selector) {}

    void examine_edge(DecayChain::edge_descriptor e, const DecayChain &g, std::size_t currentCaloParticleIndex) {
      if (!insideCluster_ && selector_(get(edge_weight, g, e))) {
        insideCluster_ = true;
        clusterRootEdge_ = e;

        // Create a new cluster
        auto edge_property = get(edge_weight, g, e);
        clusters_.emplace_back(*edge_property.simTrack);
        indexRecorder.recordParentClusterIndex(currentCaloParticleIndex);
      }

      if (insideCluster_) {
        accumulator.accumulate_edge_in_cluster(e, g);
      }
    }

    void finish_edge(DecayChain::edge_descriptor e, const DecayChain &g) {
      if (insideCluster_ && e == clusterRootEdge_) {
        insideCluster_ = false;
        accumulator.doEndCluster(clusters_.back());
      }
    }

  private:
    std::vector<SubClusterT> &clusters_;
    ClusterEnergyAccumulator<SubClusterT> accumulator;
    ClusterParentIndexRecorderT indexRecorder;
    Selector_t selector_;
    bool insideCluster_{false};
    DecayChain::edge_descriptor clusterRootEdge_;
  };

  /**
   * Creates SimCluster for every SimTrack with simHits. One SimCluster only includes the simHits from the SimTrack it was made of (depends on SimTrack saving criteria)
   */
  template <typename SubClusterT, typename ClusterParentIndexRecorderT>
  class LegacySimClusterVisitor {
  public:
    LegacySimClusterVisitor(std::vector<SubClusterT> &clusters,
                      ClusterParentIndexRecorderT parentIndexRecorder,
                      VisitorHelper const &helper)
        : clusters_(clusters), accumulator(helper), indexRecorder(parentIndexRecorder) {}

    void examine_edge(DecayChain::edge_descriptor e, const DecayChain &g, std::size_t currentCaloParticleIndex) {
      auto edge_property = get(edge_weight, g, e);
      if (edge_property.simHits != 0) {
        // Create a new cluster
        clusters_.emplace_back(*edge_property.simTrack);
        indexRecorder.recordParentClusterIndex(currentCaloParticleIndex);
        accumulator.accumulate_edge_in_cluster(e, g);
        accumulator.doEndCluster(clusters_.back());
      }
    }
    void finish_edge(DecayChain::edge_descriptor e, const DecayChain &g) {}

  private:
    std::vector<SubClusterT> &clusters_;
    ClusterEnergyAccumulator<SubClusterT> accumulator;
    ClusterParentIndexRecorderT indexRecorder;
  };

}  // namespace

CaloTruthAccumulator::CaloTruthAccumulator(const edm::ParameterSet &config,
                                           edm::ProducesCollector producesCollector,
                                           edm::ConsumesCollector &iC)
    : messageCategory_("CaloTruthAccumulator"),
      maximumPreviousBunchCrossing_(config.getParameter<unsigned int>("maximumPreviousBunchCrossing")),
      maximumSubsequentBunchCrossing_(config.getParameter<unsigned int>("maximumSubsequentBunchCrossing")),
      simTrackLabel_(config.getParameter<edm::InputTag>("simTrackCollection")),
      simVertexLabel_(config.getParameter<edm::InputTag>("simVertexCollection")),
      collectionTags_(),
      genParticleLabel_(config.getParameter<edm::InputTag>("genParticleCollection")),
      hepMCproductLabel_(config.getParameter<edm::InputTag>("HepMCProductLabel")),
      geomToken_(iC.esConsumes()),
      minEnergy_(config.getParameter<double>("MinEnergy")),
      maxPseudoRapidity_(config.getParameter<double>("MaxPseudoRapidity")),
      premixStage1_(config.getParameter<bool>("premixStage1")),
      outputCaloParticles_token_(producesCollector.produces<CaloParticleCollection>("MergedCaloTruth")),
      legacySimClusters_config_(producesCollector, "MergedCaloTruth"),
      boundarySimClusters_config_(producesCollector, "MergedCaloTruthBoundaryTrackSimCluster"),
      caloParticleSimClusters_config_(producesCollector, "MergedCaloTruthCaloParticle"),
      geometryType_(-1),
      doHGCAL(config.getParameter<bool>("doHGCAL")) {
  if (premixStage1_) {
    producesCollector.produces<std::vector<std::pair<unsigned int, float>>>(
        "MergedCaloTruth");  // (DetId, total simhit energy) pairs
  }

  iC.consumes<std::vector<SimTrack>>(simTrackLabel_);
  iC.consumes<std::vector<SimVertex>>(simVertexLabel_);
  iC.consumes<std::vector<int>>(genParticleLabel_);
  iC.consumes<std::vector<int>>(hepMCproductLabel_);

  // Fill the collection tags
  const edm::ParameterSet &simHitCollectionConfig = config.getParameterSet("simHitCollections");
  std::vector<std::string> parameterNames = simHitCollectionConfig.getParameterNames();

  for (auto const &parameterName : parameterNames) {
    std::vector<edm::InputTag> tags = simHitCollectionConfig.getParameter<std::vector<edm::InputTag>>(parameterName);
    collectionTags_.insert(collectionTags_.end(), tags.begin(), tags.end());
  }

  for (auto const &collectionTag : collectionTags_) {
    iC.consumes<std::vector<PCaloHit>>(collectionTag);
  }
}

void CaloTruthAccumulator::initializeEvent(edm::Event const &event, edm::EventSetup const &setup) {
  outputCaloParticles_.clear();
  legacySimClusters_config_.clear();
  boundarySimClusters_config_.clear();
  caloParticleSimClusters_config_.clear();

  m_detIdToTotalSimEnergy.clear();

  if (geomWatcher_.check(setup)) {
    auto const &geom = setup.getData(geomToken_);
    const HGCalGeometry *eegeom = nullptr, *fhgeom = nullptr, *bhgeomnew = nullptr;
    const HcalGeometry *bhgeom = nullptr;
    bhgeom = static_cast<const HcalGeometry *>(geom.getSubdetectorGeometry(DetId::Hcal, HcalEndcap));

    if (doHGCAL) {
      eegeom = static_cast<const HGCalGeometry *>(
          geom.getSubdetectorGeometry(DetId::HGCalEE, ForwardSubdetector::ForwardEmpty));
      // check if it's the new geometry
      if (eegeom) {
        geometryType_ = 1;
        fhgeom = static_cast<const HGCalGeometry *>(
            geom.getSubdetectorGeometry(DetId::HGCalHSi, ForwardSubdetector::ForwardEmpty));
        bhgeomnew = static_cast<const HGCalGeometry *>(
            geom.getSubdetectorGeometry(DetId::HGCalHSc, ForwardSubdetector::ForwardEmpty));
      } else {
        geometryType_ = 0;
        eegeom = static_cast<const HGCalGeometry *>(geom.getSubdetectorGeometry(DetId::Forward, HGCEE));
        fhgeom = static_cast<const HGCalGeometry *>(geom.getSubdetectorGeometry(DetId::Forward, HGCHEF));
        bhgeom = static_cast<const HcalGeometry *>(geom.getSubdetectorGeometry(DetId::Hcal, HcalEndcap));
      }
      hgtopo_[0] = &(eegeom->topology());
      hgtopo_[1] = &(fhgeom->topology());
      if (bhgeomnew)
        hgtopo_[2] = &(bhgeomnew->topology());

      for (unsigned i = 0; i < 3; ++i) {
        if (hgtopo_[i])
          hgddd_[i] = &(hgtopo_[i]->dddConstants());
      }
    }

    if (bhgeom) {
      hcddd_ = bhgeom->topology().dddConstants();
    }
  }
  // Seems const_cast is necessary here
  caloParticles_refProd_ =
      const_cast<edm::Event &>(event).getRefBeforePut<SimClusterCollection>(caloParticleSimClusters_config_.outputClusters_token);
}

/** Create handle to edm::HepMCProduct here because event.getByLabel with
    edm::HepMCProduct only works for edm::Event but not for
    PileUpEventPrincipal; PileUpEventPrincipal::getByLabel tries to call
    T::value_type and T::iterator (where T is the type of the object one wants
    to get a handle to) which is only implemented for container-like objects
    like std::vector but not for edm::HepMCProduct!
*/
void CaloTruthAccumulator::accumulate(edm::Event const &event, edm::EventSetup const &setup) {
  edm::Handle<edm::HepMCProduct> hepmc;
  event.getByLabel(hepMCproductLabel_, hepmc);

  edm::LogInfo(messageCategory_) << " CaloTruthAccumulator::accumulate (signal)";
  accumulateEvent(event, setup, hepmc);
}

void CaloTruthAccumulator::accumulate(PileUpEventPrincipal const &event,
                                      edm::EventSetup const &setup,
                                      edm::StreamID const &) {
  if (event.bunchCrossing() >= -static_cast<int>(maximumPreviousBunchCrossing_) &&
      event.bunchCrossing() <= static_cast<int>(maximumSubsequentBunchCrossing_)) {
    // simply create empty handle as we do not have a HepMCProduct in PU anyway
    edm::Handle<edm::HepMCProduct> hepmc;
    edm::LogInfo(messageCategory_) << " CaloTruthAccumulator::accumulate (pileup) bunchCrossing="
                                   << event.bunchCrossing();
    accumulateEvent(event, setup, hepmc);
  } else {
    edm::LogInfo(messageCategory_) << "Skipping pileup event for bunch crossing " << event.bunchCrossing();
  }
}

namespace {
/** Normalize the energies in the SimCluster/CaloParticle collection, from absolute SimHit energies to fraction of total simHits energies */
template <typename SimCaloCollection>
void normalizeCollection(SimCaloCollection &simClusters,
                         std::unordered_map<Index_t, float> const &detIdToTotalSimEnergy) {
  for (auto &sc : simClusters) {
    auto hitsAndEnergies = sc.hits_and_fractions();
    sc.clearHitsAndFractions();
    // sc.clearHitsEnergy();
    for (auto &hAndE : hitsAndEnergies) {
      const float totalenergy = detIdToTotalSimEnergy.at(hAndE.first);
      float fraction = 0.;
      if (totalenergy > 0)
        fraction = hAndE.second / totalenergy;
      else
        edm::LogWarning("CaloTruthAccumulator")
            << "TotalSimEnergy for hit " << hAndE.first << " is 0! The fraction for this hit cannot be computed.";
      sc.addRecHitAndFraction(hAndE.first, fraction);
      // sc.addHitEnergy(hAndE.second); // addHitEnergy is actually never used
    }
  }
}
};

void CaloTruthAccumulator::finalizeEvent(edm::Event &event, edm::EventSetup const &setup) {
  edm::LogInfo(messageCategory_) << "Adding " << legacySimClusters_config_.outputClusters.size() << " legacy SimClusters and "
                                 << outputCaloParticles_.size() << " CaloParticles to the event.";

  // We need to normalize the hits and energies into hits and fractions (since
  // we have looped over all pileup events)
  // For premixing stage1 we keep the energies, they will be normalized to
  // fractions in stage2 (in PreMixingCaloParticleWorker.cc)

  if (premixStage1_) {
    auto totalEnergies = std::make_unique<std::vector<std::pair<unsigned int, float>>>();
    totalEnergies->reserve(m_detIdToTotalSimEnergy.size());
    std::copy(m_detIdToTotalSimEnergy.begin(), m_detIdToTotalSimEnergy.end(), std::back_inserter(*totalEnergies));
    std::sort(totalEnergies->begin(), totalEnergies->end());
    event.put(std::move(totalEnergies), "MergedCaloTruth");
  } else {
    normalizeCollection(outputCaloParticles_, m_detIdToTotalSimEnergy);
    normalizeCollection(legacySimClusters_config_.outputClusters, m_detIdToTotalSimEnergy);
    normalizeCollection(boundarySimClusters_config_.outputClusters, m_detIdToTotalSimEnergy);
    normalizeCollection(caloParticleSimClusters_config_.outputClusters, m_detIdToTotalSimEnergy);
  }

  // fill the calo particles with their ref to SimCluster
  auto legacySimClustersRefProd = event.getRefBeforePut(legacySimClusters_config_.outputClusters_token);
  for (unsigned i = 0; i < legacySimClusters_config_.outputClusters.size(); ++i) {
    auto& cp = outputCaloParticles_[legacySimClusters_config_.clustersToCaloParticleMap[i].key()]; // get the key of the edm::Ref<CaloParticle>
    cp.addSimCluster(edm::Ref<SimClusterCollection>(legacySimClustersRefProd, i));
  }

  event.emplace(outputCaloParticles_token_, std::move(outputCaloParticles_));
  
  event.emplace(legacySimClusters_config_.outputClusters_token, std::move(legacySimClusters_config_.outputClusters));
  event.emplace(legacySimClusters_config_.clustersToCaloParticleMap_token, std::move(legacySimClusters_config_.clustersToCaloParticleMap));

  event.emplace(boundarySimClusters_config_.outputClusters_token, std::move(boundarySimClusters_config_.outputClusters));
  event.emplace(boundarySimClusters_config_.clustersToCaloParticleMap_token, std::move(boundarySimClusters_config_.clustersToCaloParticleMap));

  event.emplace(caloParticleSimClusters_config_.outputClusters_token, std::move(caloParticleSimClusters_config_.outputClusters));
  event.emplace(caloParticleSimClusters_config_.clustersToCaloParticleMap_token, std::move(caloParticleSimClusters_config_.clustersToCaloParticleMap)); // trivial mapping (kept for convenience)

  m_detIdToTotalSimEnergy.clear();
  m_simHitBarcodeToIndex.clear();
}

template <class T>
void CaloTruthAccumulator::accumulateEvent(const T &event,
                                           const edm::EventSetup &setup,
                                           const edm::Handle<edm::HepMCProduct> &hepMCproduct) {
  edm::Handle<std::vector<reco::GenParticle>> hGenParticles;
  edm::Handle<std::vector<int>> hGenParticleIndices;
  edm::Handle<std::vector<SimTrack>> hSimTracks;
  edm::Handle<std::vector<SimVertex>> hSimVertices;
  // We must always use getByLabel (event can be PileupEventPrincipal whihc does not have the get, getByToken, etc functions)
  event.getByLabel(simTrackLabel_, hSimTracks);
  event.getByLabel(simVertexLabel_, hSimVertices);

  event.getByLabel(genParticleLabel_, hGenParticles);
  event.getByLabel(genParticleLabel_, hGenParticleIndices);

  std::vector<std::pair<DetId, const PCaloHit *>> simHitPointers;
  std::unordered_map<int, std::map<int, float>> simTrackDetIdEnergyMap;
  fillSimHits(simHitPointers, simTrackDetIdEnergyMap, event, setup);

  // Clear maps from previous event fill them for this one
  m_simHitBarcodeToIndex.clear();
  for (unsigned int i = 0; i < simHitPointers.size(); ++i) {
    m_simHitBarcodeToIndex.emplace(simHitPointers[i].second->geantTrackId(), i);
  }

  auto const &tracks = *hSimTracks;
  auto const &vertices = *hSimVertices;
  std::unordered_map<int, int> trackid_to_track_index;
  DecayChain decay;
  int idx = 0;

  IfLogDebug(DEBUG, messageCategory_) << " TRACKS" << std::endl;
  for (auto const &t : tracks) {
    IfLogDebug(DEBUG, messageCategory_) << " " << idx << "\t" << t.trackId() << "\t" << t << std::endl;
    trackid_to_track_index[t.trackId()] = idx;
    idx++;
  }

  std::unordered_map<uint32_t, float> vertex_time_map;
  for (uint32_t i = 0; i < vertices.size(); i++) {
    // Geant4 time is in seconds, convert to ns (CLHEP::s = 1e9)
    vertex_time_map[i] = vertices[i].position().t() * CLHEP::s;
  }

  /**
  Build the main decay graph and assign the SimTrack to each edge. The graph
  built here will only contain the particles that have a decay vertex
  associated to them. In order to recover also the particles that will not
  decay, we need to keep track of the SimTrack used here and add, a-posteriori,
  the ones not used, associating a ghost vertex (starting from the highest
  simulated vertex number), in order to build the edge and identify them
  immediately as stable (i.e. not decayed).

  To take into account the multi-bremsstrahlung effects in which a single
  particle is emitting photons in different vertices **keeping the same
  track index**, we also collapsed those vertices into 1 unique vertex. The
  other approach of fully representing the decay chain keeping the same
  track index would have the problem of over-counting the contributions of
  that track, especially in terms of hits.

  The 2 auxiliary vectors are structured as follow:

  1. used_sim_tracks is a vector that has the same size as the overall
     number of simulated tracks. The associated integer is the vertexId of
     the **decaying vertex for that track**.
  2. collapsed_vertices is a vector that has the same size as the overall
     number of simulated vertices. The vector's index is the vertexId
     itself, the associated value is the vertexId of the vertex on which
     this should collapse.
  */
  idx = 0;
  std::vector<int> used_sim_tracks(tracks.size(), 0);
  std::vector<int> collapsed_vertices(vertices.size(), 0);
  IfLogDebug(DEBUG, messageCategory_) << " VERTICES" << std::endl;
  for (auto const &v : vertices) {
    IfLogDebug(DEBUG, messageCategory_) << " " << idx++ << "\t" << v << std::endl;
    if (v.parentIndex() != -1) {
      auto trk_idx = trackid_to_track_index[v.parentIndex()];
      auto origin_vtx = tracks[trk_idx].vertIndex();
      if (used_sim_tracks[trk_idx]) {
        // collapse the vertex into the original first vertex we saw associated
        // to this track. Omit adding the edge in order to avoid double
        // counting of the very same particles  and its associated hits.
        collapsed_vertices[v.vertexId()] = used_sim_tracks[trk_idx];
        continue;
      }
      // Perform the actual vertex collapsing, if needed.
      if (collapsed_vertices[origin_vtx])
        origin_vtx = collapsed_vertices[origin_vtx];
      add_edge(origin_vtx,
               v.vertexId(),
               EdgeProperty(&tracks[trk_idx], simTrackDetIdEnergyMap[v.parentIndex()].size(), 0),
               decay);
      used_sim_tracks[trk_idx] = v.vertexId();
    }
  }
  // Build the motherParticle property to each vertex
  auto const &vertexMothersProp = get(vertex_name, decay);
  // Now recover the particles that did not decay. Append them with an index
  // bigger than the size of the generated vertices.
  int offset = vertices.size();
  for (size_t i = 0; i < tracks.size(); ++i) {
    if (!used_sim_tracks[i]) {
      auto origin_vtx = tracks[i].vertIndex();
      // Perform the actual vertex collapsing, if needed.
      if (collapsed_vertices[origin_vtx])
        origin_vtx = collapsed_vertices[origin_vtx];
      add_edge(
          origin_vtx, offset, EdgeProperty(&tracks[i], simTrackDetIdEnergyMap[tracks[i].trackId()].size(), 0), decay);
      // The properties for "fake" vertices associated to stable particles have
      // to be set inside this loop, since they do not belong to the vertices
      // collection and would be skipped by that loop (coming next)
      put(vertexMothersProp, offset, VertexProperty(&tracks[i], 0));
      offset++;
    }
  }
  for (auto const &v : vertices) {
    if (v.parentIndex() != -1) {
      // Skip collapsed_vertices
      if (collapsed_vertices[v.vertexId()])
        continue;
      put(vertexMothersProp, v.vertexId(), VertexProperty(&tracks[trackid_to_track_index[v.parentIndex()]], 0));
    }
  }
  SimHitsAccumulator_dfs_visitor vis;
  depth_first_search(decay, visitor(vis));

  VisitorHelper visitorHelper(m_simHitBarcodeToIndex, simTrackDetIdEnergyMap, vertex_time_map);

  auto passGenPartSelections_lambda = [&](EdgeProperty const &edge_property) {
    return (edge_property.cumulative_simHits != 0 and !edge_property.simTrack->noGenpart() and
            edge_property.simTrack->momentum().E() > minEnergy_ and
            std::abs(edge_property.simTrack->momentum().Eta()) < maxPseudoRapidity_);
  };

  // Do the graph search for all 3 vistors at the same tim
  PrimaryVisitor primaryVisitor(
      visitorHelper,
      outputCaloParticles_,
      passGenPartSelections_lambda,
      std::make_tuple(
          // std::vector<SubClusterT>& clusters, ClusterParentIndexRecorder<ParentClusterCollectionT> parentIndexRecorder, VisitorHelper const& helper, Selector_t selector
          SubClusterVisitor(
            caloParticleSimClusters_config_.outputClusters,
            ClusterParentIndexRecorder(caloParticleSimClusters_config_.clustersToCaloParticleMap, caloParticles_refProd_), // Trivial 1-1 mapping, kept for convenience
            //ClusterParentIndexRecorderVoid(), // We don't need to make references to self
            visitorHelper,
            [&](EdgeProperty const &edge_property) -> bool {
              // Create a SimCluster for every CaloParticle (duplicates the CaloParticle for convenience of use, to have one single dataformat)
              return true;
            }),
          LegacySimClusterVisitor(
              legacySimClusters_config_.outputClusters,
              ClusterParentIndexRecorder(legacySimClusters_config_.clustersToCaloParticleMap, caloParticles_refProd_),
              visitorHelper
          ),
          SubClusterVisitor(
              boundarySimClusters_config_.outputClusters,
              ClusterParentIndexRecorder(boundarySimClusters_config_.clustersToCaloParticleMap, caloParticles_refProd_),
              visitorHelper,
              [&](EdgeProperty const &edge_property) -> bool {
                // Create SimCluster from every SimTrack crossing boundary (and which has simhits either itself as in its descendants), and that is inside a CaloParticle
                return edge_property.cumulative_simHits != 0 && edge_property.simTrack->crossedBoundary();
              })
      )
  );
  depth_first_search(decay, visitor(primaryVisitor));  //  "visitor()" is a Boost BGL named parameter

#if DEBUG
  boost::write_graphviz(std::cout,
                        decay,
                        make_label_writer(make_transform_value_property_map(&graphviz_vertex, get(vertex_name, decay))),
                        make_label_writer(make_transform_value_property_map(&graphviz_edge, get(edge_weight, decay))));
#endif
}

template <class T>
void CaloTruthAccumulator::fillSimHits(std::vector<std::pair<DetId, const PCaloHit *>> &returnValue,
                                       std::unordered_map<int, std::map<int, float>> &simTrackDetIdEnergyMap,
                                       const T &event,
                                       const edm::EventSetup &setup) {
  for (auto const &collectionTag : collectionTags_) {
    edm::Handle<std::vector<PCaloHit>> hSimHits;
    const bool isHcal = (collectionTag.instance().find("HcalHits") != std::string::npos);
    const bool isHGCal = (collectionTag.instance().find("HGCHits") != std::string::npos);
    event.getByLabel(collectionTag, hSimHits);

    for (auto const &simHit : *hSimHits) {
      DetId id(0);

      //Relabel as necessary for HGCAL
      if (isHGCal) {
        const uint32_t simId = simHit.id();
        if (geometryType_ == 1) {
          // no test numbering in new geometry
          id = simId;
        } else if (isHcal) {
          HcalDetId hid = HcalHitRelabeller::relabel(simId, hcddd_);
          if (hid.subdet() == HcalEndcap)
            id = hid;
        } else {
          int subdet, layer, cell, sec, subsec, zp;
          HGCalTestNumbering::unpackHexagonIndex(simId, subdet, zp, layer, sec, subsec, cell);
          const HGCalDDDConstants *ddd = hgddd_[subdet - 3];
          std::pair<int, int> recoLayerCell = ddd->simToReco(cell, layer, sec, hgtopo_[subdet - 3]->detectorType());
          cell = recoLayerCell.first;
          layer = recoLayerCell.second;
          // skip simhits with bad barcodes or non-existant layers
          if (layer == -1 || simHit.geantTrackId() == 0)
            continue;
          id = HGCalDetId((ForwardSubdetector)subdet, zp, layer, subsec, sec, cell);
        }
      } else {
        id = simHit.id();
        //Relabel all HCAL hits
        if (isHcal) {
          HcalDetId hid = HcalHitRelabeller::relabel(simHit.id(), hcddd_);
          id = hid;
        }
      }

      if (id == DetId(0)) {
        continue;
      }
      if (simHit.geantTrackId() == 0) {
        continue;
      }

      returnValue.emplace_back(id, &simHit);
      simTrackDetIdEnergyMap[simHit.geantTrackId()][id.rawId()] += simHit.energy();
      m_detIdToTotalSimEnergy[id.rawId()] += simHit.energy();
    }
  }  // end of loop over InputTags
}

// Register with the framework
DEFINE_DIGI_ACCUMULATOR(CaloTruthAccumulator);
