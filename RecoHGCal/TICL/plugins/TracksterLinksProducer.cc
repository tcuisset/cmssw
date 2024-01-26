// Author: Felice Pantaleo, Wahid Redjeb (CERN) - felice.pantaleo@cern.ch, wahid.redjeb@cern.ch
// Date: 12/2023
#include <memory>  // unique_ptr
#include "CommonTools/RecoAlgos/interface/MultiVectorManager.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

#include "DataFormats/Common/interface/OrphanHandle.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/Common.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"

#include "RecoHGCal/TICL/interface/TracksterLinkingAlgoBase.h"
#include "RecoHGCal/TICL/plugins/TracksterLinkingPluginFactory.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "TrackstersPCA.h"

using namespace ticl;
using cms::Ort::ONNXRuntime;

class TracksterLinksProducer : public edm::stream::EDProducer<edm::GlobalCache<ONNXRuntime>> {
public:
  explicit TracksterLinksProducer(const edm::ParameterSet &ps, const ONNXRuntime *);
  ~TracksterLinksProducer() override{};
  void produce(edm::Event &, const edm::EventSetup &) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  void beginJob();
  void endJob();

  void beginRun(edm::Run const &iEvent, edm::EventSetup const &es) override;
  static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet &iConfig);
  static void globalEndJob(const ONNXRuntime *);

private:
  void printTrackstersDebug(const std::vector<Trackster> &, const char *label) const;
  void dumpTrackster(const Trackster &) const;

  std::unique_ptr<TracksterLinkingAlgoBase> linkingAlgo_;
  std::string algoType_;

  std::vector<edm::EDGetTokenT<std::vector<Trackster>>> tracksters_tokens_;
  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_token_;
  const edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> clustersTime_token_;

  std::vector<edm::EDGetTokenT<std::vector<float>>> original_masks_tokens_;

  std::once_flag initializeGeometry_;

  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometry_token_;
  const std::string detector_;
  const std::string propName_;

  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bfield_token_;
  const edm::ESGetToken<Propagator, TrackingComponentsRecord> propagator_token_;
  const HGCalDDDConstants *hgcons_;
  hgcal::RecHitTools rhtools_;
  edm::ESGetToken<HGCalDDDConstants, IdealGeometryRecord> hdc_token_;
};

TracksterLinksProducer::TracksterLinksProducer(const edm::ParameterSet &ps, const ONNXRuntime *onnxRuntime)
    : algoType_(ps.getParameter<edm::ParameterSet>("linkingPSet").getParameter<std::string>("type")),
      clusters_token_(consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layer_clusters"))),
      clustersTime_token_(
          consumes<edm::ValueMap<std::pair<float, float>>>(ps.getParameter<edm::InputTag>("layer_clustersTime"))),
      geometry_token_(esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>()),
      detector_(ps.getParameter<std::string>("detector")),
      propName_(ps.getParameter<std::string>("propagator")),
      bfield_token_(esConsumes<MagneticField, IdealMagneticFieldRecord, edm::Transition::BeginRun>()),
      propagator_token_(
          esConsumes<Propagator, TrackingComponentsRecord, edm::Transition::BeginRun>(edm::ESInputTag("", propName_))) {
  // Loop over the edm::VInputTag and append the token to tracksters_tokens_
  for (auto const &tag : ps.getParameter<std::vector<edm::InputTag>>("tracksters_collections")) {
    tracksters_tokens_.emplace_back(consumes<std::vector<Trackster>>(tag));
  }
  //Loop over the edm::VInputTag of masks and append the token to original_masks_tokens_
  for (auto const &tag : ps.getParameter<std::vector<edm::InputTag>>("original_masks")) {
    original_masks_tokens_.emplace_back(consumes<std::vector<float>>(tag));
  }

  // New trackster collection after linking
  produces<std::vector<Trackster>>();
  produces<std::vector<Trackster>>("CLUE3D");

  // Links
  produces<std::vector<std::vector<unsigned int>>>();
  // LayerClusters Mask
  produces<std::vector<float>>();

  auto linkingPSet = ps.getParameter<edm::ParameterSet>("linkingPSet");

  if (algoType_ == "Skeletons") {
    std::string detectorName_ = (detector_ == "HFNose") ? "HGCalHFNoseSensitive" : "HGCalEESensitive";
    hdc_token_ = esConsumes<HGCalDDDConstants, IdealGeometryRecord, edm::Transition::BeginRun>(
        edm::ESInputTag("", detectorName_));
  }

  linkingAlgo_ = TracksterLinkingPluginFactory::get()->create(algoType_, linkingPSet, consumesCollector(), onnxRuntime);
}

std::unique_ptr<ONNXRuntime> TracksterLinksProducer::initializeGlobalCache(const edm::ParameterSet &iConfig) {
  auto const& pluginPset = iConfig.getParameter<edm::ParameterSet>("linkingPSet");
  if (pluginPset.exists("onnxModelPath"))
    return std::make_unique<ONNXRuntime>(pluginPset.getParameter<edm::FileInPath>("onnxModelPath").fullPath());
  else
    return std::unique_ptr<ONNXRuntime>(nullptr);
}

void TracksterLinksProducer::beginJob() {}

void TracksterLinksProducer::endJob(){};

void TracksterLinksProducer::globalEndJob(const ONNXRuntime *cache) {}

void TracksterLinksProducer::beginRun(edm::Run const &iEvent, edm::EventSetup const &es) {
  if (algoType_ == "Skeletons") {
    edm::ESHandle<HGCalDDDConstants> hdc = es.getHandle(hdc_token_);
    hgcons_ = hdc.product();
  }

  edm::ESHandle<CaloGeometry> geom = es.getHandle(geometry_token_);
  rhtools_.setGeometry(*geom);

  edm::ESHandle<MagneticField> bfield = es.getHandle(bfield_token_);
  edm::ESHandle<Propagator> propagator = es.getHandle(propagator_token_);

  linkingAlgo_->initialize(hgcons_, rhtools_, bfield, propagator);
};

void TracksterLinksProducer::dumpTrackster(const Trackster &t) const {
  auto e_over_h = (t.raw_em_pt() / ((t.raw_pt() - t.raw_em_pt()) != 0. ? (t.raw_pt() - t.raw_em_pt()) : 1.));
  LogDebug("TracksterLinksProducer")
      << "\nTrackster raw_pt: " << t.raw_pt() << " raw_em_pt: " << t.raw_em_pt() << " eoh: " << e_over_h
      << " barycenter: " << t.barycenter() << " eta,phi (baricenter): " << t.barycenter().eta() << ", "
      << t.barycenter().phi() << " eta,phi (eigen): " << t.eigenvectors(0).eta() << ", " << t.eigenvectors(0).phi()
      << " pt(eigen): " << std::sqrt(t.eigenvectors(0).Unit().perp2()) * t.raw_energy() << " seedID: " << t.seedID()
      << " seedIndex: " << t.seedIndex() << " size: " << t.vertices().size() << " average usage: "
      << (std::accumulate(std::begin(t.vertex_multiplicity()), std::end(t.vertex_multiplicity()), 0.) /
          (float)t.vertex_multiplicity().size())
      << " raw_energy: " << t.raw_energy() << " regressed energy: " << t.regressed_energy()
      << " probs(ga/e/mu/np/cp/nh/am/unk): ";
  for (auto const &p : t.id_probabilities()) {
    LogDebug("TracksterLinksProducer") << std::fixed << p << " ";
  }
  LogDebug("TracksterLinksProducer") << " sigmas: ";
  for (auto const &s : t.sigmas()) {
    LogDebug("TracksterLinksProducer") << s << " ";
  }
  LogDebug("TracksterLinksProducer") << std::endl;
}

void TracksterLinksProducer::produce(edm::Event &evt, const edm::EventSetup &es) {
  auto resultTracksters = std::make_unique<std::vector<Trackster>>();

  auto linkedResultTracksters = std::make_unique<std::vector<std::vector<unsigned int>>>();

  const auto &layerClusters = evt.get(clusters_token_);
  const auto &layerClustersTimes = evt.get(clustersTime_token_);

  // loop over the original_masks_tokens_ and get the original masks collections and multiply them
  // to get the global mask
  std::vector<float> original_global_mask(layerClusters.size(), 1.f);
  for (unsigned int i = 0; i < original_masks_tokens_.size(); ++i) {
    const auto &tmp_mask = evt.get(original_masks_tokens_[i]);
    for (unsigned int j = 0; j < tmp_mask.size(); ++j) {
      original_global_mask[j] *= tmp_mask[j];
    }
  }

  auto resultMask = std::make_unique<std::vector<float>>(original_global_mask);

  std::vector<edm::Handle<std::vector<Trackster>>> tracksters_h(tracksters_tokens_.size());
  MultiVectorManager<Trackster> trackstersManager;
  for (unsigned int i = 0; i < tracksters_tokens_.size(); ++i) {
    evt.getByToken(tracksters_tokens_[i], tracksters_h[i]);
    //Fill MultiVectorManager
    trackstersManager.addVector(*tracksters_h[i]);
  }

  // Linking
  const typename TracksterLinkingAlgoBase::Inputs input(evt, es, layerClusters, layerClustersTimes, trackstersManager);
  std::vector<std::vector<unsigned int>> linkedTracksterIdToInputTracksterId;

  // LinkTracksters will produce a vector of vector of indices of tracksters that:
  // 1) are linked together if more than one
  // 2) are isolated if only one
  // Result tracksters contains the final version of the trackster collection
  // linkedTrackstersToInputTrackstersMap contains the mapping between the linked tracksters and the input tracksters
  linkingAlgo_->linkTracksters(input, *resultTracksters, *linkedResultTracksters, linkedTracksterIdToInputTracksterId);

  // Now we need to remove the tracksters that are not linked
  // We need to emplace_back in the resultTracksters only the tracksters that are linked

  for (auto const &resultTrackster : *resultTracksters) {
    for (auto const &clusterIndex : resultTrackster.vertices()) {
      (*resultMask)[clusterIndex] = 0.f;
    }
  }

  assignPCAtoTracksters(
      *resultTracksters, layerClusters, layerClustersTimes, rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z(), true);

  evt.put(std::move(linkedResultTracksters));
  evt.put(std::move(resultMask));
  evt.put(std::move(resultTracksters));
}

void TracksterLinksProducer::printTrackstersDebug(const std::vector<Trackster> &tracksters, const char *label) const {
#ifdef EDM_ML_DEBUG
  int counter = 0;
  for (auto const &t : tracksters) {
    LogDebug("TracksterLinksProducer")
        << counter++ << " TracksterLinksProducer (" << label << ") obj barycenter: " << t.barycenter()
        << " eta,phi (baricenter): " << t.barycenter().eta() << ", " << t.barycenter().phi()
        << " eta,phi (eigen): " << t.eigenvectors(0).eta() << ", " << t.eigenvectors(0).phi()
        << " pt(eigen): " << std::sqrt(t.eigenvectors(0).Unit().perp2()) * t.raw_energy() << " seedID: " << t.seedID()
        << " seedIndex: " << t.seedIndex() << " size: " << t.vertices().size() << " average usage: "
        << (std::accumulate(std::begin(t.vertex_multiplicity()), std::end(t.vertex_multiplicity()), 0.) /
            (float)t.vertex_multiplicity().size())
        << " raw_energy: " << t.raw_energy() << " regressed energy: " << t.regressed_energy()
        << " probs(ga/e/mu/np/cp/nh/am/unk): ";
    for (auto const &p : t.id_probabilities()) {
      LogDebug("TracksterLinksProducer") << std::fixed << p << " ";
    }
    LogDebug("TracksterLinksProducer") << " sigmas: ";
    for (auto const &s : t.sigmas()) {
      LogDebug("TracksterLinksProducer") << s << " ";
    }
    LogDebug("TracksterLinksProducer") << std::endl;
  }
#endif
}

void TracksterLinksProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  edm::ParameterSetDescription linkingDesc;
  linkingDesc.addNode(edm::PluginDescription<TracksterLinkingPluginFactory>("type", "Skeletons", true));

  desc.add<edm::ParameterSetDescription>("linkingPSet", linkingDesc);
  desc.add<std::vector<edm::InputTag>>(
      "tracksters_collections", {edm::InputTag("ticlTrackstersCLUE3DEM"), edm::InputTag("ticlTrackstersCLUE3DHAD")});
  desc.add<std::vector<edm::InputTag>>("original_masks",
                                       {edm::InputTag("hgcalMergeLayerClusters", "InitialLayerClustersMask")});
  desc.add<edm::InputTag>("layer_clusters", edm::InputTag("hgcalMergeLayerClusters"));
  desc.add<edm::InputTag>("layer_clustersTime", edm::InputTag("hgcalMergeLayerClusters", "timeLayerCluster"));
  desc.add<std::string>("detector", "HGCAL");
  desc.add<std::string>("propagator", "PropagatorWithMaterial");
  descriptions.add("tracksterLinksProducer", desc);
}

DEFINE_FWK_MODULE(TracksterLinksProducer);
