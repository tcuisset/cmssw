// Author: Felice Pantaleo, felice.pantaleo@cern.ch 06/2024
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"
#include "DataFormats/Provenance/interface/ProductID.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/Common/interface/RefProdVector.h"
#include "DataFormats/Common/interface/MultiSpan.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

class HitToSimClusterCaloParticleAssociatorProducer : public edm::global::EDProducer<> {
public:
  explicit HitToSimClusterCaloParticleAssociatorProducer(const edm::ParameterSet &);
  ~HitToSimClusterCaloParticleAssociatorProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;

  const edm::EDGetTokenT<std::vector<SimCluster>> simClusterToken_;
  const edm::EDGetTokenT<std::vector<CaloParticle>> caloParticleToken_;

  const edm::EDGetTokenT<std::unordered_map<DetId, const unsigned int>> hitMapToken_;
  edm::EDGetTokenT<edm::RefProdVector<HGCRecHitCollection>> hitsToken_;
};


HitToSimClusterCaloParticleAssociatorProducer::HitToSimClusterCaloParticleAssociatorProducer(
    const edm::ParameterSet &pset)
    : simClusterToken_(consumes<std::vector<SimCluster>>(pset.getParameter<edm::InputTag>("simClusters"))),
      hitMapToken_(consumes<std::unordered_map<DetId, const unsigned int>>(pset.getParameter<edm::InputTag>("hitMap"))),
      hitsToken_(consumes<edm::RefProdVector<HGCRecHitCollection>>(pset.getParameter<edm::InputTag>("hits"))) {
  produces<ticl::AssociationMap<ticl::mapWithFraction>>();
}

void HitToSimClusterCaloParticleAssociatorProducer::produce(edm::StreamID,
                                                            edm::Event &iEvent,
                                                            const edm::EventSetup &iSetup) const {
  using namespace edm;

  std::vector<SimCluster> simClusters = iEvent.get(simClusterToken_);
  Handle<std::unordered_map<DetId, const unsigned int>> hitMap;
  iEvent.getByToken(hitMapToken_, hitMap);

  if (!iEvent.getHandle(hitsToken_).isValid()) {
    edm::LogWarning("HitToSimClusterCaloParticleAssociatorProducer")
        << "No valid HGCRecHitCollections found. Association maps will be empty.";
    // Store empty maps in the event
    iEvent.put(std::make_unique<ticl::AssociationMap<ticl::mapWithFraction>>());
    return;
  }

  // Protection against missing HGCRecHitCollection
  const auto hits = iEvent.get(hitsToken_);
  for (std::size_t index = 0; const auto &hgcRecHitCollection : hits) {
    if (hgcRecHitCollection->empty()) {
      edm::LogWarning("HitToSimClusterCaloParticleAssociatorProducer")
          << "HGCRecHitCollection #" << index << " is empty or not valid.";
    }
    index++;
  }

  edm::MultiSpan<HGCRecHit> rechitSpan(hits);
  // Check if rechitSpan is empty after processing hitsTokens_
  if (rechitSpan.size() == 0) {
    edm::LogWarning("HitToSimClusterCaloParticleAssociatorProducer")
        << "No valid HGCRecHitCollections found. Association maps will be empty.";
    // Store empty maps in the event
    iEvent.put(std::make_unique<ticl::AssociationMap<ticl::mapWithFraction>>());
    return;
  }

  // Create association maps
  auto hitToSimClusterMap = std::make_unique<ticl::AssociationMap<ticl::mapWithFraction>>(rechitSpan.size());

  for (std::size_t scId = 0; scId < simClusters.size(); ++scId) {
    // Loop over hits in simCluster
    for (const auto &hitAndFraction : simClusters[scId].hits_and_fractions()) {
      auto hitMapIter = hitMap->find(hitAndFraction.first);
      if (hitMapIter != hitMap->end()) {
        unsigned int rechitIndex = hitMapIter->second;
        float fraction = hitAndFraction.second;
        hitToSimClusterMap->insert(rechitIndex, scId, fraction);
      }
    }
  }
  iEvent.put(std::move(hitToSimClusterMap));
}

void HitToSimClusterCaloParticleAssociatorProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("simClusters", edm::InputTag("mix", "MergedCaloTruth"));

  desc.add<edm::InputTag>("hitMap", edm::InputTag("recHitMapProducer", "hgcalRecHitMap"));
  desc.add<edm::InputTag>("hits", edm::InputTag("recHitMapProducer", "RefProdVectorHGCRecHitCollection"));
  descriptions.add("hitToSimClusterCaloParticleAssociator", desc);
}

// Define this as a plug-in
DEFINE_FWK_MODULE(HitToSimClusterCaloParticleAssociatorProducer);
