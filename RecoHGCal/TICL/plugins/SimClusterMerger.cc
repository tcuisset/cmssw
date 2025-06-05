#include <numeric>
#include <memory>  // unique_ptr
#include <ranges>

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

#include "DataFormats/Common/interface/OrphanHandle.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"
#include "DataFormats/HGCalReco/interface/Common.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include <Math/PositionVector3D.h>
#include <Math/GenVector/CoordinateSystemTags.h>
#include <fastjet/ClusterSequence.hh>

#include "TrackstersPCA.h"



using namespace ticl;

class SimClusterMerger : public edm::stream::EDProducer<> {
public:
  explicit SimClusterMerger(const edm::ParameterSet &ps);
  ~SimClusterMerger() override {};
  void produce(edm::Event &, const edm::EventSetup &) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  const edm::EDGetTokenT<std::vector<SimCluster>> simclusters_token_;
  const edm::EDGetTokenT<std::vector<CaloParticle>> caloparticles_token_;
  double jetClusteringRadius_;
};

SimClusterMerger::SimClusterMerger(const edm::ParameterSet &ps)
    : simclusters_token_(consumes<std::vector<SimCluster>>(ps.getParameter<edm::InputTag>("simclusters"))),
      caloparticles_token_(consumes(ps.getParameter<edm::InputTag>("caloparticles"))),
      jetClusteringRadius_(ps.getParameter<double>("jetClusteringRadius")) {
    
    produces<SimClusterCollection>();
    produces<ticl::AssociationMap<ticl::oneToOneMapWithFraction, std::vector<SimCluster>, std::vector<CaloParticle>>>(
        "simClusterToCaloParticleMap");
}

void SimClusterMerger::produce(edm::Event &evt, const edm::EventSetup &es) {
edm::Handle<std::vector<SimCluster>> simClusters_h;
  evt.getByToken(simclusters_token_, simClusters_h);
  const std::vector<SimCluster>& simClusters = *simClusters_h;

  edm::Handle<std::vector<CaloParticle>> caloparticles_h;
  evt.getByToken(caloparticles_token_, caloparticles_h);
  const std::vector<CaloParticle>& caloparticles = *caloparticles_h;

  auto mergedSimClusters = std::make_unique<SimClusterCollection>();
  auto simClusterToCaloParticleMap = std::make_unique<
      ticl::AssociationMap<ticl::oneToOneMapWithFraction, std::vector<SimCluster>, std::vector<CaloParticle>>>(
      simClusters_h, caloparticles_h, evt);
  //auto simTrackstersMergedToSimTrackstersMap = std::make_unique<std::vector<std::vector<edm::Ref<Trackster>>>>;
//   auto simTrackstersMergedToSimTrackstersMap = std::make_unique<edm::AssociationMap<edm::OneToMany<Trackster, Trackster>>>;
//   auto resultTracksterLinks = std::make_unique<std::vector<std::vector<unsigned int>>>(); // List of components


  for (std::size_t i_cp = 0; i_cp < caloparticles.size(); i_cp++) {
    CaloParticle const& cp = caloparticles[i_cp];

    if (cp.g4Tracks()[0].crossedBoundary()) {
      // In case the CaloParticle crossed boundary, we just merge all SimClusters
      assert(!cp.simClusters().empty());
      mergedSimClusters->emplace_back(std::ranges::subrange(cp.simClusters().begin(), cp.simClusters().end()) | std::views::transform([](auto const& ref){return *ref;}));
      // mergedSimClusters->emplace_back()
      simClusterToCaloParticleMap->insert(mergedSimClusters->size()-1, i_cp, 1.); // fraction of the simCluster energy to the caloParticle energy is 1
      continue;
    }

    std::vector<fastjet::PseudoJet> fjInputs;
    for (auto const& sc : cp.simClusters()) {
        SimTrack const& simtrack = sc->g4Tracks()[0];
        if (simtrack.crossedBoundary()) {
            // positionsAtBoundary.emplace_back(simtrack.getPositionAtBoundary().x(), simtrack.getPositionAtBoundary().y());
            // momentumsAtBoundary.emplace_back(simtrack.getMomentumAtBoundary().P());

            /* Build the particle 4-vector such that the energy is the energy of SimTrack at boundary,
            the momentum 3-vector points to the boundary position, and the mass is zero */
            auto energyAtBoundary = simtrack.getMomentumAtBoundary().E();
            auto momentum3D = energyAtBoundary * simtrack.getPositionAtBoundary().Vect().Unit();
            auto& jet = fjInputs.emplace_back(momentum3D.X(), momentum3D.Y(), momentum3D.Z(), energyAtBoundary);
            jet.set_user_index(sc.index()); // Store the index in SimCLuster collection for later
        }
    }

    // Clustering
    fastjet::ClusterSequence sequence(fjInputs, fastjet::JetDefinition(fastjet::antikt_algorithm, jetClusteringRadius_));
    auto jets = sequence.inclusive_jets();

    // Merging
    // if (jets.size() == 1) {
    //   // Merge the SimTracksterSC into one (actually just take the SimTracksterCP)
    //   resultTracksters->push_back(simTs_CP);
    // } else {
      
      for (fastjet::PseudoJet const& jet : jets) {
        auto constituents = fastjet::sorted_by_E(jet.constituents());
        assert(constituents.size() >= 1);
        

        // Build the merged trackster by cloning the leading energy SimTrackster from the cluster
        // Then add into it all the other SimTrackster of the cluster

        mergedSimClusters->emplace_back(constituents | std::views::transform([&](fastjet::PseudoJet const& pseudoJet) { return simClusters[pseudoJet.user_index()];}));
        simClusterToCaloParticleMap->insert(mergedSimClusters->size()-1, i_cp, std::accumulate(constituents.begin(), constituents.end(), 0.,
        [&](double fraction, fastjet::PseudoJet const& pseudoJet)
        {
            return fraction + simClusters[pseudoJet.user_index()].energy();
        })/ cp.energy() );

        // TODO set pdgId to something
    }

  }

 

  evt.put(std::move(mergedSimClusters));
  evt.put(std::move(simClusterToCaloParticleMap));
  }

void SimClusterMerger::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("simclusters", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("caloparticles", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<double>("jetClusteringRadius", 0.05)->setComment("Distance parameter for clustering algorithm");
  descriptions.add("simClusterMerger", desc);
}

DEFINE_FWK_MODULE(SimClusterMerger);
