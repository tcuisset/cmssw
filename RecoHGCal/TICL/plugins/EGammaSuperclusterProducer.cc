// Author : Theo Cuisset - theo.cuisset@polytechnique.edu
// Date : 01/2024
/* Translates TICL supercluster to ECAL supercluster, same as RecoEcal/EgammaCLusterProducers/PFECALSuperClusterProducer */

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"

class EGammaSuperclusterProducer : public edm::stream::EDProducer<> {
public:
  EGammaSuperclusterProducer(const edm::ParameterSet&);

  void produce(edm::Event&, const edm::EventSetup&) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    edm::EDGetTokenT<ticl::TracksterCollection> ticlSuperClustersToken_;
    edm::EDGetTokenT<std::vector<std::vector<unsigned int>>> superClusterLinksToken_;
    edm::EDGetTokenT<ticl::TracksterCollection> ticlTrackstersEMToken_;
    edm::EDGetTokenT<reco::CaloClusterCollection> layerClustersToken_;
};


EGammaSuperclusterProducer::EGammaSuperclusterProducer(const edm::ParameterSet& ps)
    : ticlSuperClustersToken_(consumes<ticl::TracksterCollection>(ps.getParameter<edm::InputTag>("ticlSuperClusters"))),
    superClusterLinksToken_(consumes<std::vector<std::vector<unsigned int>>>(edm::InputTag(ps.getParameter<edm::InputTag>("ticlSuperClusters").label(), "linkedTracksterIdToInputTracksterId", ps.getParameter<edm::InputTag>("ticlSuperClusters").process()))),
    ticlTrackstersEMToken_(consumes<ticl::TracksterCollection>(ps.getParameter<edm::InputTag>("ticlTrackstersEM"))),
    layerClustersToken_(consumes<reco::CaloClusterCollection>(ps.getParameter<edm::InputTag>("layerClusters"))) {
  produces<reco::SuperClusterCollection>();
  produces<reco::CaloClusterCollection>(); // The CaloCluster corresponding to each EM trackster 
}


void EGammaSuperclusterProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto const& ticlSuperclusters = iEvent.get(ticlSuperClustersToken_);
  auto const& ticlSuperclusterLinks = iEvent.get(superClusterLinksToken_);
  auto emTracksters_h = iEvent.getHandle(ticlTrackstersEMToken_);
  auto const& emTracksters = *emTracksters_h;
  auto const& layerClusters = iEvent.get(layerClustersToken_);
  auto egammaSuperclusters = std::make_unique<reco::SuperClusterCollection>();
  auto caloClustersEM = std::make_unique<reco::CaloClusterCollection>();

  // Fill CaloCluster collection
  for (ticl::Trackster const& emTrackster : emTracksters) {
    std::vector<std::pair<DetId, float> > hitsAndFractions;
    int iLC = 0;
    std::for_each(std::begin(emTrackster.vertices()), std::end(emTrackster.vertices()), [&](unsigned int lcId) {
      const auto fraction = 1.f / emTrackster.vertex_multiplicity(iLC++);
      for (const auto& cell : layerClusters[lcId].hitsAndFractions()) {
        hitsAndFractions.emplace_back(cell.first, cell.second * fraction);
      }
    });

    reco::CaloCluster& caloCluster = caloClustersEM->emplace_back(
        emTrackster.regressed_energy(), // energy
        math::XYZPoint(emTrackster.barycenter()), // position
        reco::CaloID(reco::CaloID::DET_HGCAL_ENDCAP),
        hitsAndFractions,
        reco::CaloCluster::particleFlow, // algoID (copying from output of PFECALCSuperClusterProducer)
        hitsAndFractions.at(0).first // seedId : TODO figure out if needed or if we need to put the highest energy hit
    );
    caloCluster.setCorrectedEnergy(emTrackster.regressed_energy());

  }

  edm::OrphanHandle<reco::CaloClusterCollection> caloClustersEM_h = iEvent.put(std::move(caloClustersEM));

  assert(ticlSuperclusters.size() == ticlSuperclusterLinks.size());
  for (std::size_t sc_i = 0; sc_i < ticlSuperclusters.size(); sc_i++) {
    ticl::Trackster const& ticlSupercluster = ticlSuperclusters[sc_i];
    std::vector<unsigned int> const& superclusterLink = ticlSuperclusterLinks[sc_i];

    //if (superclusterLink.size() == 1)
    //  continue; // Drop superclusters that are formed of one trackster only

    reco::CaloClusterPtrVector trackstersEMInSupercluster; 
    double regressedEnergySum = 0.; // Sum of regressed_energy of all tracksters in supercluster
    for (unsigned int tsInSc_id : superclusterLink) {
      trackstersEMInSupercluster.push_back(reco::CaloClusterPtr(caloClustersEM_h, tsInSc_id));
      regressedEnergySum += emTracksters[tsInSc_id].regressed_energy();
    }
    reco::SuperCluster& egammaSc = egammaSuperclusters->emplace_back(
      ticlSupercluster.raw_energy(),
      reco::SuperCluster::Point(ticlSupercluster.barycenter()),
      reco::CaloClusterPtr(caloClustersEM_h, superclusterLink[0]), // seed (first trackster in superclusterLink is the seed)
      trackstersEMInSupercluster, // clusters
      0., // Epreshower
      0.046, // phiwidth (TODO placheolder value for now)
      0.017 // etawidth (TODO placheolder value for now)
    );
    egammaSc.setCorrectedEnergy(regressedEnergySum); // TODO 
    //egammaSc.setCorrectedEnergyUncertainty(0.01 * regressedEnergySum); // TODO placeholder value

  }
  
  iEvent.put(std::move(egammaSuperclusters));
}

void EGammaSuperclusterProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("ticlSuperClusters", edm::InputTag("ticlTracksterLinksSuperclustering"));
  desc.add<edm::InputTag>("ticlTrackstersEM", edm::InputTag("ticlTrackstersCLUE3DEM"))->setComment("The trackster collection used before superclustering, ie CLUE3D EM tracksters");
  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalMergeLayerClusters"))->setComment("The layer cluster collection that goes with ticlTrackstersEM");

  descriptions.add("ticlEGammaSuperClusterProducer", desc);
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EGammaSuperclusterProducer);
