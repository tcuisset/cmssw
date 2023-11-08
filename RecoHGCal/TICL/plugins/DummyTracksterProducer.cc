// -*- C++ -*-
//
// Package:    UserCode/DummyTracksterProducer
// Class:      DummyTracksterProducer
//
/**\class DummyTracksterProducer DummyTracksterProducer.cc UserCode/DummyTracksterProducer/plugins/DummyTracksterProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Theo Cuisset
//         Created:  Fri, 20 Oct 2023 13:25:33 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"


//
// class declaration
//

class DummyTracksterProducer : public edm::stream::EDProducer<> {
public:
  explicit DummyTracksterProducer(const edm::ParameterSet&);
  ~DummyTracksterProducer() override {};

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    void beginJob();
  void endJob();

  void beginRun(edm::Run const &iEvent, edm::EventSetup const &es) override;

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // ----------member data ---------------------------
  int dummyTracksterCount;
};
void DummyTracksterProducer::beginJob() {}

void DummyTracksterProducer::endJob(){};

void DummyTracksterProducer::beginRun(edm::Run const &iEvent, edm::EventSetup const &es) {
};

DummyTracksterProducer::DummyTracksterProducer(const edm::ParameterSet& iConfig) 
  : dummyTracksterCount(iConfig.getParameter<int>("dummyTracksterCount")) {
  produces<std::vector<ticl::Trackster>>();
}


// ------------ method called to produce the data  ------------
void DummyTracksterProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace ticl;
  Trackster defaultTrackster;
  defaultTrackster.setBarycenter(ticl::Trackster::Vector(1, 2, 3));
  defaultTrackster.setRawEnergy(5);
  Eigen::Vector3d eigenvals({0.5, 0.2, 1.2});
  Eigen::Matrix3d eigenvectors({{0., 0.5, 0.2}, {0.1, 0.2, 0.3}, {0.9, 0.2, 0.5}});
  Eigen::Vector3d sigmas({1., 2., 3.});
  defaultTrackster.fillPCAVariables(eigenvals, eigenvectors, sigmas, sigmas, 
       3, ticl::Trackster::PCAOrdering::ascending);
  auto tracksters = std::make_unique<std::vector<Trackster>>(dummyTracksterCount, defaultTrackster);
  
  iEvent.put(std::move(tracksters));
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DummyTracksterProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<int>("dummyTracksterCount", 3000);
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DummyTracksterProducer);
