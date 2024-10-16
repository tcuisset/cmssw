#include <ostream>

#include "IOMC/ParticleGuns/interface/MultiParticleProducer.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "CLHEP/Random/RandFlat.h"

namespace CLHEP {
  class HepRandomEngine;
}

using namespace edm;

MultiParticleProducer::MultiParticleProducer(const edm::ParameterSet& pset)
    : BaseFlatGunProducer(pset) {
  edm::ParameterSet defpset;
  edm::ParameterSet pgun_params = pset.getParameter<edm::ParameterSet>("PGunParameters");

  // doesn't seem necessary to check if pset is empty - if this
  // is the case, default values will be taken for params
  fMinE = pgun_params.getParameter<double>("MinE");
  fMaxE = pgun_params.getParameter<double>("MaxE");

  fNParticlesInEta = pgun_params.getParameter<unsigned int>("NParticlesInEta");
  fEtaMinSeparation = pgun_params.getParameter<double>("EtaMinSeparation");
  fEtaRandomSpread = pgun_params.getParameter<double>("EtaRandomSpread");
  fNParticlesInPhi = pgun_params.getParameter<unsigned int>("NParticlesInPhi");
  fPhiRandomSpread = pgun_params.getParameter<double>("PhiRandomSpread");

  produces<HepMCProduct>("unsmeared");
  produces<GenEventInfoProduct>();
}

MultiParticleProducer::~MultiParticleProducer() {}

void MultiParticleProducer::produce(edm::Event& e, const edm::EventSetup& es) {
  if (fVerbosity > 0) {
    edm::LogVerbatim("MultiParticleProducer") << "MultiParticleProducer : Begin New Event Generation";
  }

  edm::Service<edm::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine* engine = &rng->getEngine(e.streamID());

  // event loop (well, another step in it...)

  // no need to clean up GenEvent memory - done in HepMCProduct

  // here re-create fEvt (memory)
  //
  fEvt = new HepMC::GenEvent();

  // now actualy, cook up the event from PDGTable and gun parameters
  //

  // 1st, primary vertex
  //
  HepMC::GenVertex* Vtx = new HepMC::GenVertex(HepMC::FourVector(0., 0., 0.));

  // loop over particles
  //
  std::vector<double> etas;
  int barcode = 1;
  for (unsigned int i_eta = 0; i_eta < fNParticlesInEta; i_eta++) {
    double eta_central = CLHEP::RandFlat::shoot(engine, fMinEta, fMaxEta);
    if (std::all_of(etas.begin(), etas.end(), [&](double other_eta) { return std::abs(eta_central - other_eta) >= fEtaMinSeparation; })) {
      etas.push_back(eta_central);
      const double start_phi = CLHEP::RandFlat::shoot(engine, 0., (2.*M_PI));
      for (unsigned int i_phi = 0; i_phi < fNParticlesInPhi; i_phi++) {
        double average_phi = start_phi + i_phi * (2.*M_PI) / fNParticlesInPhi;
        double phi = reco::reduceRange(CLHEP::RandFlat::shoot(engine, average_phi-fPhiRandomSpread/2., average_phi+fPhiRandomSpread/2.));
        double eta = CLHEP::RandFlat::shoot(engine, eta_central-fEtaRandomSpread/2., eta_central+fEtaRandomSpread/2.);

        double energy = CLHEP::RandFlat::shoot(engine, fMinE, fMaxE);
        int PartID = fPartIDs[CLHEP::RandFlat::shootInt(engine, fPartIDs.size())];
        const HepPDT::ParticleData* PData = fPDGTable->particle(HepPDT::ParticleID(abs(PartID)));
        double mass = PData->mass().value();
        double mom2 = energy * energy - mass * mass;
        double mom = (mom2 > 0. ? std::sqrt(mom2) : 0.);
        double theta = 2. * atan(exp(-eta));
        double px = mom * sin(theta) * cos(phi);
        double py = mom * sin(theta) * sin(phi);
        double pz = mom * cos(theta);

        HepMC::FourVector p(px, py, pz, energy);
        HepMC::GenParticle* Part = new HepMC::GenParticle(p, PartID, 1);
        Part->suggest_barcode(barcode);
        barcode++;
        Vtx->add_particle_out(Part);

        if (fAddAntiParticle) {
          HepMC::FourVector ap(-px, -py, -pz, energy);
          int APartID = -PartID;
          if (PartID == 22 || PartID == 23) {
            APartID = PartID;
          }
          HepMC::GenParticle* APart = new HepMC::GenParticle(ap, APartID, 1);
          APart->suggest_barcode(barcode);
          barcode++;
          Vtx->add_particle_out(APart);
        }
      }
    }
  }
  fEvt->add_vertex(Vtx);
  fEvt->set_event_number(e.id().event());
  fEvt->set_signal_process_id(20);

  if (fVerbosity > 0) {
    fEvt->print();
  }

  std::unique_ptr<HepMCProduct> BProduct(new HepMCProduct());
  BProduct->addHepMCData(fEvt);
  e.put(std::move(BProduct), "unsmeared");

  std::unique_ptr<GenEventInfoProduct> genEventInfo(new GenEventInfoProduct(fEvt));
  e.put(std::move(genEventInfo));

  if (fVerbosity > 0) {
    edm::LogVerbatim("MultiParticleProducer") << "MultiParticleProducer : Event Generation Done";
  }
}
