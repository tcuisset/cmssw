
#include <memory>

#include "DataFormats/Math/interface/deltaPhi.h"

#include "GeneratorInterface/Core/interface/GeneratorFilter.h"
#include "GeneratorInterface/ExternalDecays/interface/ExternalDecayDriver.h"

#include "GeneratorInterface/Pythia8Interface/interface/Py8GunBase.h"

#include "CLHEP/Random/RandFlat.h"

namespace gen {

  /* Pythia-gun equivalent of MultiParticleProducer. Generates multiple particles with customisable separation between them. 
  Not mean for quarks and gluons. Suitable for taus, pi-zeros. */
  class Py8EMultiParticleGun : public Py8GunBase {
  public:
    Py8EMultiParticleGun(edm::ParameterSet const&);
    ~Py8EMultiParticleGun() override {}

    bool generatePartonsAndHadronize() override;
    const char* classname() const override;

  private:
    // EGun particle(s) characteristics
    double fMinEta;
    double fMaxEta;
    double fMinE;
    double fMaxE;
    unsigned int fNParticlesInEta;
    double fEtaMinSeparation; 
    double fEtaRandomSpread;
    unsigned int fNParticlesInPhi;
    double fPhiRandomSpread;
  };

  // implementation
  //
  Py8EMultiParticleGun::Py8EMultiParticleGun(edm::ParameterSet const& ps) : Py8GunBase(ps) {
    // ParameterSet defpset ;
    edm::ParameterSet pgun_params = ps.getParameter<edm::ParameterSet>("PGunParameters");  // , defpset ) ;
    fMinEta = pgun_params.getParameter<double>("MinEta");                                  // ,-2.2);
    fMaxEta = pgun_params.getParameter<double>("MaxEta");                                  // , 2.2);

    fMinE = pgun_params.getParameter<double>("MinE");
    fMaxE = pgun_params.getParameter<double>("MaxE");
  
    fNParticlesInEta = pgun_params.getParameter<unsigned int>("NParticlesInEta");
    fEtaMinSeparation = pgun_params.getParameter<double>("EtaMinSeparation");
    fEtaRandomSpread = pgun_params.getParameter<double>("EtaRandomSpread");
    fNParticlesInPhi = pgun_params.getParameter<unsigned int>("NParticlesInPhi");
    fPhiRandomSpread = pgun_params.getParameter<double>("PhiRandomSpread");

    if (std::find_if(fPartIDs.begin(), fPartIDs.end(), [](int particleID) { return std::abs(particleID) <= 6 || particleID == 21;}) != fPartIDs.end()) {
      throw cms::Exception("PythiaError") << "Attempting to generate quarks or gluons with Py8EMultiParticleGun. This will not handle color properly."
                                          << std::endl;
    }
  }

  bool Py8EMultiParticleGun::generatePartonsAndHadronize() {
    fMasterGen->event.reset();

    edm::Service<edm::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine* engine = &rng->getEngine(getEDMEvent().streamID());

    int NTotParticles = fPartIDs.size();

    // energy below is dummy, it is not used
    (fMasterGen->event).append(990, -11, 0, 0, 2, 1 + NTotParticles, 0, 0, 0., 0., 0., 15000., 15000.);

    std::vector<double> etas;
    for (unsigned int i_eta = 0; i_eta < fNParticlesInEta; i_eta++) {
      double eta_central = CLHEP::RandFlat::shoot(engine, fMinEta, fMaxEta);
      if (std::all_of(etas.begin(), etas.end(), [&](double other_eta) { return std::abs(eta_central - other_eta) >= fEtaMinSeparation; })) {
        etas.push_back(eta_central);
        const double start_phi = CLHEP::RandFlat::shoot(engine, 0., (2.*M_PI));
        for (unsigned int i_phi = 0; i_phi < fNParticlesInPhi; i_phi++) {
          double average_phi = start_phi + i_phi * (2.*M_PI) / fNParticlesInPhi;
          double phi = reco::reducePhiRange(CLHEP::RandFlat::shoot(engine, average_phi-fPhiRandomSpread/2., average_phi+fPhiRandomSpread/2.));
          double eta = CLHEP::RandFlat::shoot(engine, eta_central-fEtaRandomSpread/2., eta_central+fEtaRandomSpread/2.);
  
          double ee = CLHEP::RandFlat::shoot(engine, fMinE, fMaxE);
          int particleID = fPartIDs[CLHEP::RandFlat::shootInt(engine, fPartIDs.size())];  // this is PDG - need to convert to Py8 ???

          double the = 2. * atan(exp(-eta));
          double mass = (fMasterGen->particleData).m0(particleID);

          double pp = sqrt(ee * ee - mass * mass);
          double px = pp * sin(the) * cos(phi);
          double py = pp * sin(the) * sin(phi);
          double pz = pp * cos(the);

          if (!((fMasterGen->particleData).isParticle(particleID))) {
            particleID = std::fabs(particleID);
          }
          (fMasterGen->event).append(particleID, 1, 1, 0, 0, 0, 0, 0, px, py, pz, ee, mass);
          int eventSize = (fMasterGen->event).size() - 1;
          // -log(flat) = exponential distribution
          double tauTmp = -(fMasterGen->event)[eventSize].tau0() * log(randomEngine().flat());
          (fMasterGen->event)[eventSize].tau(tauTmp);
        }
      }
    }

    if (!fMasterGen->next())
      return false;
    evtGenDecay();

    event() = std::make_unique<HepMC::GenEvent>();
    return toHepMC.fill_next_event(fMasterGen->event, event().get());
  }

  const char* Py8EMultiParticleGun::classname() const { return "Py8EMultiParticleGun"; }

  typedef edm::GeneratorFilter<gen::Py8EMultiParticleGun, gen::ExternalDecayDriver> Pythia8EMultiParticleGun;

}  // namespace gen

using gen::Pythia8EMultiParticleGun;
DEFINE_FWK_MODULE(Pythia8EMultiParticleGun);
