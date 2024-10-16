#ifndef MultiParticleProducer_H
#define MultiParticleProducer_H

#include "IOMC/ParticleGuns/interface/BaseFlatGunProducer.h"

namespace edm {
  // Particle gun producing multiple particles spread in eta and phi at the vertex, with random flat energies
  class MultiParticleProducer : public BaseFlatGunProducer {
  public:
    MultiParticleProducer(const ParameterSet &);
    ~MultiParticleProducer() override;

  private:
    void produce(Event &e, const EventSetup &es) override;

  protected:
    // data members

    double fMinE;
    double fMaxE;

    unsigned int fNParticlesInEta;
    double fEtaMinSeparation; 
    double fEtaRandomSpread;
    unsigned int fNParticlesInPhi;
    double fPhiRandomSpread;
  };
}  // namespace edm

#endif
