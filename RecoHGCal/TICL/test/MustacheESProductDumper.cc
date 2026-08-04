#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "CondFormats/DataRecord/interface/EcalMustacheSCParametersRcd.h"
#include "CondFormats/DataRecord/interface/EcalSCDynamicDPhiParametersRcd.h"
#include "CondFormats/EcalObjects/interface/EcalMustacheSCParameters.h"
#include "CondFormats/EcalObjects/interface/EcalSCDynamicDPhiParameters.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "RecoEcal/EgammaCoreTools/interface/Mustache.h"

namespace {
  class MustacheESProductDumper : public edm::one::EDAnalyzer<> {
  public:
    explicit MustacheESProductDumper(edm::ParameterSet const&);
    static void fillDescriptions(edm::ConfigurationDescriptions&);

  private:
    void analyze(edm::Event const&, edm::EventSetup const&) override;

    edm::ESGetToken<EcalMustacheSCParameters, EcalMustacheSCParametersRcd> mustacheToken_;
    edm::ESGetToken<EcalSCDynamicDPhiParameters, EcalSCDynamicDPhiParametersRcd> dynamicDPhiToken_;
    double seedEta_;
    double seedPhi_;
    double clusterEt_;
    double clusterEnergy_;
    double etaHalfRange_;
    double phiHalfRange_;
    unsigned int nEta_;
    unsigned int nPhi_;
    std::string output_;
  };

  MustacheESProductDumper::MustacheESProductDumper(edm::ParameterSet const& config)
      : mustacheToken_(esConsumes()),
        dynamicDPhiToken_(esConsumes()),
        seedEta_(config.getParameter<double>("seedEta")),
        seedPhi_(config.getParameter<double>("seedPhi")),
        clusterEt_(config.getParameter<double>("clusterEt")),
        clusterEnergy_(config.getParameter<double>("clusterEnergy")),
        etaHalfRange_(config.getParameter<double>("etaHalfRange")),
        phiHalfRange_(config.getParameter<double>("phiHalfRange")),
        nEta_(config.getParameter<unsigned int>("nEta")),
        nPhi_(config.getParameter<unsigned int>("nPhi")),
        output_(config.getParameter<std::string>("output")) {
    if ((clusterEt_ > 0.) == (clusterEnergy_ > 0.)) {
      throw cms::Exception("Configuration") << "Set exactly one of clusterEt and clusterEnergy to a positive value.";
    }
    if (nEta_ < 2 || nPhi_ < 2 || etaHalfRange_ <= 0. || phiHalfRange_ <= 0.) {
      throw cms::Exception("Configuration")
          << "nEta and nPhi must be at least 2, and both half-ranges must be positive.";
    }
  }

  void MustacheESProductDumper::analyze(edm::Event const&, edm::EventSetup const& setup) {
    auto const& mustache = setup.getData(mustacheToken_);
    auto const& dynamicDPhi = setup.getData(dynamicDPhiToken_);

    std::cout << "\n========== EcalMustacheSCParameters ==========\n"
              << mustache << "\n========== EcalSCDynamicDPhiParameters ==========\n"
              << dynamicDPhi << '\n';

    const double centralEnergy = clusterEnergy_ > 0. ? clusterEnergy_ : clusterEt_ * std::cosh(seedEta_);
    auto const* parabola = mustache.parabolaParameters(std::log10(centralEnergy), std::abs(seedEta_));
    auto const* dphi = dynamicDPhi.dynamicDPhiParameters(centralEnergy, std::abs(seedEta_));
    std::cout << "Central lookup at candidate eta = seed eta = " << seedEta_ << ", candidate E = " << centralEnergy
              << " GeV";
    if (parabola) {
      std::cout << "\n  Mustache bin: log10EMin=" << parabola->log10EMin << ", etaMin=" << parabola->etaMin;
    } else {
      std::cout << "\n  Mustache bin: none";
    }
    if (dphi) {
      std::cout << "\n  dynamic-dPhi bin: eMin=" << dphi->eMin << ", etaMin=" << dphi->etaMin;
    } else {
      std::cout << "\n  dynamic-dPhi bin: none";
    }
    std::cout << "\n\n";

    std::ofstream out(output_);
    if (!out) {
      throw cms::Exception("FileOpenError") << "Cannot open output file '" << output_ << "'.";
    }
    out << std::setprecision(12);
    out << "# seed_eta=" << seedEta_ << "\n"
        << "# seed_phi=" << seedPhi_ << "\n"
        << "# cluster_et=" << clusterEt_ << "\n"
        << "# cluster_energy=" << clusterEnergy_ << "\n"
        << "# eta_half_range=" << etaHalfRange_ << "\n"
        << "# phi_half_range=" << phiHalfRange_ << "\n"
        << "eta,phi,candidate_energy,in_mustache,in_dynamic_dphi,in_combined\n";

    for (unsigned int iEta = 0; iEta < nEta_; ++iEta) {
      const double eta = seedEta_ - etaHalfRange_ + 2. * etaHalfRange_ * static_cast<double>(iEta) / (nEta_ - 1);
      const double energy = clusterEnergy_ > 0. ? clusterEnergy_ : clusterEt_ * std::cosh(eta);
      for (unsigned int iPhi = 0; iPhi < nPhi_; ++iPhi) {
        const double phi = seedPhi_ - phiHalfRange_ + 2. * phiHalfRange_ * static_cast<double>(iPhi) / (nPhi_ - 1);
        const bool inMustache = reco::MustacheKernel::inMustache(&mustache, seedEta_, seedPhi_, energy, eta, phi);
        const bool inDynamicDPhi =
            reco::MustacheKernel::inDynamicDPhiWindow(&dynamicDPhi, seedEta_, seedPhi_, energy, eta, phi);
        out << eta << ',' << phi << ',' << energy << ',' << inMustache << ',' << inDynamicDPhi << ','
            << (inMustache && inDynamicDPhi) << '\n';
      }
    }
    std::cout << "Wrote exact MustacheKernel scan to " << output_ << '\n';
  }

  void MustacheESProductDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<double>("seedEta", 2.0);
    desc.add<double>("seedPhi", 0.0);
    desc.add<double>("clusterEt", 10.0)
        ->setComment("Candidate ET in GeV; E=ET*cosh(candidate eta). Set <=0 when using clusterEnergy.");
    desc.add<double>("clusterEnergy", -1.0)
        ->setComment("Fixed candidate raw energy in GeV. Set <=0 when using clusterEt.");
    desc.add<double>("etaHalfRange", 0.3);
    desc.add<double>("phiHalfRange", 0.6);
    desc.add<unsigned int>("nEta", 301);
    desc.add<unsigned int>("nPhi", 401);
    desc.add<std::string>("output", "mustache_scan.csv");
    descriptions.add("mustacheESProductDumper", desc);
  }
}  // namespace

DEFINE_FWK_MODULE(MustacheESProductDumper);
