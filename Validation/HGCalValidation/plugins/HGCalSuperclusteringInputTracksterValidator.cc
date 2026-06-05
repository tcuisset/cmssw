// Validation plots for trackster behaviour as input to superclusters
// Author : Theo Cuisset (theo.cuisset@polytechnique.edu)

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <numeric>
#include <string>
#include <vector>

#include "DQMServices/Core/interface/DQMGlobalEDAnalyzer.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "TH3F.h"

using namespace ticl;

namespace {
  using dqm::reco::MonitorElement;

  double pidScore(Trackster const& trackster, std::vector<Trackster::ParticleType> const& pidsToConsider) {
    return std::transform_reduce(
        pidsToConsider.begin(), pidsToConsider.end(), 0., std::plus<>{}, [&trackster](Trackster::ParticleType partType) {
          return trackster.id_probability(partType);
        });
  }

  std::size_t etaBinIndex(std::vector<double> const& absEtaBins, double absEta) {
    auto upper = std::upper_bound(absEtaBins.begin(), absEtaBins.end(), absEta);
    if (upper == absEtaBins.begin() || upper == absEtaBins.end()) {
      return absEtaBins.size();
    }
    return static_cast<std::size_t>(std::distance(absEtaBins.begin(), upper) - 1);
  }

  void fillPair(MonitorElement* histogram, double deltaEta, double deltaPhi) {
    if (histogram != nullptr) {
      histogram->Fill(deltaEta, deltaPhi);
    }
  }

  void fillPair(std::vector<MonitorElement*> const& histograms, std::size_t etaBin, double deltaEta, double deltaPhi) {
    if (etaBin < histograms.size()) {
      fillPair(histograms[etaBin], deltaEta, deltaPhi);
    }
  }

  std::string etaBinSuffix(std::size_t index) { return "_absEtaBin" + std::to_string(index); }
}  // namespace

struct HistogramsSuperclusteringInputTracksters {
  MonitorElement* deltaEta_deltaPhi_toSeed_;
  MonitorElement* deltaEta_deltaPhi_toSeed_realBrem_;
  MonitorElement* deltaEta_deltaPhi_toSeed_fake_;

  MonitorElement* deltaEta_deltaPhi_toSeed_afterCandidatePID_;
  MonitorElement* deltaEta_deltaPhi_toSeed_realBrem_afterCandidatePID_;
  MonitorElement* deltaEta_deltaPhi_toSeed_fake_afterCandidatePID_;

  MonitorElement* deltaEta_deltaPhi_candidatePID_toSeed_;
  MonitorElement* deltaEta_deltaPhi_candidatePID_toSeed_realBrem_;
  MonitorElement* deltaEta_deltaPhi_candidatePID_toSeed_fake_;

  std::vector<MonitorElement*> deltaEta_deltaPhi_toSeed_byAbsEta_;
  std::vector<MonitorElement*> deltaEta_deltaPhi_toSeed_realBrem_byAbsEta_;
  std::vector<MonitorElement*> deltaEta_deltaPhi_toSeed_fake_byAbsEta_;
};

class HGCalSuperclusteringInputTracksterValidator
    : public DQMGlobalEDAnalyzer<HistogramsSuperclusteringInputTracksters> {
public:
  explicit HGCalSuperclusteringInputTracksterValidator(const edm::ParameterSet&);
  ~HGCalSuperclusteringInputTracksterValidator() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void bookHistograms(DQMStore::IBooker&,
                      edm::Run const&,
                      edm::EventSetup const&,
                      HistogramsSuperclusteringInputTracksters&) const override;

  void dqmAnalyze(edm::Event const&,
                  edm::EventSetup const&,
                  HistogramsSuperclusteringInputTracksters const&) const override;

  const std::string folder_;
  const edm::EDGetTokenT<ticl::TracksterCollection> tracksters_token_;
  const edm::EDGetTokenT<std::vector<int>> seedTracksterMask_token_;
  const edm::EDGetTokenT<std::vector<int>> candidateTracksterMask_token_;
  const edm::EDGetTokenT<std::vector<int>> fakeTracksterMask_token_;

  const double pidCut_;
  const std::vector<ticl::Trackster::ParticleType> pidsToConsider_;

  const double deltaEtaWindow_;
  const double deltaPhiWindow_;
  const unsigned int pidBins_;
  const unsigned int deltaEtaBins_;
  const unsigned int deltaPhiBins_;
  const std::vector<double> absEtaBins_;
};

HGCalSuperclusteringInputTracksterValidator::HGCalSuperclusteringInputTracksterValidator(
    const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")),
      tracksters_token_(consumes<ticl::TracksterCollection>(iConfig.getParameter<edm::InputTag>("tracksters"))),
      seedTracksterMask_token_(consumes<std::vector<int>>(iConfig.getParameter<edm::InputTag>("seedTracksterMask"))),
      candidateTracksterMask_token_(
          consumes<std::vector<int>>(iConfig.getParameter<edm::InputTag>("candidateTracksterMask"))),
      fakeTracksterMask_token_(consumes<std::vector<int>>(iConfig.getParameter<edm::InputTag>("fakeTracksterMask"))),
      pidCut_(iConfig.getParameter<double>("pidCut")),
      pidsToConsider_({ticl::Trackster::ParticleType::electron, ticl::Trackster::ParticleType::photon}),
      deltaEtaWindow_(iConfig.getParameter<double>("deltaEtaWindow")),
      deltaPhiWindow_(iConfig.getParameter<double>("deltaPhiWindow")),
      pidBins_(iConfig.getParameter<unsigned int>("pidBins")),
      deltaEtaBins_(iConfig.getParameter<unsigned int>("deltaEtaBins")),
      deltaPhiBins_(iConfig.getParameter<unsigned int>("deltaPhiBins")),
      absEtaBins_(iConfig.getParameter<std::vector<double>>("absEtaBins")) {
  if (deltaEtaWindow_ <= 0. || deltaPhiWindow_ <= 0.) {
    throw cms::Exception("Configuration") << "deltaEtaWindow and deltaPhiWindow must be positive.";
  }
  if (deltaEtaBins_ == 0 || deltaPhiBins_ == 0) {
    throw cms::Exception("Configuration") << "deltaEtaBins and deltaPhiBins must be non-zero.";
  }
  if (absEtaBins_.size() < 2 || !std::ranges::is_sorted(absEtaBins_)) {
    throw cms::Exception("Configuration") << "absEtaBins must contain at least two sorted bin edges.";
  }
  if (pidBins_ == 0) {
    throw cms::Exception("Configuration") << "pidBins must be non-zero.";
  }
}

void HGCalSuperclusteringInputTracksterValidator::dqmAnalyze(
    edm::Event const& iEvent,
    edm::EventSetup const& iSetup,
    HistogramsSuperclusteringInputTracksters const& histos) const {
  ticl::TracksterCollection const& tracksters = iEvent.get(tracksters_token_);
  std::vector<int> const& seedTracksterMask = iEvent.get(seedTracksterMask_token_);
  std::vector<int> const& candidateTracksterMask = iEvent.get(candidateTracksterMask_token_);
  std::vector<int> const& fakeTracksterMask = iEvent.get(fakeTracksterMask_token_);
  assert(seedTracksterMask.size() == tracksters.size());
  assert(candidateTracksterMask.size() == tracksters.size());
  assert(fakeTracksterMask.size() == tracksters.size());

  std::array<TICLLayerTile, 2> tracksterTilesBothEndcaps;
  for (unsigned int i = 0; i < tracksters.size(); ++i) {
    auto const& ts = tracksters[i];
    tracksterTilesBothEndcaps[ts.barycenter().eta() > 0.f].fill(ts.barycenter().eta(), ts.barycenter().phi(), i);
  }

  for (std::size_t seed_i = 0; seed_i < tracksters.size(); ++seed_i) {
    if (seedTracksterMask[seed_i] != 0) {
      continue;
    }

    ticl::Trackster const& ts_seed = tracksters[seed_i];
    auto& tiles = tracksterTilesBothEndcaps[ts_seed.barycenter().eta() > 0.f];
    const auto search_box = tiles.searchBoxEtaPhi(ts_seed.barycenter().eta() - deltaEtaWindow_,
                                                  ts_seed.barycenter().eta() + deltaEtaWindow_,
                                                  ts_seed.barycenter().phi() - deltaPhiWindow_,
                                                  ts_seed.barycenter().phi() + deltaPhiWindow_);

    for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
      for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
        const int phiBin = (phi_i % TileConstants::nPhiBins + TileConstants::nPhiBins) % TileConstants::nPhiBins;
        const auto bin = tiles.globalBin(eta_i, phiBin);
        for (unsigned int cand_i : tiles[bin]) {
          auto const& ts_cand = tracksters[cand_i];
          if (ts_cand.raw_pt() >= ts_seed.raw_pt()) {
            continue;
          }

          const double deltaEta = ts_cand.barycenter().eta() - ts_seed.barycenter().eta();
          const double deltaPhiValue = deltaPhi(ts_cand.barycenter().phi(), ts_seed.barycenter().phi());
          if (std::abs(deltaEta) >= deltaEtaWindow_ || std::abs(deltaPhiValue) >= deltaPhiWindow_) {
            continue;
          }

          const std::size_t absEtaBin = etaBinIndex(absEtaBins_, std::abs(ts_seed.barycenter().eta()));
          const double candidatePid = pidScore(ts_cand, pidsToConsider_);
          const bool candidatePassesPid = candidatePid > pidCut_;
          const bool isRealBrem = candidateTracksterMask[cand_i] == 0;
          const bool isFake = fakeTracksterMask[cand_i] == 0;

          fillPair(histos.deltaEta_deltaPhi_toSeed_, deltaEta, deltaPhiValue);
          histos.deltaEta_deltaPhi_candidatePID_toSeed_->Fill(deltaEta, deltaPhiValue, candidatePid);
          fillPair(histos.deltaEta_deltaPhi_toSeed_byAbsEta_, absEtaBin, deltaEta, deltaPhiValue);

          if (candidatePassesPid) {
            fillPair(histos.deltaEta_deltaPhi_toSeed_afterCandidatePID_, deltaEta, deltaPhiValue);
          }

          if (isRealBrem) {
            fillPair(histos.deltaEta_deltaPhi_toSeed_realBrem_, deltaEta, deltaPhiValue);
            histos.deltaEta_deltaPhi_candidatePID_toSeed_realBrem_->Fill(deltaEta, deltaPhiValue, candidatePid);
            fillPair(histos.deltaEta_deltaPhi_toSeed_realBrem_byAbsEta_, absEtaBin, deltaEta, deltaPhiValue);
            if (candidatePassesPid) {
              fillPair(histos.deltaEta_deltaPhi_toSeed_realBrem_afterCandidatePID_, deltaEta, deltaPhiValue);
            }
          }

          if (isFake) {
            fillPair(histos.deltaEta_deltaPhi_toSeed_fake_, deltaEta, deltaPhiValue);
            histos.deltaEta_deltaPhi_candidatePID_toSeed_fake_->Fill(deltaEta, deltaPhiValue, candidatePid);
            fillPair(histos.deltaEta_deltaPhi_toSeed_fake_byAbsEta_, absEtaBin, deltaEta, deltaPhiValue);
            if (candidatePassesPid) {
              fillPair(histos.deltaEta_deltaPhi_toSeed_fake_afterCandidatePID_, deltaEta, deltaPhiValue);
            }
          }
        }
      }
    }
  }
}

void HGCalSuperclusteringInputTracksterValidator::bookHistograms(
    DQMStore::IBooker& ibook,
    edm::Run const& run,
    edm::EventSetup const& iSetup,
    HistogramsSuperclusteringInputTracksters& histos) const {
  ibook.setCurrentFolder(folder_);

  auto bookDeltaEtaDeltaPhi = [&](std::string const& name, std::string const& title) {
    return ibook.book2D(
        name, title, deltaEtaBins_, -deltaEtaWindow_, deltaEtaWindow_, deltaPhiBins_, -deltaPhiWindow_, deltaPhiWindow_);
  };

  histos.deltaEta_deltaPhi_toSeed_ = bookDeltaEtaDeltaPhi(
      "deltaEta_deltaPhi_toSeed",
      "Candidate trackster position relative to seed;#Delta#eta(candidate, seed);#Delta#phi(candidate, seed)");
  histos.deltaEta_deltaPhi_toSeed_realBrem_ =
      bookDeltaEtaDeltaPhi("deltaEta_deltaPhi_toSeed_realBrem",
                           "Real-brem candidate trackster position relative to seed;#Delta#eta(candidate, "
                           "seed);#Delta#phi(candidate, seed)");
  histos.deltaEta_deltaPhi_toSeed_fake_ = bookDeltaEtaDeltaPhi(
      "deltaEta_deltaPhi_toSeed_fake",
      "Fake candidate trackster position relative to seed;#Delta#eta(candidate, seed);#Delta#phi(candidate, seed)");

  histos.deltaEta_deltaPhi_toSeed_afterCandidatePID_ =
      bookDeltaEtaDeltaPhi("deltaEta_deltaPhi_toSeed_afterCandidatePID",
                           "Candidate trackster position relative to seed after candidate PID "
                           "cut;#Delta#eta(candidate, seed);#Delta#phi(candidate, seed)");
  histos.deltaEta_deltaPhi_toSeed_realBrem_afterCandidatePID_ =
      bookDeltaEtaDeltaPhi("deltaEta_deltaPhi_toSeed_realBrem_afterCandidatePID",
                           "Real-brem candidate trackster position relative to seed after candidate PID "
                           "cut;#Delta#eta(candidate, seed);#Delta#phi(candidate, seed)");
  histos.deltaEta_deltaPhi_toSeed_fake_afterCandidatePID_ =
      bookDeltaEtaDeltaPhi("deltaEta_deltaPhi_toSeed_fake_afterCandidatePID",
                           "Fake candidate trackster position relative to seed after candidate PID "
                           "cut;#Delta#eta(candidate, seed);#Delta#phi(candidate, seed)");

  auto makeDeltaEtaDeltaPhiPid = [&](char const* name, char const* title) -> TH3F* {
    return new TH3F(name,
                    title,
                    deltaEtaBins_,
                    -deltaEtaWindow_,
                    deltaEtaWindow_,
                    deltaPhiBins_,
                    -deltaPhiWindow_,
                    deltaPhiWindow_,
                    pidBins_,
                    0.,
                    1.);
  };
  histos.deltaEta_deltaPhi_candidatePID_toSeed_ = ibook.book3D(
      "deltaEta_deltaPhi_candidatePID_toSeed",
      makeDeltaEtaDeltaPhiPid("deltaEta_deltaPhi_candidatePID_toSeed",
                              "Candidate trackster position and PID relative to seed;#Delta#eta(candidate, "
                              "seed);#Delta#phi(candidate, seed);candidate e/#gamma PID score"));
  histos.deltaEta_deltaPhi_candidatePID_toSeed_realBrem_ =
      ibook.book3D("deltaEta_deltaPhi_candidatePID_toSeed_realBrem",
                   makeDeltaEtaDeltaPhiPid("deltaEta_deltaPhi_candidatePID_toSeed_realBrem",
                                           "Real-brem candidate trackster position and PID relative to "
                                           "seed;#Delta#eta(candidate, seed);#Delta#phi(candidate, seed);candidate "
                                           "e/#gamma PID score"));
  histos.deltaEta_deltaPhi_candidatePID_toSeed_fake_ = ibook.book3D(
      "deltaEta_deltaPhi_candidatePID_toSeed_fake",
      makeDeltaEtaDeltaPhiPid("deltaEta_deltaPhi_candidatePID_toSeed_fake",
                              "Fake candidate trackster position and PID relative to seed;#Delta#eta(candidate, "
                              "seed);#Delta#phi(candidate, seed);candidate e/#gamma PID score"));

  const auto nEtaBins = absEtaBins_.size() - 1;
  histos.deltaEta_deltaPhi_toSeed_byAbsEta_.reserve(nEtaBins);
  histos.deltaEta_deltaPhi_toSeed_realBrem_byAbsEta_.reserve(nEtaBins);
  histos.deltaEta_deltaPhi_toSeed_fake_byAbsEta_.reserve(nEtaBins);
  for (std::size_t i = 0; i < nEtaBins; ++i) {
    const std::string etaRange =
        std::to_string(absEtaBins_[i]) + " < |#eta(seed)| < " + std::to_string(absEtaBins_[i + 1]);
    histos.deltaEta_deltaPhi_toSeed_byAbsEta_.push_back(
        bookDeltaEtaDeltaPhi("deltaEta_deltaPhi_toSeed" + etaBinSuffix(i),
                             "Candidate trackster position relative to seed, " + etaRange +
                                 ";#Delta#eta(candidate, seed);#Delta#phi(candidate, seed)"));
    histos.deltaEta_deltaPhi_toSeed_realBrem_byAbsEta_.push_back(
        bookDeltaEtaDeltaPhi("deltaEta_deltaPhi_toSeed_realBrem" + etaBinSuffix(i),
                             "Real-brem candidate trackster position relative to seed, " + etaRange +
                                 ";#Delta#eta(candidate, seed);#Delta#phi(candidate, seed)"));
    histos.deltaEta_deltaPhi_toSeed_fake_byAbsEta_.push_back(
        bookDeltaEtaDeltaPhi("deltaEta_deltaPhi_toSeed_fake" + etaBinSuffix(i),
                             "Fake candidate trackster position relative to seed, " + etaRange +
                                 ";#Delta#eta(candidate, seed);#Delta#phi(candidate, seed)"));
  }
}

void HGCalSuperclusteringInputTracksterValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("folder", "HGCAL/HGCalSuperclusteringInputTracksterValidator/")
      ->setComment("DQM folder. Please keep the trailing '/'.");
  desc.add<edm::InputTag>("tracksters", edm::InputTag("ticlTrackstersCLUE3DHigh"))
      ->setComment("Input trackster collection used as input to superclustering.");
  desc.add<edm::InputTag>("seedTracksterMask", edm::InputTag("ticlValidSuperclusteringSeedMask"))
      ->setComment("Mask selecting valid seed tracksters. A value of 0 means selected.");
  desc.add<edm::InputTag>("candidateTracksterMask", edm::InputTag("tracksterSuperclusteringValidCandidateMaskProducer"))
      ->setComment("Mask selecting real-brem candidate tracksters. A value of 0 means selected.");
  desc.add<edm::InputTag>("fakeTracksterMask", edm::InputTag("ticlValidSuperclusteringSeedMask", "fakes"))
      ->setComment("Mask selecting fake candidate tracksters. A value of 0 means selected.");

  desc.add<double>("pidCut", 0.2)->setComment("Cut on the candidate trackster electron+photon PID score.");
  desc.add<double>("deltaEtaWindow", 0.2)->setComment("Size of delta eta window used to select seed-candidate pairs.");
  desc.add<double>("deltaPhiWindow", 0.7)->setComment("Size of delta phi window used to select seed-candidate pairs.");
  desc.add<unsigned int>("deltaEtaBins", 80)->setComment("Number of histogram bins in delta eta.");
  desc.add<unsigned int>("deltaPhiBins", 80)->setComment("Number of histogram bins in delta phi.");
  desc.add<unsigned int>("pidBins", 50)
      ->setComment("Number of bins for the candidate trackster electron+photon PID score.");
  desc.add<std::vector<double>>("absEtaBins", {1.6, 1.8, 2.1, 2.3, 2.5, 2.7, 2.8, 2.9, 3.0})
      ->setComment("Abs(eta) bin edges for seed-eta-sliced deltaEta-deltaPhi histograms.");

  descriptions.add("hgcalSuperclusteringInputTracksterValidator", desc);
}

DEFINE_FWK_MODULE(HGCalSuperclusteringInputTracksterValidator);
