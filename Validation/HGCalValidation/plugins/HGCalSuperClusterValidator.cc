// Validation plots for superclusters in HGCal (TICL trackster dataformat)
// Distributions of multiplicity, ET, etc
// Energy resolution wrt electron/photon CaloParticle using association scores for matching
// Author : Théo Cuisset (theo.cuisset@polytechnique.edu)

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMServices/Core/interface/DQMGlobalEDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"
// #include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include <string>
#include <vector>

using dqm::reco::MonitorElement;
using edm::InputTag;
using std::vector;
using ticl::Trackster;
using ticl::TracksterCollection;
using TracksterToTracksterMap =
    ticl::AssociationMap<ticl::mapWithSharedEnergyAndScore, vector<ticl::Trackster>, vector<ticl::Trackster>>;

struct SuperClusterHistos;

class HGCalSuperClusterValidator : public DQMGlobalEDAnalyzer<SuperClusterHistos> {
public:
  explicit HGCalSuperClusterValidator(const edm::ParameterSet&);
  ~HGCalSuperClusterValidator() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&, SuperClusterHistos&) const override;

  void dqmAnalyze(edm::Event const&, edm::EventSetup const&, SuperClusterHistos const&) const override;

  const edm::EDGetTokenT<TracksterCollection> sc_tracksters_token;
  const edm::EDGetTokenT<TracksterCollection> tracksters_before_linking_token;
  const edm::EDGetTokenT<vector<vector<unsigned int>>> linkedTracksterIdToInputTracksterId_token;

  const edm::EDGetTokenT<TracksterToTracksterMap> simToRecoAssociator_token;

  const float simToRecoScoreThreshold_forEfficiency_, simToRecoScoreThreshold_forResolution_;
};

DEFINE_FWK_MODULE(HGCalSuperClusterValidator);

HGCalSuperClusterValidator::HGCalSuperClusterValidator(const edm::ParameterSet& ps)
    : sc_tracksters_token(consumes<TracksterCollection>(ps.getParameter<InputTag>("sc_tracksters"))),
      tracksters_before_linking_token(
          consumes<TracksterCollection>(ps.getParameter<InputTag>("sc_tracksters_before_linking"))),
      linkedTracksterIdToInputTracksterId_token(
          consumes<vector<vector<unsigned int>>>(ps.getParameter<InputTag>("linkedTracksterIdToInputTracksterId"))),
      simToRecoAssociator_token(
          consumes<TracksterToTracksterMap>(ps.getParameter<edm::InputTag>("associatorsimToReco"))),
      simToRecoScoreThreshold_forEfficiency_(ps.getParameter<double>("simToRecoScoreThreshold_forEfficiency")),
      simToRecoScoreThreshold_forResolution_(ps.getParameter<double>("simToRecoScoreThreshold_forResolution")) {}

void HGCalSuperClusterValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("sc_tracksters")->setComment("Supercluster collection (as Trackster dataformat)");
  desc.add<edm::InputTag>("sc_tracksters_before_linking")
      ->setComment(
          "Collection of tracksters that were used to build the superclusters (ie CLUE3D tracksters before "
          "superclustering)");
  desc.add<edm::InputTag>("linkedTracksterIdToInputTracksterId")
      ->setComment("Vector linking a supercluster back to its constituent tracksters");

  desc.add<edm::InputTag>("associatorsimToReco")
      ->setComment(
          "TICL association map mapping reco tracksters (sc_tracksters collection) to simTracksters (CaloParticle)");
  desc.add<double>("simToRecoScoreThreshold_forEfficiency", 0.3)
      ->setComment(
          "Threshold on the Sim->Reco score between SimTrackster from CaloParticle and Supercluster trackster "
          "collection to consider a supercluster correctly reconstructed (for efficiency plots)");
  desc.add<double>("simToRecoScoreThreshold_forResolution", 0.9)
      ->setComment(
          "Threshold on the Sim->Reco score between SimTrackster from CaloParticle and Supercluster trackster "
          "collection to consider a supercluster for energy scale and resolution plots (typically a looseer cut than "
          "for efficiency, otherwise the energy scale plots are one-sided)");

  descriptions.add("hgcalSuperClusterValidator", desc);
}

struct SuperClusterHistos {
  MonitorElement* sc_count;
  MonitorElement* sc_numConstituentTracksters;

  MonitorElement* sc_ET;
  MonitorElement* sc_eta;
  MonitorElement* sc_ET_vs_Eta;  // 2D

  MonitorElement* sim_ET_vs_Eta_num;    // 2D
  MonitorElement* sim_ET_vs_Eta_denom;  // 2D

  MonitorElement* sc_EoverEtruth;
  MonitorElement* sc_EoverEtruth_vs_ET;
  MonitorElement* sc_EoverEtruth_vs_ET_lowEta;
  MonitorElement* sc_EoverEtruth_vs_ET_midEta;
  MonitorElement* sc_EoverEtruth_vs_ET_highEta;
  MonitorElement* sc_EoverEtruth_vs_eta;
};

void HGCalSuperClusterValidator::bookHistograms(DQMStore::IBooker& ibook,
                                                edm::Run const& run,
                                                edm::EventSetup const& iSetup,
                                                SuperClusterHistos& histos) const {
  ibook.setCurrentFolder("HGCAL/SuperClusters/");
  // clang-format off
  const std::vector<float> sc_count_bins = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 28, 30, 35, 40, 50, 100, 500};
  const std::vector<float> ptBins = {0, 1, 2, 3, 4, 5, 7, 10, 15, 20, 30, 40, 50, 70, 100, 200, 300, 500};
  const std::vector<float> etaBins = {1.6, 1.9, 2.1, 2.3, 2.5, 2.6, 2.7, 2.8, 2.9, 3.};
  // clang-format on
  histos.sc_count = ibook.book1D("sc_count", "# SuperClusters in HGCAL", sc_count_bins.size() - 1, sc_count_bins.data()

  );
  histos.sc_numConstituentTracksters = ibook.book1D(
      "sc_numConstituentTracksters", "Number of constituent tracksters in SuperClusters in HGCAL", 20, 0, 20);

  histos.sc_ET = ibook.book1D("sc_ET", "SuperClusters ET in HGCAL", ptBins.size() - 1, ptBins.data());
  histos.sc_eta = ibook.book1D("sc_eta", "SuperClusters abs(eta) in HGCAL", etaBins.size() - 1, etaBins.data());
  histos.sc_ET_vs_Eta = ibook.book2D("sc_ET_vs_Eta",
                                     "SuperClusters in HGCAL pt-abs(eta)",
                                     ptBins.size() - 1,
                                     ptBins.data(),
                                     etaBins.size() - 1,
                                     etaBins.data());

  histos.sim_ET_vs_Eta_num =
      ibook.book2D("sim_ET_vs_Eta_num",
                   std::string("SuperClusters in HGCAL pt-abs(eta) (matched to CaloParticle, reco2sim<") +
                       simToRecoScoreThreshold_forEfficiency_ + ")",
                   ptBins.size() - 1,
                   ptBins.data(),
                   etaBins.size() - 1,
                   etaBins.data());
  histos.sim_ET_vs_Eta_num->setXTitle(
      "CaloParticle regressed energy (ie sum of sim energy of electron/photon in calorimeters) [GeV]");
  histos.sim_ET_vs_Eta_denom = ibook.book2D("sim_ET_vs_Eta_denom",
                                            "Electron/photon CaloParticle ",
                                            ptBins.size() - 1,
                                            ptBins.data(),
                                            etaBins.size() - 1,
                                            etaBins.data());
  histos.sim_ET_vs_Eta_denom->setXTitle(
      "CaloParticle regressed energy (ie sum of sim energy of electron/photon in calorimeters) [GeV]");

  histos.sc_EoverEtruth =
      ibook.book1D("sc_EoverEtruth", "SuperClusters in HGCAL E/Etruth (from CaloParticle)", 50, 0, 2);

  auto constexpr EoverEtruthBins_range = std::views::iota(0, 50)  // evenly spaced values between 0 and 2
                                         | std::views::transform([step = 2.0 / 49](int i) { return i * step; });
  std::vector<float> const EoverEtruthBins(EoverEtruthBins_range.begin(), EoverEtruthBins_range.end());
  histos.sc_EoverEtruth_vs_ET = ibook.book2D("sc_EoverEtruth_vs_ET",
                                             "SuperClusters in HGCAL E/Etruth (from CaloParticle)",
                                             ptBins.size() - 1,
                                             ptBins.data(),
                                             EoverEtruthBins.size() - 1,
                                             EoverEtruthBins.data());
  histos.sc_EoverEtruth_vs_ET->setXTitle(
      "Transverse energy from CaloParticle (ie sum of sim transverse energy of electron/photon in calorimeters) [GeV]");
  histos.sc_EoverEtruth_vs_ET->setYTitle(
      "Reconstructed raw energy of supercluster / energy from CaloParticle (Ereco/Etrue)");

  histos.sc_EoverEtruth_vs_ET_lowEta = ibook.book2D("sc_EoverEtruth_vs_ET_lowEta",
                                                    "SuperClusters in HGCAL E/Etruth (from CaloParticle) abs(eta)<2.1",
                                                    ptBins.size() - 1,
                                                    ptBins.data(),
                                                    EoverEtruthBins.size() - 1,
                                                    EoverEtruthBins.data());
  histos.sc_EoverEtruth_vs_ET_lowEta->setXTitle(
      "Transverse energy from CaloParticle (ie sum of sim transverse energy of electron/photon in calorimeters) [GeV]");
  histos.sc_EoverEtruth_vs_ET_lowEta->setYTitle(
      "Reconstructed raw energy of supercluster / energy from CaloParticle (Ereco/Etrue)");

  histos.sc_EoverEtruth_vs_ET_midEta =
      ibook.book2D("sc_EoverEtruth_vs_ET_midEta",
                   "SuperClusters in HGCAL E/Etruth (from CaloParticle) 2.1<abs(eta)<2.6",
                   ptBins.size() - 1,
                   ptBins.data(),
                   EoverEtruthBins.size() - 1,
                   EoverEtruthBins.data());
  histos.sc_EoverEtruth_vs_ET_midEta->setXTitle(
      "Transverse energy from CaloParticle (ie sum of sim transverse energy of electron/photon in calorimeters) [GeV]");
  histos.sc_EoverEtruth_vs_ET_midEta->setYTitle(
      "Reconstructed raw energy of supercluster / energy from CaloParticle (Ereco/Etrue)");

  histos.sc_EoverEtruth_vs_ET_highEta = ibook.book2D("sc_EoverEtruth_vs_ET_highEta",
                                                     "SuperClusters in HGCAL E/Etruth (from CaloParticle) abs(eta)>2.6",
                                                     ptBins.size() - 1,
                                                     ptBins.data(),
                                                     EoverEtruthBins.size() - 1,
                                                     EoverEtruthBins.data());
  histos.sc_EoverEtruth_vs_ET_highEta->setXTitle(
      "Transverse energy from CaloParticle (ie sum of sim transverse energy of electron/photon in calorimeters) [GeV]");
  histos.sc_EoverEtruth_vs_ET_highEta->setYTitle(
      "Reconstructed raw energy of supercluster / energy from CaloParticle (Ereco/Etrue)");
  histos.sc_EoverEtruth_vs_eta = ibook.book2D("sc_EoverEtruth_vs_eta",
                                              "SuperClusters in HGCAL E/Etruth (from CaloParticle)",
                                              etaBins.size() - 1,
                                              etaBins.data(),
                                              EoverEtruthBins.size() - 1,
                                              EoverEtruthBins.data());
  histos.sc_EoverEtruth_vs_eta->setYTitle(
      "Reconstructed raw energy of supercluster / energy from CaloParticle (Ereco/Etrue)");
}

void HGCalSuperClusterValidator::dqmAnalyze(edm::Event const& iEvent,
                                            edm::EventSetup const& iSetup,
                                            SuperClusterHistos const& histos) const {
  TracksterCollection const& sc_tracksters = iEvent.get(sc_tracksters_token);
  vector<vector<unsigned int>> const& linkedTracksterIdToInputTracksterId =
      iEvent.get(linkedTracksterIdToInputTracksterId_token);
  assert(sc_tracksters.size() == linkedTracksterIdToInputTracksterId.size());

  histos.sc_count->Fill(sc_tracksters.size());

  for (std::size_t i_ts = 0; i_ts < sc_tracksters.size(); i_ts++) {
    Trackster const& sc_ts = sc_tracksters[i_ts];
    histos.sc_eta->Fill(std::abs(sc_ts.barycenter().eta()));
    histos.sc_ET->Fill(sc_ts.raw_energy());
    histos.sc_ET_vs_Eta->Fill(sc_ts.raw_pt(), std::abs(sc_ts.barycenter().eta()));

    histos.sc_numConstituentTracksters->Fill(linkedTracksterIdToInputTracksterId[i_ts].size());
  }

  if (simToRecoAssociator_token.isUninitialized())
    return;

  TracksterToTracksterMap const& simToRecoMap = iEvent.get(simToRecoAssociator_token);
  for (std::size_t i_simts = 0; i_simts < simToRecoMap.size(); i_simts++) {
    Trackster const& simts = *simToRecoMap.getRefFirst(i_simts);
    if (simts.isHadronic())
      continue;  // Consider only gen photons/electrons

    const double simET = simts.regressed_energy() / std::cosh(simts.barycenter().eta());
    ;
    histos.sim_ET_vs_Eta_denom->Fill(simET, std::abs(simts.barycenter().eta()));

    // Find the supercluster with the highest shared energy with the SimTrackster
    auto best_supercluster_it = std::ranges::max_element(
        simToRecoMap[i_simts],
        [](TracksterToTracksterMap::AssociationElementType const& a,
           TracksterToTracksterMap::AssociationElementType const& b) { return a.sharedEnergy() < b.sharedEnergy(); });
    if (best_supercluster_it == simToRecoMap[i_simts].end())
      continue;  // SimTrackster without any associated supercluster

    Trackster const& bestSupercluster = *simToRecoMap.getRefSecond(best_supercluster_it->index());

    if (best_supercluster_it->score() < simToRecoScoreThreshold_forEfficiency_) {
      // Ignore in cases where the supercluster is not well reconstructed ()
      histos.sim_ET_vs_Eta_num->Fill(simET, std::abs(simts.barycenter().eta()));
    }

    if (best_supercluster_it->score() < simToRecoScoreThreshold_forResolution_) {
      // Ignore in cases where the supercluster is too badly reconstructed (will skip superclusters with very low shared energy, which are PU but just have a little sim energy from eg. wide angle scattering)
      // The SimTs regressed_energy is the sum of simhits energies in the calorimeter
      const double EoverTrue = bestSupercluster.raw_energy() / simts.regressed_energy();

      histos.sc_EoverEtruth->Fill(EoverTrue);
      histos.sc_EoverEtruth_vs_eta->Fill(bestSupercluster.barycenter().eta(), EoverTrue);
      histos.sc_EoverEtruth_vs_ET->Fill(simET, EoverTrue);
      if (std::abs(bestSupercluster.barycenter().eta()) < 2.1)
        histos.sc_EoverEtruth_vs_ET_lowEta->Fill(simET, EoverTrue);
      else if (std::abs(bestSupercluster.barycenter().eta()) < 2.6)
        histos.sc_EoverEtruth_vs_ET_midEta->Fill(simET, EoverTrue);
      else
        histos.sc_EoverEtruth_vs_ET_highEta->Fill(simET, EoverTrue);
    }
  }
}
