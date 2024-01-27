// Original Authors:  Philipp Zehetner, Wahid Redjeb

#include "TTree.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <variant>

#include <memory>  // unique_ptr
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/Common.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "SimDataFormats/Associations/interface/TracksterToSimTracksterHitLCAssociator.h"
#include "RecoHGCal/TICL/interface/commons.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Dumps a trackster collection 
class TracksterDumperHelper {
public:
  /**
   * TO be called once after tree creation. ev_event_ should be a pointer to the event counter
  */
  void initTree(TTree* trackster_tree_, unsigned int* ev_event_) {
    trackster_tree_->Branch("event", ev_event_);
    trackster_tree_->Branch("NTracksters", &nTracksters);
    trackster_tree_->Branch("NClusters", &nClusters);
    trackster_tree_->Branch("time", &trackster_time);
    trackster_tree_->Branch("timeError", &trackster_timeError);
    trackster_tree_->Branch("regressed_energy", &trackster_regressed_energy);
    trackster_tree_->Branch("raw_energy", &trackster_raw_energy);
    trackster_tree_->Branch("raw_em_energy", &trackster_raw_em_energy);
    trackster_tree_->Branch("raw_pt", &trackster_raw_pt);
    trackster_tree_->Branch("raw_em_pt", &trackster_raw_em_pt);
    trackster_tree_->Branch("barycenter_x", &trackster_barycenter_x);
    trackster_tree_->Branch("barycenter_y", &trackster_barycenter_y);
    trackster_tree_->Branch("barycenter_z", &trackster_barycenter_z);
    trackster_tree_->Branch("barycenter_eta", &trackster_barycenter_eta);
    trackster_tree_->Branch("barycenter_phi", &trackster_barycenter_phi);
    trackster_tree_->Branch("EV1", &trackster_EV1);
    trackster_tree_->Branch("EV2", &trackster_EV2);
    trackster_tree_->Branch("EV3", &trackster_EV3);
    trackster_tree_->Branch("eVector0_x", &trackster_eVector0_x);
    trackster_tree_->Branch("eVector0_y", &trackster_eVector0_y);
    trackster_tree_->Branch("eVector0_z", &trackster_eVector0_z);
    trackster_tree_->Branch("sigmaPCA1", &trackster_sigmaPCA1);
    trackster_tree_->Branch("sigmaPCA2", &trackster_sigmaPCA2);
    trackster_tree_->Branch("sigmaPCA3", &trackster_sigmaPCA3);
    trackster_tree_->Branch("id_probabilities", &trackster_id_probabilities);
    trackster_tree_->Branch("vertices_indexes", &trackster_vertices_indexes);
    trackster_tree_->Branch("vertices_x", &trackster_vertices_x);
    trackster_tree_->Branch("vertices_y", &trackster_vertices_y);
    trackster_tree_->Branch("vertices_z", &trackster_vertices_z);
    trackster_tree_->Branch("vertices_time", &trackster_vertices_time);
    trackster_tree_->Branch("vertices_timeErr", &trackster_vertices_timeErr);
    trackster_tree_->Branch("vertices_energy", &trackster_vertices_energy);
    trackster_tree_->Branch("vertices_correctedEnergy", &trackster_vertices_correctedEnergy);
    trackster_tree_->Branch("vertices_correctedEnergyUncertainty", &trackster_vertices_correctedEnergyUncertainty);
    trackster_tree_->Branch("vertices_multiplicity", &trackster_vertices_multiplicity);
  }

  void clearVariables() {
    nTracksters = 0;
    nClusters = 0;
    trackster_time.clear();
    trackster_timeError.clear();
    trackster_regressed_energy.clear();
    trackster_raw_energy.clear();
    trackster_raw_em_energy.clear();
    trackster_raw_pt.clear();
    trackster_raw_em_pt.clear();
    trackster_barycenter_x.clear();
    trackster_barycenter_y.clear();
    trackster_barycenter_z.clear();
    trackster_EV1.clear();
    trackster_EV2.clear();
    trackster_EV3.clear();
    trackster_eVector0_x.clear();
    trackster_eVector0_y.clear();
    trackster_eVector0_z.clear();
    trackster_sigmaPCA1.clear();
    trackster_sigmaPCA2.clear();
    trackster_sigmaPCA3.clear();
    trackster_barycenter_eta.clear();
    trackster_barycenter_phi.clear();
    trackster_id_probabilities.clear();
    trackster_vertices_indexes.clear();
    trackster_vertices_x.clear();
    trackster_vertices_y.clear();
    trackster_vertices_z.clear();
    trackster_vertices_time.clear();
    trackster_vertices_timeErr.clear();
    trackster_vertices_energy.clear();
    trackster_vertices_correctedEnergy.clear();
    trackster_vertices_correctedEnergyUncertainty.clear();
    trackster_vertices_multiplicity.clear();
  }

  void fillFromEvent(std::vector<ticl::Trackster> const& tracksters, std::vector<reco::CaloCluster> const& clusters,
      edm::ValueMap<std::pair<float, float>>const& layerClustersTimes) {
    nTracksters = tracksters.size();
    nClusters = clusters.size();
    for (auto trackster_iterator = tracksters.begin(); trackster_iterator != tracksters.end(); ++trackster_iterator) {
      //per-trackster analysis
      trackster_time.push_back(trackster_iterator->time());
      trackster_timeError.push_back(trackster_iterator->timeError());
      trackster_regressed_energy.push_back(trackster_iterator->regressed_energy());
      trackster_raw_energy.push_back(trackster_iterator->raw_energy());
      trackster_raw_em_energy.push_back(trackster_iterator->raw_em_energy());
      trackster_raw_pt.push_back(trackster_iterator->raw_pt());
      trackster_raw_em_pt.push_back(trackster_iterator->raw_em_pt());
      trackster_barycenter_x.push_back(trackster_iterator->barycenter().x());
      trackster_barycenter_y.push_back(trackster_iterator->barycenter().y());
      trackster_barycenter_z.push_back(trackster_iterator->barycenter().z());
      trackster_barycenter_eta.push_back(trackster_iterator->barycenter().eta());
      trackster_barycenter_phi.push_back(trackster_iterator->barycenter().phi());
      trackster_EV1.push_back(trackster_iterator->eigenvalues()[0]);
      trackster_EV2.push_back(trackster_iterator->eigenvalues()[1]);
      trackster_EV3.push_back(trackster_iterator->eigenvalues()[2]);
      trackster_eVector0_x.push_back((trackster_iterator->eigenvectors()[0]).x());
      trackster_eVector0_y.push_back((trackster_iterator->eigenvectors()[0]).y());
      trackster_eVector0_z.push_back((trackster_iterator->eigenvectors()[0]).z());
      trackster_sigmaPCA1.push_back(trackster_iterator->sigmasPCA()[0]);
      trackster_sigmaPCA2.push_back(trackster_iterator->sigmasPCA()[1]);
      trackster_sigmaPCA3.push_back(trackster_iterator->sigmasPCA()[2]);
      std::vector<float> id_probs;
      for (size_t i = 0; i < 8; i++)
        id_probs.push_back(trackster_iterator->id_probabilities(i));
      trackster_id_probabilities.push_back(id_probs);

      // Clusters
      std::vector<uint32_t> vertices_indexes;
      std::vector<float> vertices_x;
      std::vector<float> vertices_y;
      std::vector<float> vertices_z;
      std::vector<float> vertices_time;
      std::vector<float> vertices_timeErr;
      std::vector<float> vertices_energy;
      std::vector<float> vertices_correctedEnergy;
      std::vector<float> vertices_correctedEnergyUncertainty;
      for (auto idx : trackster_iterator->vertices()) {
        vertices_indexes.push_back(idx);
        auto associated_cluster = clusters[idx];
        vertices_x.push_back(associated_cluster.x());
        vertices_y.push_back(associated_cluster.y());
        vertices_z.push_back(associated_cluster.z());
        vertices_energy.push_back(associated_cluster.energy());
        vertices_correctedEnergy.push_back(associated_cluster.correctedEnergy());
        vertices_correctedEnergyUncertainty.push_back(associated_cluster.correctedEnergyUncertainty());
        vertices_time.push_back(layerClustersTimes.get(idx).first);
        vertices_timeErr.push_back(layerClustersTimes.get(idx).second);
      }
      trackster_vertices_indexes.push_back(vertices_indexes);
      trackster_vertices_x.push_back(vertices_x);
      trackster_vertices_y.push_back(vertices_y);
      trackster_vertices_z.push_back(vertices_z);
      trackster_vertices_time.push_back(vertices_time);
      trackster_vertices_timeErr.push_back(vertices_timeErr);
      trackster_vertices_energy.push_back(vertices_energy);
      trackster_vertices_correctedEnergy.push_back(vertices_correctedEnergy);
      trackster_vertices_correctedEnergyUncertainty.push_back(vertices_correctedEnergyUncertainty);

      // Multiplicity
      std::vector<float> vertices_multiplicity;
      for (auto multiplicity : trackster_iterator->vertex_multiplicity()) {
        vertices_multiplicity.push_back(multiplicity);
      }
      trackster_vertices_multiplicity.push_back(vertices_multiplicity);
    }
  }
private:
  unsigned int nTracksters;
  unsigned int nClusters;
  std::vector<float> trackster_time;
  std::vector<float> trackster_timeError;
  std::vector<float> trackster_regressed_energy;
  std::vector<float> trackster_raw_energy;
  std::vector<float> trackster_raw_em_energy;
  std::vector<float> trackster_raw_pt;
  std::vector<float> trackster_raw_em_pt;
  std::vector<float> trackster_barycenter_x;
  std::vector<float> trackster_barycenter_y;
  std::vector<float> trackster_barycenter_z;
  std::vector<float> trackster_EV1;
  std::vector<float> trackster_EV2;
  std::vector<float> trackster_EV3;
  std::vector<float> trackster_eVector0_x;
  std::vector<float> trackster_eVector0_y;
  std::vector<float> trackster_eVector0_z;
  std::vector<float> trackster_sigmaPCA1;
  std::vector<float> trackster_sigmaPCA2;
  std::vector<float> trackster_sigmaPCA3;
  std::vector<float> trackster_barycenter_eta;
  std::vector<float> trackster_barycenter_phi;
  std::vector<std::vector<float>> trackster_id_probabilities;
  std::vector<std::vector<uint32_t>> trackster_vertices_indexes;
  std::vector<std::vector<float>> trackster_vertices_x;
  std::vector<std::vector<float>> trackster_vertices_y;
  std::vector<std::vector<float>> trackster_vertices_z;
  std::vector<std::vector<float>> trackster_vertices_time;
  std::vector<std::vector<float>> trackster_vertices_timeErr;
  std::vector<std::vector<float>> trackster_vertices_energy;
  std::vector<std::vector<float>> trackster_vertices_correctedEnergy;
  std::vector<std::vector<float>> trackster_vertices_correctedEnergyUncertainty;
  std::vector<std::vector<float>> trackster_vertices_multiplicity;
};

class TracksterToSimTracksterAssociationHelper {
public:
  /**
   * TO be called once after tree creation. branchPrefix : for example tsCLUE3D. branchSuffix : SC or CP
  */
  void initTree(TTree* tree, std::string branchPrefix, std::string branchSuffix) {
    tree->Branch((branchPrefix + "_recoToSim_" + branchSuffix).c_str(), &recoToSim);
    tree->Branch((branchPrefix + "_recoToSim_" + branchSuffix + "_score").c_str(), &recoToSim_score);
    tree->Branch((branchPrefix + "_recoToSim_" + branchSuffix + "_sharedE").c_str(), &recoToSim_sharedE);
    tree->Branch((branchPrefix + "_simToReco_" + branchSuffix).c_str(), &simToReco);
    tree->Branch((branchPrefix + "_simToReco_" + branchSuffix + "_score").c_str(), &simToReco_score);
    tree->Branch((branchPrefix + "_simToReco_" + branchSuffix + "_sharedE").c_str(), &simToReco_sharedE);
  }

  void clearVariables() {
    recoToSim.clear();
    recoToSim_score.clear();
    recoToSim_sharedE.clear();
    simToReco.clear();
    simToReco_score.clear();
    simToReco_sharedE.clear();
  }

  void fillFromEvent(edm::Handle<std::vector<ticl::Trackster>> tracksters_handle, edm::Handle<std::vector<ticl::Trackster>> simTracksters_h, hgcal::RecoToSimCollectionSimTracksters const& tsRecoSimSCMap, hgcal::SimToRecoCollectionSimTracksters const& tsSimToRecoSCMap) {
    auto const& tracksters = *tracksters_handle;
    auto simTracksters = *simTracksters_h;

    // Reco -> Sim
    recoToSim.resize(tracksters.size());
    recoToSim_score.resize(tracksters.size());
    recoToSim_sharedE.resize(tracksters.size());
    for (size_t i = 0; i < tracksters.size(); ++i) {
      const edm::Ref<ticl::TracksterCollection> tsRef(tracksters_handle, i);

      const auto stsSC_iter = tsRecoSimSCMap.find(tsRef);
      if (stsSC_iter != tsRecoSimSCMap.end()) {
        const auto& stsSCassociated = stsSC_iter->val;
        for (auto& sts : stsSCassociated) {
          auto sts_id = (sts.first).get() - (edm::Ref<ticl::TracksterCollection>(simTracksters_h, 0)).get();
          recoToSim[i].push_back(sts_id);
          recoToSim_score[i].push_back(sts.second.second);
          recoToSim_sharedE[i].push_back(sts.second.first);
        }
      }
    }

    // Sim -> Reco
    simToReco.resize(simTracksters.size());
    simToReco_score.resize(simTracksters.size());
    simToReco_sharedE.resize(simTracksters.size());
    for (size_t i = 0; i < simTracksters.size(); ++i) {
      const edm::Ref<ticl::TracksterCollection> stsSCRef(simTracksters_h, i);

      // STS-SC -> CLUE3D
      const auto ts_iter = tsSimToRecoSCMap.find(stsSCRef);
      if (ts_iter != tsSimToRecoSCMap.end()) {
        const auto& tsAssociated = ts_iter->val;
        for (auto& ts : tsAssociated) {
          auto ts_idx = (ts.first).get() - (edm::Ref<ticl::TracksterCollection>(tracksters_handle, 0)).get();
          simToReco[i].push_back(ts_idx);
          simToReco_score[i].push_back(ts.second.second);
          simToReco_sharedE[i].push_back(ts.second.first);
        }
      }
    }
  }
private:
  std::vector<std::vector<uint32_t>> recoToSim;
  std::vector<std::vector<float>> recoToSim_score;
  std::vector<std::vector<float>> recoToSim_sharedE;
  std::vector<std::vector<uint32_t>> simToReco;
  std::vector<std::vector<float>> simToReco_score;
  std::vector<std::vector<float>> simToReco_sharedE;
};

class TICLDumper : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
public:
  explicit TICLDumper(const edm::ParameterSet&);
  ~TICLDumper() override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  typedef ticl::Vector Vector;
  typedef std::vector<double> Vec;

private:
  void beginJob() override;
  void beginRun(const edm::Run&, const edm::EventSetup&) override;

  void initialize(const HGCalDDDConstants* hgcons,
                  const hgcal::RecHitTools rhtools,
                  const edm::ESHandle<MagneticField> bfieldH,
                  const edm::ESHandle<Propagator> propH);
  void buildLayers();

  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override{};
  void endJob() override;

  // Define Tokens
  const std::vector<edm::ParameterSet> tracksters_parameterSets_; ///< A parameter set for each trackster collection to dump
  std::vector<edm::EDGetTokenT<std::vector<ticl::Trackster>>> tracksters_token_; ///< a token for each trackster collection to dump
  std::vector<std::vector<ticl::Trackster>> tracksters_collection; ///< the trackster collection for each trackster collection to dump
  std::vector<TracksterDumperHelper> tracksters_dumperHelpers_; ///< the dumper helpers for each trackster collection to dump
  std::vector<TTree*> tracksters_trees; ///< TTree for each trackster collection to dump

  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> layer_clusters_token_;
  const edm::EDGetTokenT<std::vector<TICLCandidate>> ticl_candidates_token_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> ticl_candidates_tracksters_token_; ///< trackster collection used by TICLCandidate
  const edm::EDGetTokenT<std::vector<reco::Track>> tracks_token_;
  const edm::EDGetTokenT<std::vector<bool>> tracks_mask_token_;
  const edm::EDGetTokenT<edm::ValueMap<float>> tracks_time_token_;
  const edm::EDGetTokenT<edm::ValueMap<float>> tracks_time_quality_token_;
  const edm::EDGetTokenT<edm::ValueMap<float>> tracks_time_err_token_;
  const edm::EDGetTokenT<edm::ValueMap<float>> tracks_beta_token_;
  const edm::EDGetTokenT<edm::ValueMap<float>> tracks_time_mtd_token_;
  const edm::EDGetTokenT<edm::ValueMap<float>> tracks_time_mtd_err_token_;
  const edm::EDGetTokenT<edm::ValueMap<GlobalPoint>> tracks_pos_mtd_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_x_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_y_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_z_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_eta_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_phi_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_px_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_py_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_pz_token_;
  const edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> clustersTime_token_;
  const edm::EDGetTokenT<std::vector<int>> tracksterSeeds_token_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> superclustering_superclusteredTracksters_token;
  const edm::EDGetTokenT<std::vector<std::vector<unsigned int>>> superclustering_linkedResultTracksters_token;
  //const edm::EDGetTokenT<ticl::SuperclusteringDNNScore> superclustering_DNNScore_token_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometry_token_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> simTracksters_SC_token_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> simTracksters_CP_token_;
  const edm::EDGetTokenT<std::vector<TICLCandidate>> simTICLCandidate_token_;

  // associators
  const std::vector<edm::ParameterSet> associations_parameterSets_; ///< A parameter set for each associator collection to dump
  std::vector<edm::EDGetTokenT<hgcal::SimToRecoCollectionSimTracksters>> associations_simToReco_token_;
  std::vector<edm::EDGetTokenT<hgcal::RecoToSimCollectionSimTracksters>> associations_recoToSim_token_;
  std::vector<TracksterToSimTracksterAssociationHelper> associations_dumperHelpers_; ///< the dumper helpers for each trackster collection to dump
  std::vector<edm::EDGetTokenT<std::vector<ticl::Trackster>>> associations_tracksterCollection_; ///< the collection of tracksters used by the associator
  std::vector<edm::EDGetTokenT<std::vector<ticl::Trackster>>> associations_simTracksterCollection_; ///< the collection of simtracksters used by the associator
  TTree* associations_tree_;

  const edm::EDGetTokenT<std::vector<SimCluster>> simclusters_token_;
  const edm::EDGetTokenT<std::vector<CaloParticle>> caloparticles_token_;

  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometry_token_;
  const std::string detector_;
  const std::string propName_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bfield_token_;
  const edm::ESGetToken<Propagator, TrackingComponentsRecord> propagator_token_;
  hgcal::RecHitTools rhtools_;
  edm::ESGetToken<HGCalDDDConstants, IdealGeometryRecord> hdc_token_;
  const HGCalDDDConstants* hgcons_;
  std::unique_ptr<GeomDet> firstDisk_[2];
  std::unique_ptr<GeomDet> interfaceDisk_[2];
  edm::ESHandle<MagneticField> bfield_;
  edm::ESHandle<Propagator> propagator_;
  bool saveLCs_;
  bool saveSuperclustering_;
  bool saveSuperclusteringDNNScore_;
  bool saveSimTrackstersSC_;
  bool saveSimTrackstersCP_;
  bool saveTICLCandidate_;
  bool saveSimTICLCandidate_;
  bool saveTracks_;

  // Output tree
  TTree* tree_;

  void clearVariables();

  // Variables for branches
  unsigned int ev_event_;
  unsigned int nclusters_;
  unsigned int stsSC_ntracksters_;
  unsigned int stsCP_ntracksters_;

  std::vector<std::vector<unsigned int>> superclustering_linkedResultTracksters; // Map of indices from superclusteredTracksters collection back into ticlTrackstersCLUE3DEM collection
  //ticl::SuperclusteringDNNScore superclusteringDNNScore;

  std::vector<float> stsSC_trackster_time;
  std::vector<float> stsSC_trackster_timeBoundary;
  std::vector<float> stsSC_trackster_timeError;
  std::vector<float> stsSC_trackster_regressed_energy;
  std::vector<float> stsSC_trackster_regressed_pt;
  std::vector<float> stsSC_trackster_raw_energy;
  std::vector<float> stsSC_trackster_raw_em_energy;
  std::vector<float> stsSC_trackster_raw_pt;
  std::vector<float> stsSC_trackster_raw_em_pt;
  std::vector<float> stsSC_trackster_barycenter_x;
  std::vector<float> stsSC_trackster_barycenter_y;
  std::vector<float> stsSC_trackster_barycenter_z;
  std::vector<float> stsSC_trackster_barycenter_eta;
  std::vector<float> stsSC_trackster_barycenter_phi;
  std::vector<float> stsSC_trackster_EV1;
  std::vector<float> stsSC_trackster_EV2;
  std::vector<float> stsSC_trackster_EV3;
  std::vector<float> stsSC_trackster_eVector0_x;
  std::vector<float> stsSC_trackster_eVector0_y;
  std::vector<float> stsSC_trackster_eVector0_z;
  std::vector<float> stsSC_trackster_sigmaPCA1;
  std::vector<float> stsSC_trackster_sigmaPCA2;
  std::vector<float> stsSC_trackster_sigmaPCA3;
  std::vector<int> stsSC_pdgID;
  std::vector<int> stsSC_trackIdx;
  std::vector<float> stsSC_trackTime;
  std::vector<float> stsSC_boundaryX;
  std::vector<float> stsSC_boundaryY;
  std::vector<float> stsSC_boundaryZ;
  std::vector<float> stsSC_boundaryEta;
  std::vector<float> stsSC_boundaryPhi;
  std::vector<float> stsSC_boundaryPx;
  std::vector<float> stsSC_boundaryPy;
  std::vector<float> stsSC_boundaryPz;
  std::vector<float> stsSC_track_boundaryX;
  std::vector<float> stsSC_track_boundaryY;
  std::vector<float> stsSC_track_boundaryZ;
  std::vector<float> stsSC_track_boundaryEta;
  std::vector<float> stsSC_track_boundaryPhi;
  std::vector<float> stsSC_track_boundaryPx;
  std::vector<float> stsSC_track_boundaryPy;
  std::vector<float> stsSC_track_boundaryPz;
  std::vector<std::vector<float>> stsSC_trackster_id_probabilities;
  std::vector<std::vector<uint32_t>> stsSC_trackster_vertices_indexes;
  std::vector<std::vector<float>> stsSC_trackster_vertices_x;
  std::vector<std::vector<float>> stsSC_trackster_vertices_y;
  std::vector<std::vector<float>> stsSC_trackster_vertices_z;
  std::vector<std::vector<float>> stsSC_trackster_vertices_time;
  std::vector<std::vector<float>> stsSC_trackster_vertices_timeErr;
  std::vector<std::vector<float>> stsSC_trackster_vertices_energy;
  std::vector<std::vector<float>> stsSC_trackster_vertices_correctedEnergy;
  std::vector<std::vector<float>> stsSC_trackster_vertices_correctedEnergyUncertainty;
  std::vector<std::vector<float>> stsSC_trackster_vertices_multiplicity;
  std::vector<float> stsCP_trackster_time;
  std::vector<float> stsCP_trackster_timeBoundary;
  std::vector<float> stsCP_trackster_timeError;
  std::vector<float> stsCP_trackster_regressed_energy;
  std::vector<float> stsCP_trackster_regressed_pt;
  std::vector<float> stsCP_trackster_raw_energy;
  std::vector<float> stsCP_trackster_raw_em_energy;
  std::vector<float> stsCP_trackster_raw_pt;
  std::vector<float> stsCP_trackster_raw_em_pt;
  std::vector<float> stsCP_trackster_barycenter_x;
  std::vector<float> stsCP_trackster_barycenter_y;
  std::vector<float> stsCP_trackster_barycenter_z;
  std::vector<float> stsCP_trackster_barycenter_eta;
  std::vector<float> stsCP_trackster_barycenter_phi;
  std::vector<float> stsCP_trackster_EV1;
  std::vector<float> stsCP_trackster_EV2;
  std::vector<float> stsCP_trackster_EV3;
  std::vector<float> stsCP_trackster_eVector0_x;
  std::vector<float> stsCP_trackster_eVector0_y;
  std::vector<float> stsCP_trackster_eVector0_z;
  std::vector<float> stsCP_trackster_sigmaPCA1;
  std::vector<float> stsCP_trackster_sigmaPCA2;
  std::vector<float> stsCP_trackster_sigmaPCA3;
  std::vector<int> stsCP_pdgID;
  std::vector<int> stsCP_trackIdx;
  std::vector<float> stsCP_trackTime;
  std::vector<float> stsCP_boundaryX;
  std::vector<float> stsCP_boundaryY;
  std::vector<float> stsCP_boundaryZ;
  std::vector<float> stsCP_boundaryEta;
  std::vector<float> stsCP_boundaryPhi;
  std::vector<float> stsCP_boundaryPx;
  std::vector<float> stsCP_boundaryPy;
  std::vector<float> stsCP_boundaryPz;
  std::vector<float> stsCP_track_boundaryX;
  std::vector<float> stsCP_track_boundaryY;
  std::vector<float> stsCP_track_boundaryZ;
  std::vector<float> stsCP_track_boundaryEta;
  std::vector<float> stsCP_track_boundaryPhi;
  std::vector<float> stsCP_track_boundaryPx;
  std::vector<float> stsCP_track_boundaryPy;
  std::vector<float> stsCP_track_boundaryPz;
  std::vector<std::vector<float>> stsCP_trackster_id_probabilities;
  std::vector<std::vector<uint32_t>> stsCP_trackster_vertices_indexes;
  std::vector<std::vector<float>> stsCP_trackster_vertices_x;
  std::vector<std::vector<float>> stsCP_trackster_vertices_y;
  std::vector<std::vector<float>> stsCP_trackster_vertices_z;
  std::vector<std::vector<float>> stsCP_trackster_vertices_time;
  std::vector<std::vector<float>> stsCP_trackster_vertices_timeErr;
  std::vector<std::vector<float>> stsCP_trackster_vertices_energy;
  std::vector<std::vector<float>> stsCP_trackster_vertices_correctedEnergy;
  std::vector<std::vector<float>> stsCP_trackster_vertices_correctedEnergyUncertainty;
  std::vector<std::vector<float>> stsCP_trackster_vertices_multiplicity;

  std::vector<float> simTICLCandidate_raw_energy;
  std::vector<float> simTICLCandidate_regressed_energy;
  std::vector<std::vector<int>> simTICLCandidate_simTracksterCPIndex;
  std::vector<float> simTICLCandidate_boundaryX;
  std::vector<float> simTICLCandidate_boundaryY;
  std::vector<float> simTICLCandidate_boundaryZ;
  std::vector<float> simTICLCandidate_boundaryPx;
  std::vector<float> simTICLCandidate_boundaryPy;
  std::vector<float> simTICLCandidate_boundaryPz;
  std::vector<float> simTICLCandidate_caloParticleMass;
  std::vector<float> simTICLCandidate_time;
  std::vector<int> simTICLCandidate_pdgId;
  std::vector<int> simTICLCandidate_charge;
  std::vector<int> simTICLCandidate_track_in_candidate;

  // from TICLCandidate, product of linking
  size_t nCandidates;
  std::vector<int> candidate_charge;
  std::vector<int> candidate_pdgId;
  std::vector<float> candidate_energy;
  std::vector<double> candidate_px;
  std::vector<double> candidate_py;
  std::vector<double> candidate_pz;
  std::vector<float> candidate_time;
  std::vector<float> candidate_time_err;
  std::vector<std::vector<float>> candidate_id_probabilities;
  std::vector<std::vector<uint32_t>> tracksters_in_candidate;
  std::vector<int> track_in_candidate;

  std::vector<uint32_t> cluster_seedID;
  std::vector<float> cluster_energy;
  std::vector<float> cluster_correctedEnergy;
  std::vector<float> cluster_correctedEnergyUncertainty;
  std::vector<float> cluster_position_x;
  std::vector<float> cluster_position_y;
  std::vector<float> cluster_position_z;
  std::vector<float> cluster_position_eta;
  std::vector<float> cluster_position_phi;
  std::vector<unsigned int> cluster_layer_id;
  std::vector<int> cluster_type;
  std::vector<float> cluster_time;
  std::vector<float> cluster_timeErr;
  std::vector<uint32_t> cluster_number_of_hits;

  std::vector<unsigned int> track_id;
  std::vector<float> track_hgcal_x;
  std::vector<float> track_hgcal_y;
  std::vector<float> track_hgcal_z;
  std::vector<float> track_hgcal_px;
  std::vector<float> track_hgcal_py;
  std::vector<float> track_hgcal_pz;
  std::vector<float> track_hgcal_eta;
  std::vector<float> track_hgcal_phi;
  std::vector<float> track_hgcal_pt;
  std::vector<float> track_pt;
  std::vector<int> track_quality;
  std::vector<int> track_missing_outer_hits;
  std::vector<int> track_charge;
  std::vector<double> track_time;
  std::vector<float> track_time_quality;
  std::vector<float> track_time_err;
  std::vector<float> track_beta;
  std::vector<float> track_time_mtd;
  std::vector<float> track_time_mtd_err;
  std::vector<GlobalPoint> track_pos_mtd;
  std::vector<int> track_nhits;

  TTree* cluster_tree_;
  TTree* candidate_tree_;
  TTree* superclustering_tree_;
  TTree* simtrackstersSC_tree_;
  TTree* simtrackstersCP_tree_;
  TTree* tracks_tree_;
  TTree* simTICLCandidate_tree;
};

void TICLDumper::clearVariables() {
  // event info
  nclusters_ = 0;

  for (TracksterDumperHelper& tsDumper : tracksters_dumperHelpers_) {
    tsDumper.clearVariables();
  }

  //superclustering_superclusteredTracksters_dumper.clearVariables();
  superclustering_linkedResultTracksters.clear();
  // superclusteringDNNScore.clear();

  stsSC_trackster_time.clear();
  stsSC_trackster_timeBoundary.clear();
  stsSC_trackster_timeError.clear();
  stsSC_trackster_regressed_energy.clear();
  stsSC_trackster_regressed_pt.clear();
  stsSC_trackster_raw_energy.clear();
  stsSC_trackster_raw_em_energy.clear();
  stsSC_trackster_raw_pt.clear();
  stsSC_trackster_raw_em_pt.clear();
  stsSC_trackster_barycenter_x.clear();
  stsSC_trackster_barycenter_y.clear();
  stsSC_trackster_barycenter_z.clear();
  stsSC_trackster_EV1.clear();
  stsSC_trackster_EV2.clear();
  stsSC_trackster_EV3.clear();
  stsSC_trackster_eVector0_x.clear();
  stsSC_trackster_eVector0_y.clear();
  stsSC_trackster_eVector0_z.clear();
  stsSC_trackster_sigmaPCA1.clear();
  stsSC_trackster_sigmaPCA2.clear();
  stsSC_trackster_sigmaPCA3.clear();
  stsSC_trackster_barycenter_eta.clear();
  stsSC_trackster_barycenter_phi.clear();
  stsSC_pdgID.clear();
  stsSC_trackIdx.clear();
  stsSC_trackTime.clear();
  stsSC_boundaryX.clear();
  stsSC_boundaryY.clear();
  stsSC_boundaryZ.clear();
  stsSC_boundaryEta.clear();
  stsSC_boundaryPhi.clear();
  stsSC_boundaryPx.clear();
  stsSC_boundaryPy.clear();
  stsSC_boundaryPz.clear();
  stsSC_track_boundaryX.clear();
  stsSC_track_boundaryY.clear();
  stsSC_track_boundaryZ.clear();
  stsSC_track_boundaryEta.clear();
  stsSC_track_boundaryPhi.clear();
  stsSC_track_boundaryPx.clear();
  stsSC_track_boundaryPy.clear();
  stsSC_track_boundaryPz.clear();
  stsSC_trackster_id_probabilities.clear();
  stsSC_trackster_vertices_indexes.clear();
  stsSC_trackster_vertices_x.clear();
  stsSC_trackster_vertices_y.clear();
  stsSC_trackster_vertices_z.clear();
  stsSC_trackster_vertices_time.clear();
  stsSC_trackster_vertices_timeErr.clear();
  stsSC_trackster_vertices_energy.clear();
  stsSC_trackster_vertices_correctedEnergy.clear();
  stsSC_trackster_vertices_correctedEnergyUncertainty.clear();
  stsSC_trackster_vertices_multiplicity.clear();

  stsCP_trackster_time.clear();
  stsCP_trackster_timeBoundary.clear();
  stsCP_trackster_timeError.clear();
  stsCP_trackster_regressed_energy.clear();
  stsCP_trackster_regressed_pt.clear();
  stsCP_trackster_raw_energy.clear();
  stsCP_trackster_raw_em_energy.clear();
  stsCP_trackster_raw_pt.clear();
  stsCP_trackster_raw_em_pt.clear();
  stsCP_trackster_barycenter_x.clear();
  stsCP_trackster_barycenter_y.clear();
  stsCP_trackster_barycenter_z.clear();
  stsCP_trackster_sigmaPCA1.clear();
  stsCP_trackster_sigmaPCA2.clear();
  stsCP_trackster_sigmaPCA3.clear();
  stsCP_trackster_barycenter_eta.clear();
  stsCP_trackster_barycenter_phi.clear();
  stsCP_pdgID.clear();
  stsCP_trackIdx.clear();
  stsCP_trackTime.clear();
  stsCP_boundaryX.clear();
  stsCP_boundaryY.clear();
  stsCP_boundaryZ.clear();
  stsCP_boundaryEta.clear();
  stsCP_boundaryPhi.clear();
  stsCP_boundaryPx.clear();
  stsCP_boundaryPy.clear();
  stsCP_boundaryPz.clear();
  stsCP_track_boundaryX.clear();
  stsCP_track_boundaryY.clear();
  stsCP_track_boundaryZ.clear();
  stsCP_track_boundaryEta.clear();
  stsCP_track_boundaryPhi.clear();
  stsCP_track_boundaryPx.clear();
  stsCP_track_boundaryPy.clear();
  stsCP_track_boundaryPz.clear();
  stsCP_trackster_id_probabilities.clear();
  stsCP_trackster_vertices_indexes.clear();
  stsCP_trackster_vertices_x.clear();
  stsCP_trackster_vertices_y.clear();
  stsCP_trackster_vertices_z.clear();
  stsCP_trackster_vertices_time.clear();
  stsCP_trackster_vertices_timeErr.clear();
  stsCP_trackster_vertices_energy.clear();
  stsCP_trackster_vertices_correctedEnergy.clear();
  stsCP_trackster_vertices_correctedEnergyUncertainty.clear();
  stsCP_trackster_vertices_multiplicity.clear();

  simTICLCandidate_raw_energy.clear();
  simTICLCandidate_regressed_energy.clear();
  simTICLCandidate_simTracksterCPIndex.clear();
  simTICLCandidate_boundaryX.clear();
  simTICLCandidate_boundaryY.clear();
  simTICLCandidate_boundaryZ.clear();
  simTICLCandidate_boundaryPx.clear();
  simTICLCandidate_boundaryPy.clear();
  simTICLCandidate_boundaryPz.clear();
  simTICLCandidate_time.clear();
  simTICLCandidate_caloParticleMass.clear();
  simTICLCandidate_pdgId.clear();
  simTICLCandidate_charge.clear();
  simTICLCandidate_track_in_candidate.clear();

  nCandidates = 0;
  candidate_charge.clear();
  candidate_pdgId.clear();
  candidate_energy.clear();
  candidate_px.clear();
  candidate_py.clear();
  candidate_pz.clear();
  candidate_time.clear();
  candidate_time_err.clear();
  candidate_id_probabilities.clear();
  tracksters_in_candidate.clear();
  track_in_candidate.clear();

  for (auto& helper : associations_dumperHelpers_) {
    helper.clearVariables();
  }

  cluster_seedID.clear();
  cluster_energy.clear();
  cluster_correctedEnergy.clear();
  cluster_correctedEnergyUncertainty.clear();
  cluster_position_x.clear();
  cluster_position_y.clear();
  cluster_position_z.clear();
  cluster_position_eta.clear();
  cluster_position_phi.clear();
  cluster_layer_id.clear();
  cluster_type.clear();
  cluster_time.clear();
  cluster_timeErr.clear();
  cluster_number_of_hits.clear();

  track_id.clear();
  track_hgcal_x.clear();
  track_hgcal_y.clear();
  track_hgcal_z.clear();
  track_hgcal_eta.clear();
  track_hgcal_phi.clear();
  track_hgcal_px.clear();
  track_hgcal_py.clear();
  track_hgcal_pz.clear();
  track_hgcal_pt.clear();
  track_quality.clear();
  track_pt.clear();
  track_missing_outer_hits.clear();
  track_charge.clear();
  track_time.clear();
  track_time_quality.clear();
  track_time_err.clear();
  track_beta.clear();
  track_time_mtd.clear();
  track_time_mtd_err.clear();
  track_pos_mtd.clear();
  track_nhits.clear();
};

TICLDumper::TICLDumper(const edm::ParameterSet& ps)
    : tracksters_parameterSets_(ps.getParameter<std::vector<edm::ParameterSet>>("tracksterCollections")),
      tracksters_token_(),
      layer_clusters_token_(consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layerClusters"))),
      ticl_candidates_token_(consumes<std::vector<TICLCandidate>>(ps.getParameter<edm::InputTag>("ticlcandidates"))),
      ticl_candidates_tracksters_token_(consumes<std::vector<ticl::Trackster>>(ps.getParameter<edm::InputTag>("ticlcandidates"))),
      tracks_token_(consumes<std::vector<reco::Track>>(ps.getParameter<edm::InputTag>("tracks"))),
      tracks_time_token_(consumes<edm::ValueMap<float>>(ps.getParameter<edm::InputTag>("tracksTime"))),
      tracks_time_quality_token_(consumes<edm::ValueMap<float>>(ps.getParameter<edm::InputTag>("tracksTimeQual"))),
      tracks_time_err_token_(consumes<edm::ValueMap<float>>(ps.getParameter<edm::InputTag>("tracksTimeErr"))),
      tracks_beta_token_(consumes<edm::ValueMap<float>>(ps.getParameter<edm::InputTag>("tracksBeta"))),
      tracks_time_mtd_token_(consumes<edm::ValueMap<float>>(ps.getParameter<edm::InputTag>("tracksTimeMtd"))),
      tracks_time_mtd_err_token_(consumes<edm::ValueMap<float>>(ps.getParameter<edm::InputTag>("tracksTimeMtdErr"))),
      tracks_pos_mtd_token_(consumes<edm::ValueMap<GlobalPoint>>(ps.getParameter<edm::InputTag>("tracksPosMtd"))),
      clustersTime_token_(
          consumes<edm::ValueMap<std::pair<float, float>>>(ps.getParameter<edm::InputTag>("layer_clustersTime"))),
      superclustering_superclusteredTracksters_token(
          consumes<std::vector<ticl::Trackster>>(ps.getParameter<edm::InputTag>("superclustering"))),
      superclustering_linkedResultTracksters_token(
          consumes<std::vector<std::vector<unsigned int>>>(ps.getParameter<edm::InputTag>("superclustering"))),
      //superclustering_DNNScore_token_(consumes<ticl::SuperclusteringDNNScore>(
      //    ps.getParameter<edm::InputTag>("superclusteringDNNScore"))),
      caloGeometry_token_(esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>()),
      simTracksters_SC_token_(
          consumes<std::vector<ticl::Trackster>>(ps.getParameter<edm::InputTag>("simtrackstersSC"))),
      simTracksters_CP_token_(
          consumes<std::vector<ticl::Trackster>>(ps.getParameter<edm::InputTag>("simtrackstersCP"))),
      simTICLCandidate_token_(
          consumes<std::vector<TICLCandidate>>(ps.getParameter<edm::InputTag>("simTICLCandidates"))),
      associations_parameterSets_(ps.getParameter<std::vector<edm::ParameterSet>>("associators")),
      simclusters_token_(consumes(ps.getParameter<edm::InputTag>("simclusters"))),
      caloparticles_token_(consumes(ps.getParameter<edm::InputTag>("caloparticles"))),
      geometry_token_(esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>()),
      detector_(ps.getParameter<std::string>("detector")),
      propName_(ps.getParameter<std::string>("propagator")),
      bfield_token_(esConsumes<MagneticField, IdealMagneticFieldRecord, edm::Transition::BeginRun>()),
      propagator_token_(
          esConsumes<Propagator, TrackingComponentsRecord, edm::Transition::BeginRun>(edm::ESInputTag("", propName_))),
      saveLCs_(ps.getParameter<bool>("saveLCs")),
      saveSuperclustering_(ps.getParameter<bool>("saveSuperclustering")),
      //saveSuperclusteringDNNScore_(ps.getParameter<bool>("saveSuperclusteringDNNScore")),
      saveSimTrackstersSC_(ps.getParameter<bool>("saveSimTrackstersSC")),
      saveSimTrackstersCP_(ps.getParameter<bool>("saveSimTrackstersCP")),
      saveTICLCandidate_(ps.getParameter<bool>("saveSimTICLCandidate")),
      saveSimTICLCandidate_(ps.getParameter<bool>("saveSimTICLCandidate")),
      saveTracks_(ps.getParameter<bool>("saveTracks")) {
  std::string detectorName_ = (detector_ == "HFNose") ? "HGCalHFNoseSensitive" : "HGCalEESensitive";
  hdc_token_ =
      esConsumes<HGCalDDDConstants, IdealGeometryRecord, edm::Transition::BeginRun>(edm::ESInputTag("", detectorName_));
  
  for (edm::ParameterSet const& tracksterPset : tracksters_parameterSets_) {
    tracksters_token_.push_back(consumes<std::vector<ticl::Trackster>>(tracksterPset.getParameter<edm::InputTag>("tracksterInputTag")));
  }

  for (edm::ParameterSet const& associationPset : associations_parameterSets_) {
    edm::InputTag const& inputTag = associationPset.getParameter<edm::InputTag>("associatorInputTag");
    associations_recoToSim_token_.push_back(consumes<hgcal::RecoToSimCollectionSimTracksters>(edm::InputTag(inputTag.label(), "recoToSim", inputTag.process())));
    associations_simToReco_token_.push_back(consumes<hgcal::SimToRecoCollectionSimTracksters>(edm::InputTag(inputTag.label(), "simToReco", inputTag.process())));
    associations_tracksterCollection_.push_back(consumes<std::vector<ticl::Trackster>>(associationPset.getParameter<edm::InputTag>("tracksterCollection")));
    associations_simTracksterCollection_.push_back(consumes<std::vector<ticl::Trackster>>(associationPset.getParameter<edm::InputTag>("simTracksterCollection")));
  }
};

TICLDumper::~TICLDumper() { clearVariables(); };

void TICLDumper::beginRun(edm::Run const&, edm::EventSetup const& es) {
  const CaloGeometry& geom = es.getData(caloGeometry_token_);
  rhtools_.setGeometry(geom);

  edm::ESHandle<HGCalDDDConstants> hdc = es.getHandle(hdc_token_);
  hgcons_ = hdc.product();
  edm::ESHandle<MagneticField> bfield_ = es.getHandle(bfield_token_);
  edm::ESHandle<Propagator> propagator = es.getHandle(propagator_token_);
  initialize(hgcons_, rhtools_, bfield_, propagator);
}

// Define tree and branches
void TICLDumper::beginJob() {
  edm::Service<TFileService> fs;
  for (unsigned int i = 0; i < tracksters_parameterSets_.size(); i++) {
    edm::ParameterSet const& tracksterPset = tracksters_parameterSets_[i];
    tracksters_dumperHelpers_.emplace_back();
    TTree* tree = fs->make<TTree>(tracksterPset.getParameter<std::string>("treeName").c_str(), ("Tracksters : " + tracksterPset.getParameter<std::string>("treeName") + 
       " (InputTag : " + tracksterPset.getParameter<edm::InputTag>("tracksterInputTag").encode() + ")").c_str());
    tracksters_trees.push_back(tree);
    tracksters_dumperHelpers_[i].initTree(tree, &ev_event_);
  }
  if (saveLCs_) {
    cluster_tree_ = fs->make<TTree>("clusters", "TICL tracksters");
    cluster_tree_->Branch("seedID", &cluster_seedID);
    cluster_tree_->Branch("energy", &cluster_energy);
    cluster_tree_->Branch("correctedEnergy", &cluster_correctedEnergy);
    cluster_tree_->Branch("correctedEnergyUncertainty", &cluster_correctedEnergyUncertainty);
    cluster_tree_->Branch("position_x", &cluster_position_x);
    cluster_tree_->Branch("position_y", &cluster_position_y);
    cluster_tree_->Branch("position_z", &cluster_position_z);
    cluster_tree_->Branch("position_eta", &cluster_position_eta);
    cluster_tree_->Branch("position_phi", &cluster_position_phi);
    cluster_tree_->Branch("cluster_layer_id", &cluster_layer_id);
    cluster_tree_->Branch("cluster_type", &cluster_type);
    cluster_tree_->Branch("cluster_time", &cluster_time);
    cluster_tree_->Branch("cluster_timeErr", &cluster_timeErr);
    cluster_tree_->Branch("cluster_number_of_hits", &cluster_number_of_hits);
  }
  if (saveTICLCandidate_) {
    candidate_tree_ = fs->make<TTree>("candidates", "TICL candidates");
    candidate_tree_->Branch("NCandidates", &nCandidates);
    candidate_tree_->Branch("candidate_charge", &candidate_charge);
    candidate_tree_->Branch("candidate_pdgId", &candidate_pdgId);
    candidate_tree_->Branch("candidate_id_probabilities", &candidate_id_probabilities);
    candidate_tree_->Branch("candidate_time", &candidate_time);
    candidate_tree_->Branch("candidate_timeErr", &candidate_time_err);
    candidate_tree_->Branch("candidate_energy", &candidate_energy);
    candidate_tree_->Branch("candidate_px", &candidate_px);
    candidate_tree_->Branch("candidate_py", &candidate_py);
    candidate_tree_->Branch("candidate_pz", &candidate_pz);
    candidate_tree_->Branch("track_in_candidate", &track_in_candidate);
    candidate_tree_->Branch("tracksters_in_candidate", &tracksters_in_candidate);
  }

  if (saveSuperclustering_) {
    superclustering_tree_ = fs->make<TTree>("superclustering", "Superclustering in HGCAL CE-E");
    //superclustering_superclusteredTracksters_dumper.initTree(superclustering_tree_, &ev_event_);
    superclustering_tree_->Branch("linkedResultTracksters", &superclustering_linkedResultTracksters);
  }
  // if (saveSuperclusteringDNNScore_) {
  //   superclustering_tree_->Branch("superclusteringDNNScore", &superclusteringDNNScore);
  // }

  if (associations_dumperHelpers_.size()>0)
    associations_tree_ = fs->make<TTree>("associations", "Associations");
  for (unsigned int i = 0; i < associations_parameterSets_.size(); i++) {
    associations_dumperHelpers_.emplace_back();
    associations_dumperHelpers_[i].initTree(associations_tree_, associations_parameterSets_[i].getParameter<std::string>("branchName"), associations_parameterSets_[i].getParameter<std::string>("suffix"));
  }

  if (saveSimTrackstersSC_) {
    simtrackstersSC_tree_ = fs->make<TTree>("simtrackstersSC", "TICL simTracksters SC");
    simtrackstersSC_tree_->Branch("event", &ev_event_);
    simtrackstersSC_tree_->Branch("NTracksters", &stsSC_ntracksters_);
    simtrackstersSC_tree_->Branch("time", &stsSC_trackster_time);
    simtrackstersSC_tree_->Branch("timeBoundary", &stsSC_trackster_timeBoundary);
    simtrackstersSC_tree_->Branch("timeError", &stsSC_trackster_timeError);
    simtrackstersSC_tree_->Branch("regressed_energy", &stsSC_trackster_regressed_energy);
    simtrackstersSC_tree_->Branch("regressed_pt", &stsSC_trackster_regressed_pt);
    simtrackstersSC_tree_->Branch("raw_energy", &stsSC_trackster_raw_energy);
    simtrackstersSC_tree_->Branch("raw_em_energy", &stsSC_trackster_raw_em_energy);
    simtrackstersSC_tree_->Branch("raw_pt", &stsSC_trackster_raw_pt);
    simtrackstersSC_tree_->Branch("raw_em_pt", &stsSC_trackster_raw_em_pt);
    simtrackstersSC_tree_->Branch("barycenter_x", &stsSC_trackster_barycenter_x);
    simtrackstersSC_tree_->Branch("barycenter_y", &stsSC_trackster_barycenter_y);
    simtrackstersSC_tree_->Branch("barycenter_z", &stsSC_trackster_barycenter_z);
    simtrackstersSC_tree_->Branch("barycenter_eta", &stsSC_trackster_barycenter_eta);
    simtrackstersSC_tree_->Branch("barycenter_phi", &stsSC_trackster_barycenter_phi);
    simtrackstersSC_tree_->Branch("EV1", &stsSC_trackster_EV1);
    simtrackstersSC_tree_->Branch("EV2", &stsSC_trackster_EV2);
    simtrackstersSC_tree_->Branch("EV3", &stsSC_trackster_EV3);
    simtrackstersSC_tree_->Branch("eVector0_x", &stsSC_trackster_eVector0_x);
    simtrackstersSC_tree_->Branch("eVector0_y", &stsSC_trackster_eVector0_y);
    simtrackstersSC_tree_->Branch("eVector0_z", &stsSC_trackster_eVector0_z);
    simtrackstersSC_tree_->Branch("sigmaPCA1", &stsSC_trackster_sigmaPCA1);
    simtrackstersSC_tree_->Branch("sigmaPCA2", &stsSC_trackster_sigmaPCA2);
    simtrackstersSC_tree_->Branch("sigmaPCA3", &stsSC_trackster_sigmaPCA3);
    simtrackstersSC_tree_->Branch("pdgID", &stsSC_pdgID);
    simtrackstersSC_tree_->Branch("trackIdx", &stsSC_trackIdx);
    simtrackstersSC_tree_->Branch("trackTime", &stsSC_trackTime);
    simtrackstersSC_tree_->Branch("boundaryX", &stsSC_boundaryX);
    simtrackstersSC_tree_->Branch("boundaryY", &stsSC_boundaryY);
    simtrackstersSC_tree_->Branch("boundaryZ", &stsSC_boundaryZ);
    simtrackstersSC_tree_->Branch("boundaryEta", &stsSC_boundaryEta);
    simtrackstersSC_tree_->Branch("boundaryPhi", &stsSC_boundaryPhi);
    simtrackstersSC_tree_->Branch("boundaryPx", &stsSC_boundaryPx);
    simtrackstersSC_tree_->Branch("boundaryPy", &stsSC_boundaryPy);
    simtrackstersSC_tree_->Branch("boundaryPz", &stsSC_boundaryPz);
    simtrackstersSC_tree_->Branch("track_boundaryX", &stsSC_track_boundaryX);
    simtrackstersSC_tree_->Branch("track_boundaryY", &stsSC_track_boundaryY);
    simtrackstersSC_tree_->Branch("track_boundaryZ", &stsSC_track_boundaryZ);
    simtrackstersSC_tree_->Branch("track_boundaryEta", &stsSC_track_boundaryEta);
    simtrackstersSC_tree_->Branch("track_boundaryPhi", &stsSC_track_boundaryPhi);
    simtrackstersSC_tree_->Branch("track_boundaryPx", &stsSC_track_boundaryPx);
    simtrackstersSC_tree_->Branch("track_boundaryPy", &stsSC_track_boundaryPy);
    simtrackstersSC_tree_->Branch("track_boundaryPz", &stsSC_track_boundaryPz);
    simtrackstersSC_tree_->Branch("id_probabilities", &stsSC_trackster_id_probabilities);
    simtrackstersSC_tree_->Branch("vertices_indexes", &stsSC_trackster_vertices_indexes);
    simtrackstersSC_tree_->Branch("vertices_x", &stsSC_trackster_vertices_x);
    simtrackstersSC_tree_->Branch("vertices_y", &stsSC_trackster_vertices_y);
    simtrackstersSC_tree_->Branch("vertices_z", &stsSC_trackster_vertices_z);
    simtrackstersSC_tree_->Branch("vertices_time", &stsSC_trackster_vertices_time);
    simtrackstersSC_tree_->Branch("vertices_timeErr", &stsSC_trackster_vertices_timeErr);
    simtrackstersSC_tree_->Branch("vertices_energy", &stsSC_trackster_vertices_energy);
    simtrackstersSC_tree_->Branch("vertices_correctedEnergy", &stsSC_trackster_vertices_correctedEnergy);
    simtrackstersSC_tree_->Branch("vertices_correctedEnergyUncertainty",
                                  &stsSC_trackster_vertices_correctedEnergyUncertainty);
    simtrackstersSC_tree_->Branch("vertices_multiplicity", &stsSC_trackster_vertices_multiplicity);
    simtrackstersSC_tree_->Branch("NsimTrackstersSC", &stsSC_ntracksters_); // same thing as NTracksters, not sure why it is duplicated
  }
  if (saveSimTrackstersCP_) {
    simtrackstersCP_tree_ = fs->make<TTree>("simtrackstersCP", "TICL simTracksters CP");
    simtrackstersCP_tree_->Branch("event", &ev_event_);
    simtrackstersCP_tree_->Branch("NTracksters", &stsCP_ntracksters_);
    simtrackstersCP_tree_->Branch("time", &stsCP_trackster_time);
    simtrackstersCP_tree_->Branch("timeBoundary", &stsCP_trackster_timeBoundary);
    simtrackstersCP_tree_->Branch("timeError", &stsCP_trackster_timeError);
    simtrackstersCP_tree_->Branch("regressed_energy", &stsCP_trackster_regressed_energy);
    simtrackstersCP_tree_->Branch("regressed_pt", &stsCP_trackster_regressed_pt);
    simtrackstersCP_tree_->Branch("raw_energy", &stsCP_trackster_raw_energy);
    simtrackstersCP_tree_->Branch("raw_em_energy", &stsCP_trackster_raw_em_energy);
    simtrackstersCP_tree_->Branch("raw_pt", &stsCP_trackster_raw_pt);
    simtrackstersCP_tree_->Branch("raw_em_pt", &stsCP_trackster_raw_em_pt);
    simtrackstersCP_tree_->Branch("barycenter_x", &stsCP_trackster_barycenter_x);
    simtrackstersCP_tree_->Branch("barycenter_y", &stsCP_trackster_barycenter_y);
    simtrackstersCP_tree_->Branch("barycenter_z", &stsCP_trackster_barycenter_z);
    simtrackstersCP_tree_->Branch("barycenter_eta", &stsCP_trackster_barycenter_eta);
    simtrackstersCP_tree_->Branch("barycenter_phi", &stsCP_trackster_barycenter_phi);
    simtrackstersCP_tree_->Branch("pdgID", &stsCP_pdgID);
    simtrackstersCP_tree_->Branch("trackIdx", &stsCP_trackIdx);
    simtrackstersCP_tree_->Branch("trackTime", &stsCP_trackTime);
    simtrackstersCP_tree_->Branch("boundaryX", &stsCP_boundaryX);
    simtrackstersCP_tree_->Branch("boundaryY", &stsCP_boundaryY);
    simtrackstersCP_tree_->Branch("boundaryZ", &stsCP_boundaryZ);
    simtrackstersCP_tree_->Branch("boundaryEta", &stsCP_boundaryEta);
    simtrackstersCP_tree_->Branch("boundaryPhi", &stsCP_boundaryPhi);
    simtrackstersCP_tree_->Branch("boundaryPx", &stsCP_boundaryPx);
    simtrackstersCP_tree_->Branch("boundaryPy", &stsCP_boundaryPy);
    simtrackstersCP_tree_->Branch("boundaryPz", &stsCP_boundaryPz);
    simtrackstersCP_tree_->Branch("track_boundaryX", &stsCP_track_boundaryX);
    simtrackstersCP_tree_->Branch("track_boundaryY", &stsCP_track_boundaryY);
    simtrackstersCP_tree_->Branch("track_boundaryZ", &stsCP_track_boundaryZ);
    simtrackstersCP_tree_->Branch("track_boundaryEta", &stsCP_track_boundaryEta);
    simtrackstersCP_tree_->Branch("track_boundaryPhi", &stsCP_track_boundaryPhi);
    simtrackstersCP_tree_->Branch("track_boundaryPx", &stsCP_track_boundaryPx);
    simtrackstersCP_tree_->Branch("track_boundaryPy", &stsCP_track_boundaryPy);
    simtrackstersCP_tree_->Branch("track_boundaryPz", &stsCP_track_boundaryPz);
    simtrackstersCP_tree_->Branch("EV1", &stsCP_trackster_EV1);
    simtrackstersCP_tree_->Branch("EV2", &stsCP_trackster_EV2);
    simtrackstersCP_tree_->Branch("EV3", &stsCP_trackster_EV3);
    simtrackstersCP_tree_->Branch("eVector0_x", &stsCP_trackster_eVector0_x);
    simtrackstersCP_tree_->Branch("eVector0_y", &stsCP_trackster_eVector0_y);
    simtrackstersCP_tree_->Branch("eVector0_z", &stsCP_trackster_eVector0_z);
    simtrackstersCP_tree_->Branch("sigmaPCA1", &stsCP_trackster_sigmaPCA1);
    simtrackstersCP_tree_->Branch("sigmaPCA2", &stsCP_trackster_sigmaPCA2);
    simtrackstersCP_tree_->Branch("sigmaPCA3", &stsCP_trackster_sigmaPCA3);
    simtrackstersCP_tree_->Branch("id_probabilities", &stsCP_trackster_id_probabilities);
    simtrackstersCP_tree_->Branch("vertices_indexes", &stsCP_trackster_vertices_indexes);
    simtrackstersCP_tree_->Branch("vertices_x", &stsCP_trackster_vertices_x);
    simtrackstersCP_tree_->Branch("vertices_y", &stsCP_trackster_vertices_y);
    simtrackstersCP_tree_->Branch("vertices_z", &stsCP_trackster_vertices_z);
    simtrackstersCP_tree_->Branch("vertices_time", &stsCP_trackster_vertices_time);
    simtrackstersCP_tree_->Branch("vertices_timeErr", &stsCP_trackster_vertices_timeErr);
    simtrackstersCP_tree_->Branch("vertices_energy", &stsCP_trackster_vertices_energy);
    simtrackstersCP_tree_->Branch("vertices_correctedEnergy", &stsCP_trackster_vertices_correctedEnergy);
    simtrackstersCP_tree_->Branch("vertices_correctedEnergyUncertainty",
                                  &stsCP_trackster_vertices_correctedEnergyUncertainty);
    simtrackstersCP_tree_->Branch("vertices_multiplicity", &stsCP_trackster_vertices_multiplicity);
  }

  if (saveTracks_) {
    tracks_tree_ = fs->make<TTree>("tracks", "Tracks");
    tracks_tree_->Branch("event", &ev_event_);
    tracks_tree_->Branch("track_id", &track_id);
    tracks_tree_->Branch("track_hgcal_x", &track_hgcal_x);
    tracks_tree_->Branch("track_hgcal_y", &track_hgcal_y);
    tracks_tree_->Branch("track_hgcal_z", &track_hgcal_z);
    tracks_tree_->Branch("track_hgcal_eta", &track_hgcal_eta);
    tracks_tree_->Branch("track_hgcal_phi", &track_hgcal_phi);
    tracks_tree_->Branch("track_hgcal_pt", &track_hgcal_pt);
    tracks_tree_->Branch("track_pt", &track_pt);
    tracks_tree_->Branch("track_missing_outer_hits", &track_missing_outer_hits);
    tracks_tree_->Branch("track_quality", &track_quality);
    tracks_tree_->Branch("track_charge", &track_charge);
    tracks_tree_->Branch("track_time", &track_time);
    tracks_tree_->Branch("track_time_quality", &track_time_quality);
    tracks_tree_->Branch("track_time_err", &track_time_err);
    tracks_tree_->Branch("track_beta", &track_beta);
    tracks_tree_->Branch("track_time_mtd", &track_time_mtd);
    tracks_tree_->Branch("track_time_mtd_err", &track_time_mtd_err);
    tracks_tree_->Branch("track_pos_mtd", &track_pos_mtd);
    tracks_tree_->Branch("track_nhits", &track_nhits);
  }

  if (saveSimTICLCandidate_) {
    simTICLCandidate_tree = fs->make<TTree>("simTICLCandidate", "Sim TICL Candidate");
    simTICLCandidate_tree->Branch("simTICLCandidate_raw_energy", &simTICLCandidate_raw_energy);
    simTICLCandidate_tree->Branch("simTICLCandidate_regressed_energy", &simTICLCandidate_regressed_energy);
    simTICLCandidate_tree->Branch("simTICLCandidate_simTracksterCPIndex", &simTICLCandidate_simTracksterCPIndex);
    simTICLCandidate_tree->Branch("simTICLCandidate_boundaryX", &simTICLCandidate_boundaryX);
    simTICLCandidate_tree->Branch("simTICLCandidate_boundaryY", &simTICLCandidate_boundaryY);
    simTICLCandidate_tree->Branch("simTICLCandidate_boundaryZ", &simTICLCandidate_boundaryZ);
    simTICLCandidate_tree->Branch("simTICLCandidate_boundaryPx", &simTICLCandidate_boundaryPx);
    simTICLCandidate_tree->Branch("simTICLCandidate_boundaryPy", &simTICLCandidate_boundaryPy);
    simTICLCandidate_tree->Branch("simTICLCandidate_boundaryPz", &simTICLCandidate_boundaryPz);
    simTICLCandidate_tree->Branch("simTICLCandidate_time", &simTICLCandidate_time);
    simTICLCandidate_tree->Branch("simTICLCandidate_caloParticleMass", &simTICLCandidate_caloParticleMass);
    simTICLCandidate_tree->Branch("simTICLCandidate_pdgId", &simTICLCandidate_pdgId);
    simTICLCandidate_tree->Branch("simTICLCandidate_charge", &simTICLCandidate_charge);
    simTICLCandidate_tree->Branch("simTICLCandidate_track_in_candidate", &simTICLCandidate_track_in_candidate);
  }
}

void TICLDumper::buildLayers() {
  // build disks at HGCal front & EM-Had interface for track propagation

  float zVal = hgcons_->waferZ(1, true);
  std::pair<float, float> rMinMax = hgcons_->rangeR(zVal, true);

  float zVal_interface = rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z();
  std::pair<float, float> rMinMax_interface = hgcons_->rangeR(zVal_interface, true);

  for (int iSide = 0; iSide < 2; ++iSide) {
    float zSide = (iSide == 0) ? (-1. * zVal) : zVal;
    firstDisk_[iSide] =
        std::make_unique<GeomDet>(Disk::build(Disk::PositionType(0, 0, zSide),
                                              Disk::RotationType(),
                                              SimpleDiskBounds(rMinMax.first, rMinMax.second, zSide - 0.5, zSide + 0.5))
                                      .get());

    zSide = (iSide == 0) ? (-1. * zVal_interface) : zVal_interface;
    interfaceDisk_[iSide] = std::make_unique<GeomDet>(
        Disk::build(Disk::PositionType(0, 0, zSide),
                    Disk::RotationType(),
                    SimpleDiskBounds(rMinMax_interface.first, rMinMax_interface.second, zSide - 0.5, zSide + 0.5))
            .get());
  }
}

void TICLDumper::initialize(const HGCalDDDConstants* hgcons,
                            const hgcal::RecHitTools rhtools,
                            const edm::ESHandle<MagneticField> bfieldH,
                            const edm::ESHandle<Propagator> propH) {
  hgcons_ = hgcons;
  rhtools_ = rhtools;
  buildLayers();

  bfield_ = bfieldH;
  propagator_ = propH;
}

void TICLDumper::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  ev_event_ += 1;
  clearVariables();
  auto bFieldProd = bfield_.product();
  const Propagator& prop = (*propagator_);

  //get all the layer clusters
  edm::Handle<std::vector<reco::CaloCluster>> layer_clusters_h;
  event.getByToken(layer_clusters_token_, layer_clusters_h);
  const auto& clusters = *layer_clusters_h;

  edm::Handle<edm::ValueMap<std::pair<float, float>>> clustersTime_h;
  event.getByToken(clustersTime_token_, clustersTime_h);
  const auto& layerClustersTimes = *clustersTime_h;

  //TICL Candidate
  edm::Handle<std::vector<TICLCandidate>> candidates_h;
  event.getByToken(ticl_candidates_token_, candidates_h);
  const auto& ticlcandidates = *candidates_h;
  edm::Handle<std::vector<ticl::Trackster>> ticlcandidates_tracksters_h = event.getHandle(ticl_candidates_tracksters_token_);

  //Track
  edm::Handle<std::vector<reco::Track>> tracks_h;
  event.getByToken(tracks_token_, tracks_h);
  const auto& tracks = *tracks_h;

  edm::Handle<edm::ValueMap<float>> trackTime_h;
  event.getByToken(tracks_time_token_, trackTime_h);
  const auto& trackTime = *trackTime_h;

  edm::Handle<edm::ValueMap<float>> trackTimeErr_h;
  event.getByToken(tracks_time_err_token_, trackTimeErr_h);
  const auto& trackTimeErr = *trackTimeErr_h;

  edm::Handle<edm::ValueMap<float>> trackBeta_h;
  event.getByToken(tracks_beta_token_, trackBeta_h);
  const auto& trackBeta = *trackBeta_h;

  edm::Handle<edm::ValueMap<float>> trackTimeQual_h;
  event.getByToken(tracks_time_quality_token_, trackTimeQual_h);
  const auto& trackTimeQual = *trackTimeQual_h;

  edm::Handle<edm::ValueMap<float>> trackTimeMtd_h;
  event.getByToken(tracks_time_mtd_token_, trackTimeMtd_h);
  const auto& trackTimeMtd = *trackTimeMtd_h;

  edm::Handle<edm::ValueMap<float>> trackTimeMtdErr_h;
  event.getByToken(tracks_time_mtd_err_token_, trackTimeMtdErr_h);
  const auto& trackTimeMtdErr = *trackTimeMtdErr_h;

  edm::Handle<edm::ValueMap<GlobalPoint>> trackPosMtd_h;
  event.getByToken(tracks_pos_mtd_token_, trackPosMtd_h);
  const auto& trackPosMtd = *trackPosMtd_h;

  // superclustering
  superclustering_linkedResultTracksters = event.get(superclustering_linkedResultTracksters_token);
  /* edm::Handle<ticl::SuperclusteringDNNScore> superclustering_DNNScore_h;
  event.getByToken(superclustering_DNNScore_token_, superclustering_DNNScore_h);
  // superclustering DNN score is not always saved
  if (superclustering_DNNScore_h.isValid())
    superclusteringDNNScore = *superclustering_DNNScore_h; 
  */

  // simTracksters from SC
  edm::Handle<std::vector<ticl::Trackster>> simTrackstersSC_h;
  event.getByToken(simTracksters_SC_token_, simTrackstersSC_h);
  const auto& simTrackstersSC = *simTrackstersSC_h;

  // simTracksters from CP
  edm::Handle<std::vector<ticl::Trackster>> simTrackstersCP_h;
  event.getByToken(simTracksters_CP_token_, simTrackstersCP_h);
  const auto& simTrackstersCP = *simTrackstersCP_h;

  edm::Handle<std::vector<TICLCandidate>> simTICLCandidates_h;
  event.getByToken(simTICLCandidate_token_, simTICLCandidates_h);
  const auto& simTICLCandidates = *simTICLCandidates_h;

  edm::Handle<std::vector<CaloParticle>> caloparticles_h;
  event.getByToken(caloparticles_token_, caloparticles_h);
  const auto& caloparticles = *caloparticles_h;

  const auto& simclusters = event.get(simclusters_token_);

  nclusters_ = clusters.size();

  // Save all the trackster collections
  for (unsigned int i = 0; i < tracksters_dumperHelpers_.size(); i++) {
    edm::Handle<std::vector<ticl::Trackster>> tracksters_handle;
    std::vector<ticl::Trackster> const& tracksters = event.get<std::vector<ticl::Trackster>>(tracksters_token_[i]);
    tracksters_dumperHelpers_[i].fillFromEvent(tracksters, clusters, layerClustersTimes);
    tracksters_trees[i]->Fill();
    tracksters_collection[i] = tracksters;
  }

  stsSC_ntracksters_ = simTrackstersSC.size();
  using CaloObjectVariant = std::variant<CaloParticle, SimCluster>;
  for (auto trackster_iterator = simTrackstersSC.begin(); trackster_iterator != simTrackstersSC.end();
       ++trackster_iterator) {
    //per-trackster analysis
    stsSC_trackster_time.push_back(trackster_iterator->time());
    stsSC_trackster_timeBoundary.push_back(trackster_iterator->boundaryTime());
    stsSC_trackster_timeError.push_back(trackster_iterator->timeError());
    stsSC_trackster_regressed_energy.push_back(trackster_iterator->regressed_energy());
    stsSC_trackster_raw_energy.push_back(trackster_iterator->raw_energy());
    stsSC_trackster_raw_em_energy.push_back(trackster_iterator->raw_em_energy());
    stsSC_trackster_raw_pt.push_back(trackster_iterator->raw_pt());
    stsSC_trackster_raw_em_pt.push_back(trackster_iterator->raw_em_pt());
    stsSC_trackster_barycenter_x.push_back(trackster_iterator->barycenter().x());
    stsSC_trackster_barycenter_y.push_back(trackster_iterator->barycenter().y());
    stsSC_trackster_barycenter_z.push_back(trackster_iterator->barycenter().z());
    stsSC_trackster_barycenter_eta.push_back(trackster_iterator->barycenter().eta());
    stsSC_trackster_barycenter_phi.push_back(trackster_iterator->barycenter().phi());
    stsSC_trackster_EV1.push_back(trackster_iterator->eigenvalues()[0]);
    stsSC_trackster_EV2.push_back(trackster_iterator->eigenvalues()[1]);
    stsSC_trackster_EV3.push_back(trackster_iterator->eigenvalues()[2]);
    stsSC_trackster_eVector0_x.push_back((trackster_iterator->eigenvectors()[0]).x());
    stsSC_trackster_eVector0_y.push_back((trackster_iterator->eigenvectors()[0]).y());
    stsSC_trackster_eVector0_z.push_back((trackster_iterator->eigenvectors()[0]).z());
    stsSC_trackster_sigmaPCA1.push_back(trackster_iterator->sigmasPCA()[0]);
    stsSC_trackster_sigmaPCA2.push_back(trackster_iterator->sigmasPCA()[1]);
    stsSC_trackster_sigmaPCA3.push_back(trackster_iterator->sigmasPCA()[2]);
    stsSC_pdgID.push_back(simclusters[trackster_iterator->seedIndex()].pdgId());

    CaloObjectVariant caloObj;
    if (trackster_iterator->seedID() == caloparticles_h.id()) {
      caloObj = caloparticles[trackster_iterator->seedIndex()];
    } else {
      caloObj = simclusters[trackster_iterator->seedIndex()];
    }

    auto const& simTrack = std::visit([](auto&& obj) { return obj.g4Tracks()[0]; }, caloObj);
    auto const& caloPt = std::visit([](auto&& obj) { return obj.pt(); }, caloObj);
    stsSC_trackster_regressed_pt.push_back(caloPt);
    if (simTrack.crossedBoundary()) {
      stsSC_boundaryX.push_back(simTrack.getPositionAtBoundary().x());
      stsSC_boundaryY.push_back(simTrack.getPositionAtBoundary().y());
      stsSC_boundaryZ.push_back(simTrack.getPositionAtBoundary().z());
      stsSC_boundaryEta.push_back(simTrack.getPositionAtBoundary().eta());
      stsSC_boundaryPhi.push_back(simTrack.getPositionAtBoundary().phi());
      stsSC_boundaryPx.push_back(simTrack.getMomentumAtBoundary().x());
      stsSC_boundaryPy.push_back(simTrack.getMomentumAtBoundary().y());
      stsSC_boundaryPz.push_back(simTrack.getMomentumAtBoundary().z());
    } else {
      stsSC_boundaryX.push_back(-999);
      stsSC_boundaryY.push_back(-999);
      stsSC_boundaryZ.push_back(-999);
      stsSC_boundaryEta.push_back(-999);
      stsSC_boundaryPhi.push_back(-999);
      stsSC_boundaryPx.push_back(-999);
      stsSC_boundaryPy.push_back(-999);
      stsSC_boundaryPz.push_back(-999);
    }
    auto const trackIdx = trackster_iterator->trackIdx();
    stsSC_trackIdx.push_back(trackIdx);
    if (trackIdx != -1) {
      const auto& track = tracks[trackIdx];

      int iSide = int(track.eta() > 0);

      const auto& fts = trajectoryStateTransform::outerFreeState((track), bFieldProd);
      // to the HGCal front
      const auto& tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
      if (tsos.isValid()) {
        const auto& globalPos = tsos.globalPosition();
        const auto& globalMom = tsos.globalMomentum();
        stsSC_track_boundaryX.push_back(globalPos.x());
        stsSC_track_boundaryY.push_back(globalPos.y());
        stsSC_track_boundaryZ.push_back(globalPos.z());
        stsSC_track_boundaryEta.push_back(globalPos.eta());
        stsSC_track_boundaryPhi.push_back(globalPos.phi());
        stsSC_track_boundaryPx.push_back(globalMom.x());
        stsSC_track_boundaryPy.push_back(globalMom.y());
        stsSC_track_boundaryPz.push_back(globalMom.z());
        stsSC_trackTime.push_back(track.t0());
      } else {
        stsSC_track_boundaryX.push_back(-999);
        stsSC_track_boundaryY.push_back(-999);
        stsSC_track_boundaryZ.push_back(-999);
        stsSC_track_boundaryEta.push_back(-999);
        stsSC_track_boundaryPhi.push_back(-999);
        stsSC_track_boundaryPx.push_back(-999);
        stsSC_track_boundaryPy.push_back(-999);
        stsSC_track_boundaryPz.push_back(-999);
      }
    } else {
      stsSC_track_boundaryX.push_back(-999);
      stsSC_track_boundaryY.push_back(-999);
      stsSC_track_boundaryZ.push_back(-999);
      stsSC_track_boundaryEta.push_back(-999);
      stsSC_track_boundaryPhi.push_back(-999);
      stsSC_track_boundaryPx.push_back(-999);
      stsSC_track_boundaryPy.push_back(-999);
      stsSC_track_boundaryPz.push_back(-999);
    }

    std::vector<float> id_probs;
    for (size_t i = 0; i < 8; i++)
      id_probs.push_back(trackster_iterator->id_probabilities(i));
    stsSC_trackster_id_probabilities.push_back(id_probs);

    // Clusters
    std::vector<uint32_t> vertices_indexes;
    std::vector<float> vertices_x;
    std::vector<float> vertices_y;
    std::vector<float> vertices_z;
    std::vector<float> vertices_time;
    std::vector<float> vertices_timeErr;
    std::vector<float> vertices_energy;
    std::vector<float> vertices_correctedEnergy;
    std::vector<float> vertices_correctedEnergyUncertainty;
    for (auto idx : trackster_iterator->vertices()) {
      vertices_indexes.push_back(idx);
      auto associated_cluster = (*layer_clusters_h)[idx];
      vertices_x.push_back(associated_cluster.x());
      vertices_y.push_back(associated_cluster.y());
      vertices_z.push_back(associated_cluster.z());
      vertices_energy.push_back(associated_cluster.energy());
      vertices_correctedEnergy.push_back(associated_cluster.correctedEnergy());
      vertices_correctedEnergyUncertainty.push_back(associated_cluster.correctedEnergyUncertainty());
      vertices_time.push_back(layerClustersTimes.get(idx).first);
      vertices_timeErr.push_back(layerClustersTimes.get(idx).second);
    }
    stsSC_trackster_vertices_indexes.push_back(vertices_indexes);
    stsSC_trackster_vertices_x.push_back(vertices_x);
    stsSC_trackster_vertices_y.push_back(vertices_y);
    stsSC_trackster_vertices_z.push_back(vertices_z);
    stsSC_trackster_vertices_time.push_back(vertices_time);
    stsSC_trackster_vertices_timeErr.push_back(vertices_timeErr);
    stsSC_trackster_vertices_energy.push_back(vertices_energy);
    stsSC_trackster_vertices_correctedEnergy.push_back(vertices_correctedEnergy);
    stsSC_trackster_vertices_correctedEnergyUncertainty.push_back(vertices_correctedEnergyUncertainty);

    // Multiplicity
    std::vector<float> vertices_multiplicity;
    for (auto multiplicity : trackster_iterator->vertex_multiplicity()) {
      vertices_multiplicity.push_back(multiplicity);
    }
    stsSC_trackster_vertices_multiplicity.push_back(vertices_multiplicity);
  }

  stsCP_ntracksters_ = simTrackstersCP.size();

  for (auto trackster_iterator = simTrackstersCP.begin(); trackster_iterator != simTrackstersCP.end();
       ++trackster_iterator) {
    //per-trackster analysis
    stsCP_trackster_time.push_back(trackster_iterator->time());
    stsCP_trackster_timeBoundary.push_back(trackster_iterator->boundaryTime());
    stsCP_trackster_timeError.push_back(trackster_iterator->timeError());
    stsCP_trackster_regressed_energy.push_back(trackster_iterator->regressed_energy());
    stsCP_trackster_raw_energy.push_back(trackster_iterator->raw_energy());
    stsCP_trackster_raw_em_energy.push_back(trackster_iterator->raw_em_energy());
    stsCP_trackster_raw_pt.push_back(trackster_iterator->raw_pt());
    stsCP_trackster_raw_em_pt.push_back(trackster_iterator->raw_em_pt());
    stsCP_trackster_barycenter_x.push_back(trackster_iterator->barycenter().x());
    stsCP_trackster_barycenter_y.push_back(trackster_iterator->barycenter().y());
    stsCP_trackster_barycenter_z.push_back(trackster_iterator->barycenter().z());
    stsCP_trackster_barycenter_eta.push_back(trackster_iterator->barycenter().eta());
    stsCP_trackster_barycenter_phi.push_back(trackster_iterator->barycenter().phi());
    stsCP_trackster_EV1.push_back(trackster_iterator->eigenvalues()[0]);
    stsCP_trackster_EV2.push_back(trackster_iterator->eigenvalues()[1]);
    stsCP_trackster_EV3.push_back(trackster_iterator->eigenvalues()[2]);
    stsCP_trackster_eVector0_x.push_back((trackster_iterator->eigenvectors()[0]).x());
    stsCP_trackster_eVector0_y.push_back((trackster_iterator->eigenvectors()[0]).y());
    stsCP_trackster_eVector0_z.push_back((trackster_iterator->eigenvectors()[0]).z());
    stsCP_trackster_sigmaPCA1.push_back(trackster_iterator->sigmasPCA()[0]);
    stsCP_trackster_sigmaPCA2.push_back(trackster_iterator->sigmasPCA()[1]);
    stsCP_trackster_sigmaPCA3.push_back(trackster_iterator->sigmasPCA()[2]);
    stsCP_pdgID.push_back(caloparticles[trackster_iterator->seedIndex()].pdgId());
    CaloObjectVariant caloObj;
    if (trackster_iterator->seedID() == caloparticles_h.id()) {
      caloObj = caloparticles[trackster_iterator->seedIndex()];
    } else {
      caloObj = simclusters[trackster_iterator->seedIndex()];
    }

    auto const& simTrack = std::visit([](auto&& obj) { return obj.g4Tracks()[0]; }, caloObj);
    auto const& caloPt = std::visit([](auto&& obj) { return obj.pt(); }, caloObj);
    stsCP_trackster_regressed_pt.push_back(caloPt);

    if (simTrack.crossedBoundary()) {
      stsCP_boundaryX.push_back(simTrack.getPositionAtBoundary().x());
      stsCP_boundaryY.push_back(simTrack.getPositionAtBoundary().y());
      stsCP_boundaryZ.push_back(simTrack.getPositionAtBoundary().z());
      stsCP_boundaryEta.push_back(simTrack.getPositionAtBoundary().eta());
      stsCP_boundaryPhi.push_back(simTrack.getPositionAtBoundary().phi());
      stsCP_boundaryPx.push_back(simTrack.getMomentumAtBoundary().x());
      stsCP_boundaryPy.push_back(simTrack.getMomentumAtBoundary().y());
      stsCP_boundaryPz.push_back(simTrack.getMomentumAtBoundary().z());
    } else {
      stsCP_boundaryX.push_back(-999);
      stsCP_boundaryY.push_back(-999);
      stsCP_boundaryZ.push_back(-999);
      stsCP_boundaryEta.push_back(-999);
      stsCP_boundaryPhi.push_back(-999);
      stsCP_boundaryPx.push_back(-999);
      stsCP_boundaryPy.push_back(-999);
      stsCP_boundaryPz.push_back(-999);
    }
    auto const trackIdx = trackster_iterator->trackIdx();
    stsCP_trackIdx.push_back(trackIdx);
    if (trackIdx != -1) {
      const auto& track = tracks[trackIdx];

      int iSide = int(track.eta() > 0);

      const auto& fts = trajectoryStateTransform::outerFreeState((track), bFieldProd);
      // to the HGCal front
      const auto& tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
      if (tsos.isValid()) {
        const auto& globalPos = tsos.globalPosition();
        const auto& globalMom = tsos.globalMomentum();
        stsCP_track_boundaryX.push_back(globalPos.x());
        stsCP_track_boundaryY.push_back(globalPos.y());
        stsCP_track_boundaryZ.push_back(globalPos.z());
        stsCP_track_boundaryEta.push_back(globalPos.eta());
        stsCP_track_boundaryPhi.push_back(globalPos.phi());
        stsCP_track_boundaryPx.push_back(globalMom.x());
        stsCP_track_boundaryPy.push_back(globalMom.y());
        stsCP_track_boundaryPz.push_back(globalMom.z());
        stsCP_trackTime.push_back(track.t0());
      } else {
        stsCP_track_boundaryX.push_back(-999);
        stsCP_track_boundaryY.push_back(-999);
        stsCP_track_boundaryZ.push_back(-999);
        stsCP_track_boundaryEta.push_back(-999);
        stsCP_track_boundaryPhi.push_back(-999);
        stsCP_track_boundaryPx.push_back(-999);
        stsCP_track_boundaryPy.push_back(-999);
        stsCP_track_boundaryPz.push_back(-999);
      }
    } else {
      stsCP_track_boundaryX.push_back(-999);
      stsCP_track_boundaryY.push_back(-999);
      stsCP_track_boundaryZ.push_back(-999);
      stsCP_track_boundaryEta.push_back(-999);
      stsCP_track_boundaryPhi.push_back(-999);
      stsCP_track_boundaryPx.push_back(-999);
      stsCP_track_boundaryPy.push_back(-999);
      stsCP_track_boundaryPz.push_back(-999);
    }
    std::vector<float> id_probs;
    for (size_t i = 0; i < 8; i++)
      id_probs.push_back(trackster_iterator->id_probabilities(i));
    stsCP_trackster_id_probabilities.push_back(id_probs);

    // Clusters
    std::vector<uint32_t> vertices_indexes;
    std::vector<float> vertices_x;
    std::vector<float> vertices_y;
    std::vector<float> vertices_z;
    std::vector<float> vertices_time;
    std::vector<float> vertices_timeErr;
    std::vector<float> vertices_energy;
    std::vector<float> vertices_correctedEnergy;
    std::vector<float> vertices_correctedEnergyUncertainty;
    for (auto idx : trackster_iterator->vertices()) {
      vertices_indexes.push_back(idx);
      auto associated_cluster = (*layer_clusters_h)[idx];
      vertices_x.push_back(associated_cluster.x());
      vertices_y.push_back(associated_cluster.y());
      vertices_z.push_back(associated_cluster.z());
      vertices_energy.push_back(associated_cluster.energy());
      vertices_correctedEnergy.push_back(associated_cluster.correctedEnergy());
      vertices_correctedEnergyUncertainty.push_back(associated_cluster.correctedEnergyUncertainty());
      vertices_time.push_back(layerClustersTimes.get(idx).first);
      vertices_timeErr.push_back(layerClustersTimes.get(idx).second);
    }
    stsCP_trackster_vertices_indexes.push_back(vertices_indexes);
    stsCP_trackster_vertices_x.push_back(vertices_x);
    stsCP_trackster_vertices_y.push_back(vertices_y);
    stsCP_trackster_vertices_z.push_back(vertices_z);
    stsCP_trackster_vertices_time.push_back(vertices_time);
    stsCP_trackster_vertices_timeErr.push_back(vertices_timeErr);
    stsCP_trackster_vertices_energy.push_back(vertices_energy);
    stsCP_trackster_vertices_correctedEnergy.push_back(vertices_correctedEnergy);
    stsCP_trackster_vertices_correctedEnergyUncertainty.push_back(vertices_correctedEnergyUncertainty);

    // Multiplicity
    std::vector<float> vertices_multiplicity;
    for (auto multiplicity : trackster_iterator->vertex_multiplicity()) {
      vertices_multiplicity.push_back(multiplicity);
    }
    stsCP_trackster_vertices_multiplicity.push_back(vertices_multiplicity);
  }

  simTICLCandidate_track_in_candidate.resize(simTICLCandidates.size(), -1);
  for (size_t i = 0; i < simTICLCandidates.size(); ++i) {
    auto const& cand = simTICLCandidates[i];

    simTICLCandidate_raw_energy.push_back(cand.rawEnergy());
    simTICLCandidate_regressed_energy.push_back(cand.p4().energy());
    simTICLCandidate_pdgId.push_back(cand.pdgId());
    simTICLCandidate_charge.push_back(cand.charge());
    simTICLCandidate_time.push_back(cand.time());
    std::vector<int> tmpIdxVec;
    for (auto const& simTS : cand.tracksters()) {
      auto trackster_idx = simTS.get() - (edm::Ptr<ticl::Trackster>(simTrackstersSC_h, 0)).get();
      tmpIdxVec.push_back(trackster_idx);
    }
    simTICLCandidate_simTracksterCPIndex.push_back(tmpIdxVec);
    tmpIdxVec.clear();
    auto const& trackPtr = cand.trackPtr();
    if (!trackPtr.isNull()) {
      auto const& track = *trackPtr;
      int iSide = int(track.eta() > 0);
      int tk_idx = trackPtr.get() - (edm::Ptr<reco::Track>(tracks_h, 0)).get();
      simTICLCandidate_track_in_candidate[i] = tk_idx;

      const auto& fts = trajectoryStateTransform::outerFreeState((track), bFieldProd);
      // to the HGCal front
      const auto& tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
      if (tsos.isValid()) {
        const auto& globalPos = tsos.globalPosition();
        const auto& globalMom = tsos.globalMomentum();
        simTICLCandidate_boundaryX.push_back(globalPos.x());
        simTICLCandidate_boundaryY.push_back(globalPos.y());
        simTICLCandidate_boundaryZ.push_back(globalPos.z());
        simTICLCandidate_boundaryPx.push_back(globalMom.x());
        simTICLCandidate_boundaryPy.push_back(globalMom.y());
        simTICLCandidate_boundaryPz.push_back(globalMom.z());
      } else {
        simTICLCandidate_boundaryX.push_back(-999);
        simTICLCandidate_boundaryY.push_back(-999);
        simTICLCandidate_boundaryZ.push_back(-999);
        simTICLCandidate_boundaryPx.push_back(-999);
        simTICLCandidate_boundaryPy.push_back(-999);
        simTICLCandidate_boundaryPz.push_back(-999);
      }
    } else {
      simTICLCandidate_boundaryX.push_back(-999);
      simTICLCandidate_boundaryY.push_back(-999);
      simTICLCandidate_boundaryZ.push_back(-999);
      simTICLCandidate_boundaryPx.push_back(-999);
      simTICLCandidate_boundaryPy.push_back(-999);
      simTICLCandidate_boundaryPz.push_back(-999);
    }
  }

  int c_id = 0;

  for (auto cluster_iterator = clusters.begin(); cluster_iterator != clusters.end(); ++cluster_iterator) {
    auto lc_seed = cluster_iterator->seed();
    cluster_seedID.push_back(lc_seed);
    cluster_energy.push_back(cluster_iterator->energy());
    cluster_correctedEnergy.push_back(cluster_iterator->correctedEnergy());
    cluster_correctedEnergyUncertainty.push_back(cluster_iterator->correctedEnergyUncertainty());
    cluster_position_x.push_back(cluster_iterator->x());
    cluster_position_y.push_back(cluster_iterator->y());
    cluster_position_z.push_back(cluster_iterator->z());
    cluster_position_eta.push_back(cluster_iterator->eta());
    cluster_position_phi.push_back(cluster_iterator->phi());
    auto haf = cluster_iterator->hitsAndFractions();
    auto layerId = rhtools_.getLayerWithOffset(haf[0].first);
    cluster_layer_id.push_back(layerId);
    uint32_t number_of_hits = cluster_iterator->hitsAndFractions().size();
    cluster_number_of_hits.push_back(number_of_hits);
    cluster_type.push_back(rhtools_.getCellType(lc_seed));
    cluster_timeErr.push_back(layerClustersTimes.get(c_id).second);
    cluster_time.push_back(layerClustersTimes.get(c_id).first);
    c_id += 1;
  }

  tracksters_in_candidate.resize(ticlcandidates.size());
  track_in_candidate.resize(ticlcandidates.size(), -1);
  nCandidates = ticlcandidates.size();
  for (int i = 0; i < static_cast<int>(ticlcandidates.size()); ++i) {
    const auto& candidate = ticlcandidates[i];
    candidate_charge.push_back(candidate.charge());
    candidate_pdgId.push_back(candidate.pdgId());
    candidate_energy.push_back(candidate.energy());
    candidate_px.push_back(candidate.px());
    candidate_py.push_back(candidate.py());
    candidate_pz.push_back(candidate.pz());
    candidate_time.push_back(candidate.time());
    candidate_time_err.push_back(candidate.timeError());
    std::vector<float> id_probs;
    for (int j = 0; j < 8; j++) {
      ticl::Trackster::ParticleType type = static_cast<ticl::Trackster::ParticleType>(j);
      id_probs.push_back(candidate.id_probability(type));
    }
    candidate_id_probabilities.push_back(id_probs);

    auto trackster_ptrs = candidate.tracksters();
    auto track_ptr = candidate.trackPtr();
    for (const auto& ts_ptr : trackster_ptrs) {
      auto ts_idx = ts_ptr.get() - (edm::Ptr<ticl::Trackster>(ticlcandidates_tracksters_h, 0)).get();
      tracksters_in_candidate[i].push_back(ts_idx);
    }
    if (track_ptr.isNull())
      continue;
    int tk_idx = track_ptr.get() - (edm::Ptr<reco::Track>(tracks_h, 0)).get();
    track_in_candidate[i] = tk_idx;
  }


  // trackster to simTrackster associations
  for (unsigned int i = 0; i < associations_dumperHelpers_.size(); i++) {
    associations_dumperHelpers_[i].fillFromEvent(event.getHandle(associations_tracksterCollection_[i]), event.getHandle(associations_simTracksterCollection_[i]), 
      event.get(associations_recoToSim_token_[i]), event.get(associations_simToReco_token_[i]));
  }
  if (associations_dumperHelpers_.size() > 0)
    associations_tree_->Fill();


  //Tracks
  for (size_t i = 0; i < tracks.size(); i++) {
    const auto& track = tracks[i];
    reco::TrackRef trackref = reco::TrackRef(tracks_h, i);
    int iSide = int(track.eta() > 0);
    const auto& fts = trajectoryStateTransform::outerFreeState((track), bFieldProd);
    // to the HGCal front
    const auto& tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
    if (tsos.isValid()) {
      const auto& globalPos = tsos.globalPosition();
      const auto& globalMom = tsos.globalMomentum();
      track_id.push_back(i);
      track_hgcal_x.push_back(globalPos.x());
      track_hgcal_y.push_back(globalPos.y());
      track_hgcal_z.push_back(globalPos.z());
      track_hgcal_eta.push_back(globalPos.eta());
      track_hgcal_phi.push_back(globalPos.phi());
      track_hgcal_px.push_back(globalMom.x());
      track_hgcal_py.push_back(globalMom.y());
      track_hgcal_pz.push_back(globalMom.z());
      track_hgcal_pt.push_back(globalMom.perp());
      track_pt.push_back(track.pt());
      track_quality.push_back(track.quality(reco::TrackBase::highPurity));
      track_missing_outer_hits.push_back(track.missingOuterHits());
      track_charge.push_back(track.charge());
      track_time.push_back(trackTime[trackref]);
      track_time_quality.push_back(trackTimeQual[trackref]);
      track_time_err.push_back(trackTimeErr[trackref]);
      track_beta.push_back(trackBeta[trackref]);
      track_time_mtd.push_back(trackTimeMtd[trackref]);
      track_time_mtd_err.push_back(trackTimeMtdErr[trackref]);
      track_pos_mtd.push_back(trackPosMtd[trackref]);
      track_nhits.push_back(tracks[i].recHitsSize());
    }
  }

  if (saveLCs_)
    cluster_tree_->Fill();
  if (saveTICLCandidate_)
    candidate_tree_->Fill();
  if (saveSuperclustering_ /*|| saveSuperclusteringDNNScore_ */)
    superclustering_tree_->Fill();
  if (saveSimTrackstersSC_)
    simtrackstersSC_tree_->Fill();
  if (saveSimTrackstersCP_)
    simtrackstersCP_tree_->Fill();
  if (saveTracks_)
    tracks_tree_->Fill();
  if (saveSimTICLCandidate_)
    simTICLCandidate_tree->Fill();
}

void TICLDumper::endJob() {}

void TICLDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  edm::ParameterSetDescription tracksterDescValidator;
  tracksterDescValidator.add<std::string>("treeName")->setComment("Name of the output tree for the trackster collection");
  tracksterDescValidator.add<edm::InputTag>("tracksterInputTag")->setComment("Input tag for the trackster collection to write");
  desc.addVPSet("tracksterCollections", tracksterDescValidator)->setComment("Trackster collections to dump");

  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalMergeLayerClusters"));
  desc.add<edm::InputTag>("layer_clustersTime", edm::InputTag("hgcalMergeLayerClusters", "timeLayerCluster"));
  desc.add<edm::InputTag>("ticlcandidates", edm::InputTag("ticlTrackstersMerge"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("tracksTime", edm::InputTag("tofPID:t0"));
  desc.add<edm::InputTag>("tracksTimeQual", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  desc.add<edm::InputTag>("tracksTimeErr", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("tracksBeta", edm::InputTag("trackExtenderWithMTD:generalTrackBeta"));
  desc.add<edm::InputTag>("tracksTimeMtd", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
  desc.add<edm::InputTag>("tracksTimeMtdErr", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"));
  desc.add<edm::InputTag>("tracksPosMtd", edm::InputTag("trackExtenderWithMTD:generalTrackmtdpos"));
  desc.add<edm::InputTag>("superclustering", edm::InputTag("ticlTracksterLinksSuperclustering"));
  //desc.add<edm::InputTag>("superclusteringDNNScore", edm::InputTag("ticlTrackstersSuperclustering", "superclusteringTracksterDNNScore"));
  desc.add<edm::InputTag>("simtrackstersSC", edm::InputTag("ticlSimTracksters"));
  desc.add<edm::InputTag>("simtrackstersCP", edm::InputTag("ticlSimTracksters", "fromCPs"));
  desc.add<edm::InputTag>("simTICLCandidates", edm::InputTag("ticlSimTracksters"));

  edm::ParameterSetDescription associatorDescValidator;
  associatorDescValidator.add<std::string>("branchName")->setComment("Name of the output branches in the tree");
  associatorDescValidator.add<std::string>("suffix")->setComment("Should be CP or SC (for the output branch name)");
  associatorDescValidator.add<edm::InputTag>("associatorInputTag")->setComment("Input tag for the associator (do not put the instance)");
  associatorDescValidator.add<edm::InputTag>("tracksterCollection")->setComment("Collection of tracksters used by the associator");
  associatorDescValidator.add<edm::InputTag>("simTracksterCollection")->setComment("Collection of SimTrackster used by the associator");
  desc.addVPSet("associators", associatorDescValidator)->setComment("Associators to dump");

  desc.add<edm::InputTag>("simclusters", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("caloparticles", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<std::string>("detector", "HGCAL");
  desc.add<std::string>("propagator", "PropagatorWithMaterial");

  desc.add<bool>("saveLCs", true);
  desc.add<bool>("saveSimTrackstersSC", true);
  desc.add<bool>("saveSimTrackstersCP", true);
  desc.add<bool>("saveTICLCandidate", true);
  desc.add<bool>("saveSimTICLCandidate", true);
  desc.add<bool>("saveTracks", true);
  desc.add<bool>("saveSuperclustering", true);
  //desc.add<bool>("saveSuperclusteringDNNScore", false)
  //  ->setComment("Save the superclustering DNN score for all the evaluations made. Takes a large amount of disk space. Make sure the DNN score is also saved into the event as well (otherwise it won't dump anything)");
  descriptions.add("ticlDumper", desc);
}

DEFINE_FWK_MODULE(TICLDumper);
