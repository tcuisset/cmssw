#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/ComputeClusterTime.h"
#include "TrackstersPCA.h"


#include <iostream>
#include <set>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <functional>



void ticl::assignPCAtoTracksters(std::vector<Trackster> &tracksters,
                                 const std::vector<reco::CaloCluster> &layerClusters,
                                 const edm::ValueMap<std::pair<float, float>> &layerClustersTime,
				 double z_limit_em,
				 bool energyWeight,
				 const hgcal::RecHitTools rhtools,
				 int minLayer,
				 int maxLayer,
				 bool clean) {

  LogDebug("TrackstersPCA_Eigen") << "------- Eigen -------" << std::endl;

  for (auto &trackster : tracksters) {

    std::cout<<"start testing teackster with size:"<<trackster.vertices().size()<<std::endl;

    Eigen::Vector3d point;
    point << 0., 0., 0.;
    Eigen::Vector3d barycenter;
    barycenter << 0., 0., 0.;
    Eigen::Vector3d filtered_barycenter;
    filtered_barycenter << 0., 0., 0.;

    auto fillPoint = [&](const reco::CaloCluster &c, const float weight = 1.f) {
      point[0] = weight * c.x();
      point[1] = weight * c.y();
      point[2] = weight * c.z();
    };

    // Initialize this trackster with default, dummy values
    trackster.setRawEnergy(0.f);
    trackster.setRawEmEnergy(0.f);
    trackster.setRawPt(0.f);
    trackster.setRawEmPt(0.f);

    size_t N = trackster.vertices().size();
    float weight = 1.f / N;
    float weights2_sum = 0.f;
    Eigen::Vector3d sigmas;
    sigmas << 0., 0., 0.;
    Eigen::Vector3d sigmasEigen;
    sigmasEigen << 0., 0., 0.;
    Eigen::Matrix3d covM = Eigen::Matrix3d::Zero();

    std::vector<float> times;
    std::vector<float> timeErrors;
    std::set<uint32_t> usedLC;
    std::vector<double> layerClusterEnergies;

    for (size_t i = 0; i < N; ++i) {
      auto fraction = 1.f / trackster.vertex_multiplicity(i);
      trackster.addToRawEnergy(layerClusters[trackster.vertices(i)].energy() * fraction);
      if (std::abs(layerClusters[trackster.vertices(i)].z()) <= z_limit_em)
        trackster.addToRawEmEnergy(layerClusters[trackster.vertices(i)].energy() * fraction);

      // Compute the weighted barycenter.
      if (energyWeight)
        weight = layerClusters[trackster.vertices(i)].energy() * fraction;
      fillPoint(layerClusters[trackster.vertices(i)], weight);
      for (size_t j = 0; j < 3; ++j)
        barycenter[j] += point[j];

      // Add timing from layerClusters not already used
      if ((usedLC.insert(trackster.vertices(i))).second) {
        float timeE = layerClustersTime.get(trackster.vertices(i)).second;
        if (timeE > 0.f) {
          times.push_back(layerClustersTime.get(trackster.vertices(i)).first);
          timeErrors.push_back(1. / pow(timeE, 2));
        }
      }
      layerClusterEnergies.push_back(layerClusters[trackster.vertices(i)].energy());
    }
    if (energyWeight && trackster.raw_energy())
      barycenter /= trackster.raw_energy();

    hgcalsimclustertime::ComputeClusterTime timeEstimator;
    std::pair<float, float> timeTrackster = timeEstimator.fixSizeHighestDensity(times, timeErrors);
    std::cout<<"cleaning is  :"<<clean<<std::endl;

    std::vector<unsigned> filtered_idx;
    double filtered_energy = 0;
    if (clean) {
      // Filter layerclusters for the cleaned PCA
      auto maxE_vertex = std::distance(layerClusterEnergies.begin(), std::max_element(layerClusterEnergies.begin(), layerClusterEnergies.end()));
      auto maxE_layer = getLayerFromLC(layerClusters[trackster.vertices(maxE_vertex)], rhtools);

      auto vertices_by_layer = sortByLayer(trackster, layerClusters, rhtools);
      
      for (unsigned i = 1; i <= rhtools.lastLayer(); ++i) {
	auto vertices_in_layer = vertices_by_layer[i];
	if (vertices_in_layer.empty())
	  continue;
	
	std::vector<double> energies_in_layer;
	for (auto vrt : vertices_in_layer)
	  energies_in_layer.push_back(layerClusters[trackster.vertices(vrt)].energy());
	
	unsigned maxEid_inLayer = std::distance(energies_in_layer.begin(), std::max_element(energies_in_layer.begin(), energies_in_layer.end()));
	
	// layer based filtering of what goes into the PCA
	if ((int)i >= (int)maxE_layer - minLayer && (int)i <= (int)maxE_layer + maxLayer) {
	  auto filtered_vert = vertices_in_layer[maxEid_inLayer];
	  filtered_idx.push_back(filtered_vert);
	  
	  const auto maxE_LC = layerClusters[trackster.vertices(filtered_vert)];
	  fillPoint(maxE_LC, maxE_LC.energy() * (1.f / trackster.vertex_multiplicity(filtered_vert)));
	  for (size_t j = 0; j < 3; ++j)
	    filtered_barycenter[j] += point[j];
	  filtered_energy += maxE_LC.energy();
	}
      }
      
      filtered_barycenter /= filtered_energy;
    }
    std::cout<<"min, max "<<minLayer<<"  "<<maxLayer<<std::endl;
    
    
    auto calc_covM = [&](size_t i) {
      fillPoint(layerClusters[trackster.vertices(i)]);
      if (energyWeight && trackster.raw_energy()) {
        weight = (layerClusters[trackster.vertices(i)].energy() / trackster.vertex_multiplicity(i)) / (clean ? filtered_energy : trackster.raw_energy());
	if (trackster.vertex_multiplicity(i) > 1 )
	  std::cout<<"trackster.vertex_multiplicity(i)   :"<<trackster.vertex_multiplicity(i)<<std::endl;
      }
      weights2_sum += weight * weight;
      for (size_t x = 0; x < 3; ++x)
        for (size_t y = 0; y <= x; ++y) {
	  covM(x, y) += weight * (point[x] - (clean ? filtered_barycenter[x] : barycenter[x])) * (point[y] - (clean ? filtered_barycenter[y] : barycenter[y]));
	  covM(y, x) = covM(x, y);
        }
    };


    if (clean) {
      for (size_t i : filtered_idx) {
        calc_covM(i);
      }
    } else {
      for (size_t i = 0; i < N; ++i) {
        calc_covM(i);
      }
    }
    
    


    covM *= 1. / (1. - weights2_sum);
    
    // Perform the actual decomposition
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>::RealVectorType eigenvalues_fromEigen;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>::EigenvectorsType eigenvectors_fromEigen;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(covM);
    if (eigensolver.info() != Eigen::Success) {
      eigenvalues_fromEigen = eigenvalues_fromEigen.Zero();
      eigenvectors_fromEigen = eigenvectors_fromEigen.Zero();
    } else {
      eigenvalues_fromEigen = eigensolver.eigenvalues();
      eigenvectors_fromEigen = eigensolver.eigenvectors();
    }


    auto calc_spread = [&](size_t i) {
      fillPoint(layerClusters[trackster.vertices(i)]);
      sigmas += weight * (point - (clean ? filtered_barycenter : barycenter)).cwiseAbs2();
      Eigen::Vector3d point_transformed = eigenvectors_fromEigen * (point - (clean ? filtered_barycenter : barycenter));
      if (energyWeight && trackster.raw_energy())
        weight = (layerClusters[trackster.vertices(i)].energy() / trackster.vertex_multiplicity(i)) / (clean ? filtered_energy : trackster.raw_energy());
      sigmasEigen += weight * (point_transformed.cwiseAbs2());
    };

    if (clean) {
      for (size_t i : filtered_idx) {
        calc_spread(i);
      }
    } else {
      for (size_t i = 0; i < N; ++i) {
        calc_spread(i);
      }
    }





    sigmas /= (1. - weights2_sum);
    sigmasEigen /= (1. - weights2_sum);

    // Add trackster attributes
    trackster.setBarycenter(ticl::Trackster::Vector(barycenter));
    trackster.setTimeAndError(timeTrackster.first, timeTrackster.second);
    trackster.fillPCAVariables(
        eigenvalues_fromEigen, eigenvectors_fromEigen, sigmas, sigmasEigen, 3, ticl::Trackster::PCAOrdering::ascending);

    std::cout<<"covM:"<<covM <<" "<<covM.norm()<< std::endl;
    std::cout<<"<-----------tested this trackster------->"<<std::endl;

    LogDebug("TrackstersPCA") << "Use energy weighting: " << energyWeight << std::endl;
    LogDebug("TrackstersPCA") << "\nTrackster characteristics: " << std::endl;
    LogDebug("TrackstersPCA") << "Size: " << N << std::endl;
    LogDebug("TrackstersPCA") << "Energy: " << trackster.raw_energy() << std::endl;
    LogDebug("TrackstersPCA") << "raw_pt: " << trackster.raw_pt() << std::endl;
    LogDebug("TrackstersPCA") << "Means:          " << barycenter[0] << ", " << barycenter[1] << ", " << barycenter[2]
                              << std::endl;
    LogDebug("TrackstersPCA") << "Time:          " << trackster.time() << " +/- " << trackster.timeError() << std::endl;
    LogDebug("TrackstersPCA") << "EigenValues from Eigen/Tr(cov): " << eigenvalues_fromEigen[2] / covM.trace() << ", "
                              << eigenvalues_fromEigen[1] / covM.trace() << ", "
                              << eigenvalues_fromEigen[0] / covM.trace() << std::endl;
    LogDebug("TrackstersPCA") << "EigenValues from Eigen:         " << eigenvalues_fromEigen[2] << ", "
                              << eigenvalues_fromEigen[1] << ", " << eigenvalues_fromEigen[0] << std::endl;
    LogDebug("TrackstersPCA") << "EigenVector 3 from Eigen: " << eigenvectors_fromEigen(0, 2) << ", "
                              << eigenvectors_fromEigen(1, 2) << ", " << eigenvectors_fromEigen(2, 2) << std::endl;
    LogDebug("TrackstersPCA") << "EigenVector 2 from Eigen: " << eigenvectors_fromEigen(0, 1) << ", "
                              << eigenvectors_fromEigen(1, 1) << ", " << eigenvectors_fromEigen(2, 1) << std::endl;
    LogDebug("TrackstersPCA") << "EigenVector 1 from Eigen: " << eigenvectors_fromEigen(0, 0) << ", "
                              << eigenvectors_fromEigen(1, 0) << ", " << eigenvectors_fromEigen(2, 0) << std::endl;
    LogDebug("TrackstersPCA") << "Original sigmas:          " << sigmas[0] << ", " << sigmas[1] << ", " << sigmas[2]
                              << std::endl;
    LogDebug("TrackstersPCA") << "SigmasEigen in PCA space: " << sigmasEigen[2] << ", " << sigmasEigen[1] << ", "
                              << sigmasEigen[0] << std::endl;
    LogDebug("TrackstersPCA") << "covM:     \n" << covM << std::endl;
  }
}


