// Author: Theo Cuisset - theo.cuisset@polytechnique.edu
// Date: 11/2023

#ifndef __RecoHGCal_TICL_SuperclusteringDNNInputs_H__
#define __RecoHGCal_TICL_SuperclusteringDNNInputs_H__

#include <array>
#include <algorithm>
#include <memory>
#include <string>

#include <Math/Vector2D.h>
#include <Math/Vector3D.h>
#include <Math/Rotation3D.h>
#include <Math/VectorUtil.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"


// Abstract base class for DNN input preparation.
class AbstractDNNInput {
public:
  virtual ~AbstractDNNInput() = default;
  virtual unsigned int featureCount() const = 0;
  // Default implementation of featureNames. Override for custom names
  virtual std::vector<std::string> featureNames() const { 
    std::vector<std::string> defaultNames;
    defaultNames.reserve(featureCount());
    for (unsigned int i = 1; i <= featureCount(); i++) {
        defaultNames.push_back(std::string("nb_") + std::to_string(i));
    }
    return defaultNames;
  }
/* Fill features in input tensor from pair of seed and candidate tracksters.
  inputTensor is a 2D tensor. This function should fill inputTensor(batchIndex, 0...featureCount) with features 
  */
  virtual void fillTensor(tensorflow::Tensor& inputTensor, int batchIndex, ticl::Trackster const& ts_base, ticl::Trackster const& ts_toCluster) = 0;
  virtual std::vector<float> computeVector(ticl::Trackster const& ts_base, ticl::Trackster const& ts_toCluster) = 0;
};


// Helper class to create fillTensor function for different feature lengths
template<typename Concrete, unsigned int Nfeatures>
class DNNInputHelper : public AbstractDNNInput {
    virtual unsigned int featureCount() const { return Nfeatures; }

    virtual void fillTensor(tensorflow::Tensor& inputTensor, int batchIndex, ticl::Trackster const& ts_base, ticl::Trackster const& ts_toCluster) override {
        assert(inputTensor.dims() == 2 && inputTensor.dim_size(1) == Nfeatures);
        assert(batchIndex < inputTensor.dim_size(0));
        auto eigenTensor = inputTensor.tensor<float, 2>();
        auto features = static_cast<Concrete*>(this)->computeFeatures(ts_base, ts_toCluster);
        for (unsigned int i = 0; i < features.size(); i++) {
            eigenTensor(batchIndex, i) = features[i];
        }
    }
    virtual std::vector<float> computeVector(ticl::Trackster const& ts_base, ticl::Trackster const& ts_toCluster) override {
        auto features = static_cast<Concrete*>(this)->computeFeatures(ts_base, ts_toCluster);
        return std::vector<float>(features.begin(), features.end());
    }
};



/*
DNN by Alessandro Tarabini comes in 2 versions with the following observables : 
 'v1_dEtadPhi': ['DeltaEta', 'DeltaPhi', 'multi_en', 'multi_eta', 'multi_pt', 'seedEta','seedPhi','seedEn', 'seedPt'],
  'v2_dEtadPhi': ['DeltaEta', 'DeltaPhi', 'multi_en', 'multi_eta', 'multi_pt', 'seedEta','seedPhi','seedEn', 'seedPt', 'theta', 'theta_xz_seedFrame', 'theta_yz_seedFrame', 'theta_xy_cmsFrame', 'theta_yz_cmsFrame', 'theta_xz_cmsFrame', 'explVar', 'explVarRatio'],
*/

/* First version of DNN by Alessandro Tarabini. Meant as a DNN equivalent of Moustache algorithm (superclustering algo in ECAL)
Uses features : ['DeltaEta', 'DeltaPhi', 'multi_en', 'multi_eta', 'multi_pt', 'seedEta','seedPhi','seedEn', 'seedPt']
*/
class DNNInputAlessandroV1 : public DNNInputHelper<DNNInputAlessandroV1, 9> {
public:
  std::array<float, 9> computeFeatures(ticl::Trackster const& ts_base, ticl::Trackster const& ts_toCluster) {
    /*  We use the barycenter for most of the variables below as that is what seems to have been used by Alessandro Tarabini, 
      but using PCA might be better. 
     (It would need retraining of the DNN)
    */
    return {{
        std::abs(ts_toCluster.barycenter().Eta()) - std::abs(ts_base.barycenter().Eta()), //DeltaEtaBaryc
        ts_toCluster.barycenter().Phi() - ts_base.barycenter().phi(), //DeltaPhiBaryc
        ts_toCluster.raw_energy(), //multi_en
        ts_toCluster.barycenter().Eta(), //multi_eta
        (ts_toCluster.raw_energy() * std::sin(ts_toCluster.barycenter().Theta())), //multi_pt
        ts_base.barycenter().Eta(), //seedEta
        ts_base.barycenter().Phi(), //seedPhi
        ts_base.raw_energy(), //seedEn
        (ts_base.raw_energy() * std::sin(ts_toCluster.barycenter().Theta())), //seedPt
    }};
  }

  std::vector<std::string> featureNames() const override {
    return {"DeltaEtaBaryc", "DeltaPhiBaryc", "multi_en", "multi_eta", "multi_pt", "seedEta", "seedPhi", "seedEn", "seedPt"};
  }
};


// Helper functions for angles. Adapted from ROOT (3D vectors -> 2D vectors)
template <class Vector1, class Vector2>
float CosTheta2D( const Vector1 &  v1, const Vector2  & v2) {
  float arg;
  float v1_r2 = v1.X()*v1.X() + v1.Y()*v1.Y();
  float v2_r2 = v2.X()*v2.X() + v2.Y()*v2.Y();
  float ptot2 = v1_r2*v2_r2;
  if(ptot2 <= 0) {
      arg = 0.0;
  }else{
      float pdot = v1.X()*v2.X() + v1.Y()*v2.Y();
      using std::sqrt;
      arg = pdot/sqrt(ptot2);
      if(arg >  1.0) arg =  1.0;
      if(arg < -1.0) arg = -1.0;
  }
  return arg;
}
template <class Vector1, class Vector2>
inline float Angle2D( const  Vector1 & v1, const Vector2 & v2) {
  return static_cast<float>(std::acos( CosTheta2D(v1, v2) ));
}

/* Second version of DNN by Alessandro Tarabini. 
Uses features : ['DeltaEta', 'DeltaPhi', 'multi_en', 'multi_eta', 'multi_pt', 'seedEta','seedPhi','seedEn', 'seedPt', 'theta', 'theta_xz_seedFrame', 'theta_yz_seedFrame', 'theta_xy_cmsFrame', 'theta_yz_cmsFrame', 'theta_xz_cmsFrame', 'explVar', 'explVarRatio']
*/
class DNNInputAlessandroV2 : public DNNInputHelper<DNNInputAlessandroV2, 17> {
public:
  std::array<float, 17> computeFeatures(ticl::Trackster const& ts_base, ticl::Trackster const& ts_toCluster) {
    using ROOT::Math::VectorUtil::Angle;
    using ROOT::Math::XYZVectorF;
    using ROOT::Math::XYVectorF;
    XYZVectorF const& pca_seed_cmsFrame(ts_base.eigenvectors(0));
    XYZVectorF const& pca_cand_cmsFrame(ts_toCluster.eigenvectors(0));
    XYZVectorF xs(pca_seed_cmsFrame.Cross(XYZVectorF(0, 0, 1)).Unit());
    ROOT::Math::Rotation3D rot(
        xs, 
        xs.Cross(pca_seed_cmsFrame).Unit(),
        pca_seed_cmsFrame);
      
    XYZVectorF pca_cand_seedFrame = rot(pca_cand_cmsFrame); // seed coordinates

    float explVar_denominator = std::accumulate(std::begin(ts_toCluster.eigenvalues()), std::end(ts_toCluster.eigenvalues()), 0.f, std::plus<float>());
    float explVarRatio = 0.;
    if (explVar_denominator != 0.) {
        explVarRatio = ts_toCluster.eigenvalues()[0] / explVar_denominator; // explVarRatio
    } else {
      // TODO Study what would be best in this case (or if it could ever happen)
      LogDebug("HGCalTICLSuperclustering") << "Sum of eigenvalues was zero for trackster. Could not compute explained variance ratio.";
    }

    return {{
        std::abs(ts_toCluster.barycenter().Eta()) - std::abs(ts_base.barycenter().Eta()), //DeltaEtaBaryc
        ts_toCluster.barycenter().Phi() - ts_base.barycenter().phi(), //DeltaPhiBaryc
        ts_toCluster.raw_energy(), //multi_en
        ts_toCluster.barycenter().Eta(), //multi_eta
        (ts_toCluster.raw_energy() * std::sin(ts_toCluster.barycenter().Theta())), //multi_pt
        ts_base.barycenter().Eta(), //seedEta
        ts_base.barycenter().Phi(), //seedPhi
        ts_base.raw_energy(), //seedEn
        (ts_base.raw_energy() * std::sin(ts_toCluster.barycenter().Theta())), //seedPt
        static_cast<float>(Angle(pca_cand_cmsFrame, pca_seed_cmsFrame)), // theta : angle between seed and candidate
        Angle2D(XYVectorF(pca_cand_seedFrame.y(), pca_cand_seedFrame.z()), XYVectorF(0, 1)), // theta_xz_seedFrame
        Angle2D(XYVectorF(pca_cand_seedFrame.y(), pca_cand_seedFrame.z()), XYVectorF(0, 1)), // theta_yz_seedFrame
        Angle2D(XYVectorF(pca_cand_cmsFrame.x(), pca_cand_cmsFrame.y()), XYVectorF(pca_seed_cmsFrame.x(), pca_seed_cmsFrame.y())), // theta_xy_cmsFrame
        Angle2D(XYVectorF(pca_cand_cmsFrame.y(), pca_cand_cmsFrame.z()), XYVectorF(pca_seed_cmsFrame.y(), pca_seed_cmsFrame.z())), // theta_yz_cmsFrame
        Angle2D(XYVectorF(pca_cand_cmsFrame.x(), pca_cand_cmsFrame.z()), XYVectorF(pca_seed_cmsFrame.x(), pca_seed_cmsFrame.z())), // theta_xz_cmsFrame
        ts_toCluster.eigenvalues()[0], // explVar
        explVarRatio // explVarRatio
    }};
  }

  std::vector<std::string> featureNames() const override {
    return {"DeltaEtaBaryc", "DeltaPhiBaryc", "multi_en", "multi_eta", "multi_pt", "seedEta", "seedPhi", "seedEn", "seedPt", "theta", "theta_xz_seedFrame", "theta_yz_seedFrame", "theta_xy_cmsFrame", "theta_yz_cmsFrame", "theta_xz_cmsFrame", "explVar", "explVarRatio"};
  }
};


std::unique_ptr<AbstractDNNInput> makeDNNInputFromString(std::string dnnVersion); // defined in SuperclusteringProducer.cc

#endif