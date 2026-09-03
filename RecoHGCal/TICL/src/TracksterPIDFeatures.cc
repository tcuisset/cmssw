#include "RecoHGCal/TICL/interface/TracksterPIDFeatures.h"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "DataFormats/Math/interface/deltaPhi.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/TICLGeomTools.h"

namespace ticl {

  void TracksterPIDFeatures::fillSummary(float* out,
                                         Trackster const& ts,
                                         std::vector<reco::CaloCluster> const& layerClusters,
                                         ticlgeom::Tools const& rhtools) const {
    std::vector<float> layerEnergy(nLayers_, 0.f);
    int firstLayer = nLayers_;
    int lastLayer = -1;
    float weightedLayer = 0.f;
    float weightedLayer2 = 0.f;
    float totalEnergy = 0.f;
    float maximumLayerEnergy = 0.f;
    float totalHits = 0.f;

    for (int k = 0, size = ts.vertices().size(); k < size; ++k) {
      auto const& cluster = layerClusters[ts.vertices(k)];
      if (cluster.hitsAndFractions().empty())
        continue;
      const int layer = rhtools.getLayerWithOffset(cluster.hitsAndFractions()[0].first) - 1;
      if (layer < 0 || layer >= nLayers_)
        continue;
      const float energy = cluster.energy() / ts.vertex_multiplicity(k);
      layerEnergy[layer] += energy;
      totalEnergy += energy;
      weightedLayer += energy * layer;
      weightedLayer2 += energy * layer * layer;
      totalHits += cluster.hitsAndFractions().size();
      firstLayer = std::min(firstLayer, layer);
      lastLayer = std::max(lastLayer, layer);
    }
    for (float energy : layerEnergy)
      maximumLayerEnergy = std::max(maximumLayerEnergy, energy);

    const float inverseEnergy = totalEnergy > 0.f ? 1.f / totalEnergy : 0.f;
    const float meanLayer = weightedLayer * inverseEnergy;
    const float layerVariance = std::max(0.f, weightedLayer2 * inverseEnergy - meanLayer * meanLayer);
    const int split1 = std::min(10, nLayers_);
    const int split2 = std::min(28, nLayers_);
    const float earlyEnergy = std::accumulate(layerEnergy.begin(), layerEnergy.begin() + split1, 0.f);
    const float electromagneticEnergy = std::accumulate(layerEnergy.begin(), layerEnergy.begin() + split2, 0.f);
    const float lateEnergy = std::accumulate(layerEnergy.begin() + split2, layerEnergy.end(), 0.f);

    out[0] = std::log1p(std::max(0.f, static_cast<float>(ts.raw_energy())));
    out[1] = ts.raw_energy() > 0.f ? ts.raw_em_energy() / ts.raw_energy() : 0.f;
    out[2] = std::abs(ts.barycenter().eta());
    out[3] = std::sin(ts.barycenter().phi());
    out[4] = std::cos(ts.barycenter().phi());
    out[5] = std::log1p(static_cast<float>(ts.vertices().size()));
    out[6] = firstLayer < nLayers_ ? static_cast<float>(firstLayer) / nLayers_ : 0.f;
    out[7] = lastLayer >= 0 ? static_cast<float>(lastLayer) / nLayers_ : 0.f;
    out[8] = lastLayer >= firstLayer ? static_cast<float>(lastLayer - firstLayer + 1) / nLayers_ : 0.f;
    out[9] = meanLayer / nLayers_;
    out[10] = std::sqrt(layerVariance) / nLayers_;
    out[11] = maximumLayerEnergy * inverseEnergy;
    out[12] = earlyEnergy * inverseEnergy;
    out[13] = electromagneticEnergy * inverseEnergy;
    out[14] = lateEnergy * inverseEnergy;
    out[15] = ts.vertices().empty() ? 0.f : totalHits / ts.vertices().size();
  }

  void TracksterPIDFeatures::fillLayer(float* out,
                                       Trackster const& ts,
                                       std::vector<reco::CaloCluster> const& layerClusters,
                                       ticlgeom::Tools const& rhtools) const {
    const float normalisation = ts.raw_energy() > 0.f ? 1.f / ts.raw_energy() : 0.f;
    for (int k = 0, size = ts.vertices().size(); k < size; ++k) {
      auto const& cluster = layerClusters[ts.vertices(k)];
      if (cluster.hitsAndFractions().empty())
        continue;
      const int layer = rhtools.getLayerWithOffset(cluster.hitsAndFractions()[0].first) - 1;
      if (layer < 0 || layer >= nLayers_)
        continue;
      const float energy = cluster.energy() / ts.vertex_multiplicity(k);
      const float deta = cluster.eta() - ts.barycenter().eta();
      const float dphi = reco::deltaPhi(cluster.phi(), ts.barycenter().phi());
      float* features = out + static_cast<std::size_t>(layer) * layerFeatureCount;
      features[0] += energy;
      features[1] += 1.f;
      features[2] += energy * deta;
      features[3] += energy * dphi;
      features[4] += energy * (deta * deta + dphi * dphi);
      features[5] += energy * cluster.hitsAndFractions().size();
      features[6] = std::max(features[6], energy);
      features[7] = 1.f;
    }
    for (int layer = 0; layer < nLayers_; ++layer) {
      float* features = out + static_cast<std::size_t>(layer) * layerFeatureCount;
      const float layerEnergy = features[0];
      const float inverseLayerEnergy = layerEnergy > 0.f ? 1.f / layerEnergy : 0.f;
      features[0] *= normalisation;
      features[1] = std::log1p(features[1]);
      features[2] *= inverseLayerEnergy;
      features[3] *= inverseLayerEnergy;
      features[4] = std::sqrt(std::max(0.f, features[4] * inverseLayerEnergy));
      features[5] *= inverseLayerEnergy;
      features[6] *= normalisation;
    }
  }

  void TracksterPIDFeatures::fillCluster(float* out,
                                         Trackster const& ts,
                                         std::vector<reco::CaloCluster> const& layerClusters,
                                         ticlgeom::Tools const& rhtools,
                                         std::vector<int>& order) const {
    order.resize(ts.vertices().size());
    std::iota(order.begin(), order.end(), 0);
    const int selected = std::min<int>(maxClusters_, order.size());
    std::partial_sort(order.begin(), order.begin() + selected, order.end(), [&layerClusters, &ts](int a, int b) {
      return layerClusters[ts.vertices(a)].energy() > layerClusters[ts.vertices(b)].energy();
    });
    const float normalisation = ts.raw_energy() > 0.f ? 1.f / ts.raw_energy() : 0.f;
    for (int index = 0; index < selected; ++index) {
      const int k = order[index];
      auto const& cluster = layerClusters[ts.vertices(k)];
      if (cluster.hitsAndFractions().empty())
        continue;
      const int layer = rhtools.getLayerWithOffset(cluster.hitsAndFractions()[0].first) - 1;
      if (layer < 0 || layer >= nLayers_)
        continue;
      const float energy = cluster.energy() / ts.vertex_multiplicity(k);
      float* features = out + static_cast<std::size_t>(index) * clusterFeatureCount;
      features[0] = energy * normalisation;
      features[1] = std::log1p(std::max(0.f, energy));
      features[2] = static_cast<float>(layer) / nLayers_;
      features[3] = cluster.eta() - ts.barycenter().eta();
      features[4] = reco::deltaPhi(cluster.phi(), ts.barycenter().phi());
      features[5] = std::log1p(static_cast<float>(cluster.hitsAndFractions().size()));
      features[6] = 1.f / ts.vertex_multiplicity(k);
      features[7] = 1.f;
    }
  }

  void TracksterPIDFeatures::fillImage(float* out,
                                       Trackster const& ts,
                                       std::vector<reco::CaloCluster> const& layerClusters,
                                       ticlgeom::Tools const& rhtools,
                                       std::vector<int>& order,
                                       std::vector<int>& seenClusters) const {
    order.resize(ts.vertices().size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&layerClusters, &ts](int a, int b) {
      return layerClusters[ts.vertices(a)].energy() > layerClusters[ts.vertices(b)].energy();
    });
    std::fill(seenClusters.begin(), seenClusters.end(), 0);
    const float normalisation = ts.raw_energy() > 0.f ? 1.f / ts.raw_energy() : 0.f;
    const std::size_t planeSize = static_cast<std::size_t>(nLayers_) * maxClusters_;
    for (int k : order) {
      auto const& cluster = layerClusters[ts.vertices(k)];
      if (cluster.hitsAndFractions().empty())
        continue;
      const int layer = rhtools.getLayerWithOffset(cluster.hitsAndFractions()[0].first) - 1;
      if (layer < 0 || layer >= nLayers_ || seenClusters[layer] >= maxClusters_)
        continue;
      const std::size_t cell = static_cast<std::size_t>(layer) * maxClusters_ + seenClusters[layer]++;
      out[cell] = cluster.energy() / ts.vertex_multiplicity(k) * normalisation;
      out[planeSize + cell] = cluster.eta() - ts.barycenter().eta();
      out[2 * planeSize + cell] = reco::deltaPhi(cluster.phi(), ts.barycenter().phi());
    }
  }

}  // namespace ticl
