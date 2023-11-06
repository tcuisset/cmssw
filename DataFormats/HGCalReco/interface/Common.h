#ifndef DataFormats_HGCalReco_Common_h
#define DataFormats_HGCalReco_Common_h

#include <RtypesCore.h>

#include <vector>
#include <array>
#include <tuple>
#include <cstdint>

namespace ticl {
  struct TileConstants {
    static constexpr float minEta = 1.5f;
    static constexpr float maxEta = 3.2f;
    static constexpr int nEtaBins = 34;
    static constexpr int nPhiBins = 126;
    static constexpr int nLayers = 104;
    static constexpr int iterations = 4;
    static constexpr int nBins = nEtaBins * nPhiBins;
  };

  struct TileConstantsHFNose {
    static constexpr float minEta = 3.0f;
    static constexpr float maxEta = 4.2f;
    static constexpr int nEtaBins = 24;
    static constexpr int nPhiBins = 126;
    static constexpr int nLayers = 16;  // 8x2
    static constexpr int iterations = 4;
    static constexpr int nBins = nEtaBins * nPhiBins;
  };

}  // namespace ticl

namespace ticl {
  typedef std::vector<std::pair<unsigned int, float> > TICLClusterFilterMask;
  typedef std::vector<std::vector<std::size_t>> SuperclusteringResult; // Each inner vector is a list of trackster ids to supercluster
  struct SuperclusteringDNNScoreElement {
    std::size_t tsSeedId_; // Trackster id of the seed
    std::size_t tsCandId_; // Trackster id of the candidate to be superclustered with the seed
    Float16_t dnnScore_; // score of the dnn

    SuperclusteringDNNScoreElement(std::size_t tsSeedId=0, std::size_t tsCandId=0, float dnnScore=0) :
      tsSeedId_(tsSeedId), tsCandId_(tsCandId), dnnScore_(dnnScore)
    {}
  };
  typedef std::vector<SuperclusteringDNNScoreElement> SuperclusteringDNNScore; 
}  // namespace ticl

#endif  // DataFormats_HGCalReco_Common_h
