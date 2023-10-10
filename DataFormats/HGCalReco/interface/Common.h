#ifndef DataFormats_HGCalReco_Common_h
#define DataFormats_HGCalReco_Common_h

#include <vector>
#include <array>
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
  typedef std::vector<std::vector<std::size_t>> SuperclusteringResult;
  typedef unsigned short SuperclusteringDNNScoreValuePacked;
  // the DNN score is mapped from float in [0, 1] to unsigned char for storage space. Divide by 2**16=65536 to get the float value back
  // Use std::numeric_limits<SuperclusteringDNNScoreValuePacked>::max()
  typedef std::vector<std::tuple<std::size_t, std::size_t, SuperclusteringDNNScoreValuePacked>> SuperclusteringDNNScore; 
}  // namespace ticl

#endif  // DataFormats_HGCalReco_Common_h
