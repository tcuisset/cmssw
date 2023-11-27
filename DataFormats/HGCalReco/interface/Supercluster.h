// Author: Theo Cuisset - theo.cuisset@polytechnique.edu
// Date: 11/2023

#ifndef DataFormats_HGCalReco_Supercluster_h
#define DataFormats_HGCalReco_Supercluster_h

#include <vector>

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"

namespace ticl {
  typedef edm::RefVector<std::vector<Trackster>> Supercluster; // EM supercluster : list of indices into CLUE3D trackster collection. FIrst one is the seed
  typedef std::vector<Supercluster> SuperclusteringResult; // List of superclusters
  typedef unsigned short SuperclusteringDNNScoreValuePacked; // DNN score mapped from float [0, 1] to integer
  /* the DNN score is mapped from float in [0, 1] to unsigned short for storage space. Divide by 2**16=65536 to get the float value back
  Use std::numeric_limits<SuperclusteringDNNScoreValuePacked>::max()
  tuple[0] = seed index, tuple[1] = candidate index, tuple[2] = score (as integer)
  */
  typedef std::vector<std::tuple<unsigned int, unsigned int, SuperclusteringDNNScoreValuePacked>> SuperclusteringDNNScore; 
}


#endif
