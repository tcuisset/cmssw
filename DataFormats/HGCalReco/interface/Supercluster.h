// Author: Theo Cuisset - theo.cuisset@polytechnique.edu
// Date: 11/2023

#ifndef DataFormats_HGCalReco_Supercluster_h
#define DataFormats_HGCalReco_Supercluster_h

#include <vector>

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"

namespace ticl {
  typedef edm::RefVector<std::vector<Trackster>> Supercluster;
  typedef std::vector<Supercluster> SuperclusteringResult;
  typedef unsigned short SuperclusteringDNNScoreValuePacked;
  // the DNN score is mapped from float in [0, 1] to unsigned char for storage space. Divide by 2**16=65536 to get the float value back
  // Use std::numeric_limits<SuperclusteringDNNScoreValuePacked>::max()
  typedef std::vector<std::tuple<unsigned int, unsigned int, SuperclusteringDNNScoreValuePacked>> SuperclusteringDNNScore; 
}


#endif
