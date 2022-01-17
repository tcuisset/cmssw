/* Sushil Dubey, Shashi Dugad, TIFR, July 2017
 *
 * File Name: RawToClusterGPU.cu
 * Description: It converts Raw data into Digi Format on GPU
 * Finaly the Output of RawToDigi data is given to pixelClusterizer
**/

// C++ includes
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

// CUDA includes
#include <cuda_runtime.h>

// CMSSW includes
#include "CUDADataFormats/SiPixelCluster/interface/gpuClusteringConstants.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelROCsStatusAndMapping.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/SiPixelDigi/interface/SiPixelDigiConstants.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "HeterogeneousCore/CUDAUtilities/interface/device_unique_ptr.h"
#include "HeterogeneousCore/CUDAUtilities/interface/host_unique_ptr.h"
#include "RecoLocalTracker/SiPixelClusterizer/plugins/gpuCalibPixel.h"
#include "RecoLocalTracker/SiPixelClusterizer/plugins/gpuClusterChargeCut.h"
#include "RecoLocalTracker/SiPixelClusterizer/plugins/gpuClustering.h"
// local includes
#include "SiPixelRawToClusterGPUKernel.h"

namespace pixelgpudetails {

  SiPixelRawToClusterGPUKernel::WordFedAppender::WordFedAppender(uint32_t maxFedWords) {
    word_ = cms::cuda::make_host_noncached_unique<unsigned int[]>(maxFedWords, cudaHostAllocWriteCombined);
    fedId_ = cms::cuda::make_host_noncached_unique<unsigned char[]>(maxFedWords, cudaHostAllocWriteCombined);
  }

  void SiPixelRawToClusterGPUKernel::WordFedAppender::initializeWordFed(int fedId,
                                                                        unsigned int wordCounterGPU,
                                                                        const cms_uint32_t *src,
                                                                        unsigned int length) {
    std::memcpy(word_.get() + wordCounterGPU, src, sizeof(cms_uint32_t) * length);
    std::memset(fedId_.get() + wordCounterGPU / 2, fedId - FEDNumbering::MINSiPixeluTCAFEDID, length / 2);
  }

  ////////////////////

  __device__ bool isBarrel(uint32_t rawId) {
    return (PixelSubdetector::PixelBarrel == ((rawId >> DetId::kSubdetOffset) & DetId::kSubdetMask));
  }

  __device__ pixelgpudetails::DetIdGPU getRawId(const SiPixelROCsStatusAndMapping *cablingMap,
                                                uint8_t fed,
                                                uint32_t link,
                                                uint32_t roc) {
    uint32_t index = fed * MAX_LINK * MAX_ROC + (link - 1) * MAX_ROC + roc;
    pixelgpudetails::DetIdGPU detId = {
        cablingMap->rawId[index], cablingMap->rocInDet[index], cablingMap->moduleId[index]};
    return detId;
  }

  //reference http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_9_2_0/doc/html/dd/d31/FrameConversion_8cc_source.html
  //http://cmslxr.fnal.gov/source/CondFormats/SiPixelObjects/src/PixelROC.cc?v=CMSSW_9_2_0#0071
  // Convert local pixel to pixelgpudetails::global pixel
  __device__ pixelgpudetails::Pixel frameConversion(
      bool bpix, int side, uint32_t layer, uint32_t rocIdInDetUnit, pixelgpudetails::Pixel local) {
    int slopeRow = 0, slopeCol = 0;
    int rowOffset = 0, colOffset = 0;

    if (bpix) {
      if (side == -1 && layer != 1) {  // -Z side: 4 non-flipped modules oriented like 'dddd', except Layer 1
        if (rocIdInDetUnit < 8) {
          slopeRow = 1;
          slopeCol = -1;
          rowOffset = 0;
          colOffset = (8 - rocIdInDetUnit) * pixelgpudetails::numColsInRoc - 1;
        } else {
          slopeRow = -1;
          slopeCol = 1;
          rowOffset = 2 * pixelgpudetails::numRowsInRoc - 1;
          colOffset = (rocIdInDetUnit - 8) * pixelgpudetails::numColsInRoc;
        }       // if roc
      } else {  // +Z side: 4 non-flipped modules oriented like 'pppp', but all 8 in layer1
        if (rocIdInDetUnit < 8) {
          slopeRow = -1;
          slopeCol = 1;
          rowOffset = 2 * pixelgpudetails::numRowsInRoc - 1;
          colOffset = rocIdInDetUnit * pixelgpudetails::numColsInRoc;
        } else {
          slopeRow = 1;
          slopeCol = -1;
          rowOffset = 0;
          colOffset = (16 - rocIdInDetUnit) * pixelgpudetails::numColsInRoc - 1;
        }
      }

    } else {             // fpix
      if (side == -1) {  // pannel 1
        if (rocIdInDetUnit < 8) {
          slopeRow = 1;
          slopeCol = -1;
          rowOffset = 0;
          colOffset = (8 - rocIdInDetUnit) * pixelgpudetails::numColsInRoc - 1;
        } else {
          slopeRow = -1;
          slopeCol = 1;
          rowOffset = 2 * pixelgpudetails::numRowsInRoc - 1;
          colOffset = (rocIdInDetUnit - 8) * pixelgpudetails::numColsInRoc;
        }
      } else {  // pannel 2
        if (rocIdInDetUnit < 8) {
          slopeRow = 1;
          slopeCol = -1;
          rowOffset = 0;
          colOffset = (8 - rocIdInDetUnit) * pixelgpudetails::numColsInRoc - 1;
        } else {
          slopeRow = -1;
          slopeCol = 1;
          rowOffset = 2 * pixelgpudetails::numRowsInRoc - 1;
          colOffset = (rocIdInDetUnit - 8) * pixelgpudetails::numColsInRoc;
        }

      }  // side
    }

    uint32_t gRow = rowOffset + slopeRow * local.row;
    uint32_t gCol = colOffset + slopeCol * local.col;
    // inside frameConversion row: gRow, column: gCol
    pixelgpudetails::Pixel global = {gRow, gCol};
    return global;
  }

  // error decoding and handling copied from EventFilter/SiPixelRawToDigi/src/ErrorChecker.cc
  __device__ uint8_t conversionError(uint8_t fedId, uint8_t status, bool debug = false) {
    uint8_t errorType = 0;

    switch (status) {
      case (1): {
        if (debug)
          printf("Error in Fed: %i, invalid channel Id (errorType = 35\n)", fedId);
        errorType = 35;
        break;
      }
      case (2): {
        if (debug)
          printf("Error in Fed: %i, invalid ROC Id (errorType = 36)\n", fedId);
        errorType = 36;
        break;
      }
      case (3): {
        if (debug)
          printf("Error in Fed: %i, invalid dcol/pixel value (errorType = 37)\n", fedId);
        errorType = 37;
        break;
      }
      case (4): {
        if (debug)
          printf("Error in Fed: %i, dcol/pixel read out of order (errorType = 38)\n", fedId);
        errorType = 38;
        break;
      }
      default:
        if (debug)
          printf("Cabling check returned unexpected result, status = %i\n", status);
    };

    return errorType;
  }

  __device__ bool rocRowColIsValid(uint32_t rocRow, uint32_t rocCol) {
    /// row and column in ROC representation
    return ((rocRow < pixelgpudetails::numRowsInRoc) & (rocCol < pixelgpudetails::numColsInRoc));
  }

  __device__ bool dcolIsValid(uint32_t dcol, uint32_t pxid) { return ((dcol < 26) & (2 <= pxid) & (pxid < 162)); }

  // error decoding and handling copied from EventFilter/SiPixelRawToDigi/src/ErrorChecker.cc
  __device__ uint8_t checkROC(uint32_t errorWord,
                              uint8_t fedId,
                              uint32_t link,
                              const SiPixelROCsStatusAndMapping *cablingMap,
                              bool debug = false) {
    uint8_t errorType = (errorWord >> sipixelconstants::ROC_shift) & sipixelconstants::ERROR_mask;
    if (errorType < 25)
      return 0;
    bool errorFound = false;

    switch (errorType) {
      case (25): {
        errorFound = true;
        uint32_t index = fedId * MAX_LINK * MAX_ROC + (link - 1) * MAX_ROC + 1;
        if (index > 1 && index <= cablingMap->size) {
          if (!(link == cablingMap->link[index] && 1 == cablingMap->roc[index]))
            errorFound = false;
        }
        if (debug and errorFound)
          printf("Invalid ROC = 25 found (errorType = 25)\n");
        break;
      }
      case (26): {
        if (debug)
          printf("Gap word found (errorType = 26)\n");
        errorFound = true;
        break;
      }
      case (27): {
        if (debug)
          printf("Dummy word found (errorType = 27)\n");
        errorFound = true;
        break;
      }
      case (28): {
        if (debug)
          printf("Error fifo nearly full (errorType = 28)\n");
        errorFound = true;
        break;
      }
      case (29): {
        if (debug)
          printf("Timeout on a channel (errorType = 29)\n");
        if ((errorWord >> sipixelconstants::OMIT_ERR_shift) & sipixelconstants::OMIT_ERR_mask) {
          if (debug)
            printf("...first errorType=29 error, this gets masked out\n");
        }
        errorFound = true;
        break;
      }
      case (30): {
        if (debug)
          printf("TBM error trailer (errorType = 30)\n");
        int stateMatch_bits = 4;
        int stateMatch_shift = 8;
        uint32_t stateMatch_mask = ~(~uint32_t(0) << stateMatch_bits);
        int stateMatch = (errorWord >> stateMatch_shift) & stateMatch_mask;
        if (stateMatch != 1 && stateMatch != 8) {
          if (debug)
            printf("FED error 30 with unexpected State Bits (errorType = 30)\n");
        }
        if (stateMatch == 1)
          errorType = 40;  // 1=Overflow -> 40, 8=number of ROCs -> 30
        errorFound = true;
        break;
      }
      case (31): {
        if (debug)
          printf("Event number error (errorType = 31)\n");
        errorFound = true;
        break;
      }
      default:
        errorFound = false;
    };

    return errorFound ? errorType : 0;
  }

  // error decoding and handling copied from EventFilter/SiPixelRawToDigi/src/ErrorChecker.cc
  __device__ uint32_t getErrRawID(uint8_t fedId,
                                  uint32_t errWord,
                                  uint32_t errorType,
                                  const SiPixelROCsStatusAndMapping *cablingMap,
                                  bool debug = false) {
    uint32_t rID = 0xffffffff;

    switch (errorType) {
      case 25:
      case 30:
      case 31:
      case 36:
      case 40: {
        uint32_t roc = 1;
        uint32_t link = sipixelconstants::getLink(errWord);
        uint32_t rID_temp = getRawId(cablingMap, fedId, link, roc).rawId;
        if (rID_temp != gpuClustering::invalidModuleId)
          rID = rID_temp;
        break;
      }
      case 29: {
        int chanNmbr = 0;
        const int DB0_shift = 0;
        const int DB1_shift = DB0_shift + 1;
        const int DB2_shift = DB1_shift + 1;
        const int DB3_shift = DB2_shift + 1;
        const int DB4_shift = DB3_shift + 1;
        const uint32_t DataBit_mask = ~(~uint32_t(0) << 1);

        int CH1 = (errWord >> DB0_shift) & DataBit_mask;
        int CH2 = (errWord >> DB1_shift) & DataBit_mask;
        int CH3 = (errWord >> DB2_shift) & DataBit_mask;
        int CH4 = (errWord >> DB3_shift) & DataBit_mask;
        int CH5 = (errWord >> DB4_shift) & DataBit_mask;
        int BLOCK_bits = 3;
        int BLOCK_shift = 8;
        uint32_t BLOCK_mask = ~(~uint32_t(0) << BLOCK_bits);
        int BLOCK = (errWord >> BLOCK_shift) & BLOCK_mask;
        int localCH = 1 * CH1 + 2 * CH2 + 3 * CH3 + 4 * CH4 + 5 * CH5;
        if (BLOCK % 2 == 0)
          chanNmbr = (BLOCK / 2) * 9 + localCH;
        else
          chanNmbr = ((BLOCK - 1) / 2) * 9 + 4 + localCH;
        if ((chanNmbr < 1) || (chanNmbr > 36))
          break;  // signifies unexpected result

        uint32_t roc = 1;
        uint32_t link = chanNmbr;
        uint32_t rID_temp = getRawId(cablingMap, fedId, link, roc).rawId;
        if (rID_temp != gpuClustering::invalidModuleId)
          rID = rID_temp;
        break;
      }
      case 37:
      case 38: {
        uint32_t roc = sipixelconstants::getROC(errWord);
        uint32_t link = sipixelconstants::getLink(errWord);
        uint32_t rID_temp = getRawId(cablingMap, fedId, link, roc).rawId;
        if (rID_temp != gpuClustering::invalidModuleId)
          rID = rID_temp;
        break;
      }
      default:
        break;
    };

    return rID;
  }

  // Kernel to perform Raw to Digi conversion
  __global__ void RawToDigi_kernel(const SiPixelROCsStatusAndMapping *cablingMap,
                                   const unsigned char *modToUnp,
                                   const uint32_t wordCounter,
                                   const uint32_t *word,
                                   const uint8_t *fedIds,
                                   uint16_t *xx,
                                   uint16_t *yy,
                                   uint16_t *adc,
                                   uint32_t *pdigi,
                                   uint32_t *rawIdArr,
                                   uint16_t *moduleId,
                                   cms::cuda::SimpleVector<SiPixelErrorCompact> *err,
                                   bool useQualityInfo,
                                   bool includeErrors,
                                   bool debug) {
    //if (threadIdx.x==0) printf("Event: %u blockIdx.x: %u start: %u end: %u\n", eventno, blockIdx.x, begin, end);

    int32_t first = threadIdx.x + blockIdx.x * blockDim.x;
    for (int32_t iloop = first, nend = wordCounter; iloop < nend; iloop += blockDim.x * gridDim.x) {
      auto gIndex = iloop;
      xx[gIndex] = 0;
      yy[gIndex] = 0;
      adc[gIndex] = 0;
      bool skipROC = false;

      uint8_t fedId = fedIds[gIndex / 2];  // +1200;

      // initialize (too many coninue below)
      pdigi[gIndex] = 0;
      rawIdArr[gIndex] = 0;
      moduleId[gIndex] = gpuClustering::invalidModuleId;

      uint32_t ww = word[gIndex];  // Array containing 32 bit raw data
      if (ww == 0) {
        // 0 is an indicator of a noise/dead channel, skip these pixels during clusterization
        continue;
      }

      uint32_t link = sipixelconstants::getLink(ww);  // Extract link
      uint32_t roc = sipixelconstants::getROC(ww);    // Extract Roc in link
      pixelgpudetails::DetIdGPU detId = getRawId(cablingMap, fedId, link, roc);

      uint8_t errorType = checkROC(ww, fedId, link, cablingMap, debug);
      skipROC = (roc < pixelgpudetails::maxROCIndex) ? false : (errorType != 0);
      if (includeErrors and skipROC) {
        uint32_t rID = getErrRawID(fedId, ww, errorType, cablingMap, debug);
        err->push_back(SiPixelErrorCompact{rID, ww, errorType, fedId});
        continue;
      }

      uint32_t rawId = detId.rawId;
      uint32_t rocIdInDetUnit = detId.rocInDet;
      bool barrel = isBarrel(rawId);

      uint32_t index = fedId * MAX_LINK * MAX_ROC + (link - 1) * MAX_ROC + roc;
      if (useQualityInfo) {
        skipROC = cablingMap->badRocs[index];
        if (skipROC)
          continue;
      }
      skipROC = modToUnp[index];
      if (skipROC)
        continue;

      uint32_t layer = 0;
      int side = 0, panel = 0, module = 0;

      if (barrel) {
        layer = (rawId >> pixelgpudetails::layerStartBit) & pixelgpudetails::layerMask;
        module = (rawId >> pixelgpudetails::moduleStartBit) & pixelgpudetails::moduleMask;
        side = (module < 5) ? -1 : 1;
      } else {
        // endcap ids
        layer = 0;
        panel = (rawId >> pixelgpudetails::panelStartBit) & pixelgpudetails::panelMask;
        side = (panel == 1) ? -1 : 1;
      }

      // ***special case of layer to 1 be handled here
      pixelgpudetails::Pixel localPix;
      if (layer == 1) {
        uint32_t col = sipixelconstants::getCol(ww);
        uint32_t row = sipixelconstants::getRow(ww);
        localPix.row = row;
        localPix.col = col;
        if (includeErrors) {
          if (not rocRowColIsValid(row, col)) {
            uint8_t error = conversionError(fedId, 3, debug);  //use the device function and fill the arrays
            err->push_back(SiPixelErrorCompact{rawId, ww, error, fedId});
            if (debug)
              printf("BPIX1  Error status: %i\n", error);
            continue;
          }
        }
      } else {
        // ***conversion rules for dcol and pxid
        uint32_t dcol = sipixelconstants::getDCol(ww);
        uint32_t pxid = sipixelconstants::getPxId(ww);
        uint32_t row = pixelgpudetails::numRowsInRoc - pxid / 2;
        uint32_t col = dcol * 2 + pxid % 2;
        localPix.row = row;
        localPix.col = col;
        if (includeErrors and not dcolIsValid(dcol, pxid)) {
          uint8_t error = conversionError(fedId, 3, debug);
          err->push_back(SiPixelErrorCompact{rawId, ww, error, fedId});
          if (debug)
            printf("Error status: %i %d %d %d %d\n", error, dcol, pxid, fedId, roc);
          continue;
        }
      }

      pixelgpudetails::Pixel globalPix = frameConversion(barrel, side, layer, rocIdInDetUnit, localPix);
      xx[gIndex] = globalPix.row;  // origin shifting by 1 0-159
      yy[gIndex] = globalPix.col;  // origin shifting by 1 0-415
      adc[gIndex] = sipixelconstants::getADC(ww);
      pdigi[gIndex] = pixelgpudetails::pack(globalPix.row, globalPix.col, adc[gIndex]);
      moduleId[gIndex] = detId.moduleId;
      rawIdArr[gIndex] = rawId;
    }  // end of loop (gIndex < end)

  }  // end of Raw to Digi kernel
  template <bool isPhase2>
  __global__ void fillHitsModuleStart(uint32_t const *__restrict__ clusInModule,
                                      uint32_t *__restrict__ moduleStart,
                                      uint32_t const *__restrict__ nModules,
                                      uint32_t *__restrict__ nModules_Clusters) {
    constexpr int nMaxModules = isPhase2 ? phase2PixelTopology::numberOfModules : phase1PixelTopology::numberOfModules;
    constexpr int startBPIX2 = isPhase2 ? phase2PixelTopology::layerStart[1] : phase1PixelTopology::layerStart[1];

    assert(nMaxModules < phase2PixelTopology::numberOfModules);
    assert(startBPIX2 < nMaxModules);
    assert(nMaxModules < 4096);  // easy to extend at least till 32*1024
    assert(nMaxModules > 1024);

    assert(1 == gridDim.x);
    assert(0 == blockIdx.x);

    int first = threadIdx.x;

    // limit to MaxHitsInModule;
    for (int i = first, iend = nMaxModules; i < iend; i += blockDim.x) {
      moduleStart[i + 1] = std::min(gpuClustering::maxHitsInModule(), clusInModule[i]);
    }

    __shared__ uint32_t ws[64];
    cms::cuda::blockPrefixScan(moduleStart + 1, moduleStart + 1, 1024, ws);
    constexpr int lastModules = isPhase2 ? 1024 : nMaxModules - 1024;
    cms::cuda::blockPrefixScan(moduleStart + 1024 + 1, moduleStart + 1024 + 1, lastModules, ws);

    if constexpr (isPhase2) {
      cms::cuda::blockPrefixScan(moduleStart + 2048 + 1, moduleStart + 2048 + 1, 1024, ws);
      cms::cuda::blockPrefixScan(moduleStart + 3072 + 1, moduleStart + 3072 + 1, nMaxModules - 3072, ws);
    }

    for (int i = first + 1025, iend = isPhase2 ? 2049 : nMaxModules + 1; i < iend; i += blockDim.x) {
      moduleStart[i] += moduleStart[1024];
    }
    __syncthreads();

    if constexpr (isPhase2) {
      for (int i = first + 2049, iend = 3073; i < iend; i += blockDim.x) {
        moduleStart[i] += moduleStart[2048];
      }
      __syncthreads();
      for (int i = first + 3073, iend = nMaxModules + 1; i < iend; i += blockDim.x) {
        moduleStart[i] += moduleStart[3072];
      }
      __syncthreads();
    }

    if (threadIdx.x == 0) {
      // copy the number of modules
      nModules_Clusters[0] = *nModules;
      // last element holds the number of all clusters
      nModules_Clusters[1] = moduleStart[nMaxModules];
      // element 96 is the start of BPIX2 (i.e. the number of clusters in BPIX1)
      nModules_Clusters[2] = moduleStart[startBPIX2];
    }

#ifdef GPU_DEBUG
    uint16_t maxH = isPhase2 ? 3027 : 1024;
    assert(0 == moduleStart[0]);
    auto c0 = std::min(gpuClustering::maxHitsInModule(), clusInModule[0]);
    assert(c0 == moduleStart[1]);
    assert(moduleStart[maxH] >= moduleStart[maxH - 1]);
    assert(moduleStart[maxH + 1] >= moduleStart[maxH]);
    assert(moduleStart[nMaxModules] >= moduleStart[maxH + 1]);

    constexpr int startFP1 =
        isPhase2 ? phase2PixelTopology::numberOfModulesInBarrel : phase1PixelTopology::numberOfModulesInBarrel;
    constexpr int startLastFwd = isPhase2 ? phase2PixelTopology::layerStart[phase2PixelTopology::numberOfLayers]
                                          : phase1PixelTopology::layerStart[phase1PixelTopology::numberOfLayers];
    for (int i = first, iend = nMaxModules + 1; i < iend; i += blockDim.x) {
      if (0 != i)
        assert(moduleStart[i] >= moduleStart[i - i]);
      // [BPX1, BPX2, BPX3, BPX4,  FP1,  FP2,  FP3,  FN1,  FN2,  FN3, LAST_VALID]
      // [   0,   96,  320,  672, 1184, 1296, 1408, 1520, 1632, 1744,       1856]
      if (i == startBPIX2 || i == startFP1 || i == startLastFwd || i == nMaxModules)
        printf("moduleStart %d %d\n", i, moduleStart[i]);
    }

#endif
  }

  // Interface to outside
  void SiPixelRawToClusterGPUKernel::makeClustersAsync(bool isRun2,
                                                       const SiPixelClusterThresholds clusterThresholds,
                                                       const SiPixelROCsStatusAndMapping *cablingMap,
                                                       const unsigned char *modToUnp,
                                                       const SiPixelGainForHLTonGPU *gains,
                                                       const WordFedAppender &wordFed,
                                                       SiPixelFormatterErrors &&errors,
                                                       const uint32_t wordCounter,
                                                       const uint32_t fedCounter,
                                                       const uint32_t maxFedWords,
                                                       bool useQualityInfo,
                                                       bool includeErrors,
                                                       bool debug,
                                                       cudaStream_t stream) {
    // we're not opting for calling this function in case of early events
    assert(wordCounter != 0);
    nDigis = wordCounter;

#ifdef GPU_DEBUG
    std::cout << "decoding " << wordCounter << " digis. Max is " << maxFedWords << std::endl;
#endif

    // since wordCounter != 0 we're not allocating 0 bytes,
    digis_d = SiPixelDigisCUDA(wordCounter, stream);
    if (includeErrors) {
      digiErrors_d = SiPixelDigiErrorsCUDA(wordCounter, std::move(errors), stream);
    }
    clusters_d = SiPixelClustersCUDA(phase1PixelTopology::numberOfModules, stream);

    // Begin Raw2Digi block
    {
      const int threadsPerBlock = 512;
      const int blocks = (wordCounter + threadsPerBlock - 1) / threadsPerBlock;  // fill it all

      assert(0 == wordCounter % 2);
      // wordCounter is the total no of words in each event to be trasfered on device
      auto word_d = cms::cuda::make_device_unique<uint32_t[]>(wordCounter, stream);
      auto fedId_d = cms::cuda::make_device_unique<uint8_t[]>(wordCounter, stream);

      cudaCheck(
          cudaMemcpyAsync(word_d.get(), wordFed.word(), wordCounter * sizeof(uint32_t), cudaMemcpyDefault, stream));
      cudaCheck(cudaMemcpyAsync(
          fedId_d.get(), wordFed.fedId(), wordCounter * sizeof(uint8_t) / 2, cudaMemcpyDefault, stream));

      // Launch rawToDigi kernel
      RawToDigi_kernel<<<blocks, threadsPerBlock, 0, stream>>>(
          cablingMap,
          modToUnp,
          wordCounter,
          word_d.get(),
          fedId_d.get(),
          digis_d.view().xx(),
          digis_d.view().yy(),
          digis_d.view().adc(),
          digis_d.view().pdigi(),
          digis_d.view().rawIdArr(),
          digis_d.view().moduleInd(),
          digiErrors_d.error(),  // returns nullptr if default-constructed
          useQualityInfo,
          includeErrors,
          debug);
      cudaCheck(cudaGetLastError());
#ifdef GPU_DEBUG
      cudaDeviceSynchronize();
      cudaCheck(cudaGetLastError());
#endif

      if (includeErrors) {
        digiErrors_d.copyErrorToHostAsync(stream);
      }
    }
    // End of Raw2Digi and passing data for clustering

    {
      // clusterizer ...
      using namespace gpuClustering;
      int threadsPerBlock = 256;
      int blocks = (std::max(int(wordCounter), int(phase1PixelTopology::numberOfModules)) + threadsPerBlock - 1) /
                   threadsPerBlock;

      if (isRun2)
        gpuCalibPixel::calibDigis<true><<<blocks, threadsPerBlock, 0, stream>>>(digis_d.view().moduleInd(),
                                                                                digis_d.view().xx(),
                                                                                digis_d.view().yy(),
                                                                                digis_d.view().adc(),
                                                                                gains,
                                                                                wordCounter,
                                                                                clusters_d.moduleStart(),
                                                                                clusters_d.clusInModule(),
                                                                                clusters_d.clusModuleStart());
      else
        gpuCalibPixel::calibDigis<false><<<blocks, threadsPerBlock, 0, stream>>>(digis_d.view().moduleInd(),
                                                                                 digis_d.view().xx(),
                                                                                 digis_d.view().yy(),
                                                                                 digis_d.view().adc(),
                                                                                 gains,
                                                                                 wordCounter,
                                                                                 clusters_d.moduleStart(),
                                                                                 clusters_d.clusInModule(),
                                                                                 clusters_d.clusModuleStart());

      cudaCheck(cudaGetLastError());
#ifdef GPU_DEBUG
      cudaDeviceSynchronize();
      cudaCheck(cudaGetLastError());
#endif

#ifdef GPU_DEBUG
      std::cout << "CUDA countModules kernel launch with " << blocks << " blocks of " << threadsPerBlock
                << " threads\n";
#endif

      countModules<false><<<blocks, threadsPerBlock, 0, stream>>>(
          digis_d.view().moduleInd(), clusters_d.moduleStart(), digis_d.view().clus(), wordCounter);
      cudaCheck(cudaGetLastError());

      threadsPerBlock = 256 + 128;  /// should be larger than 6000/16 aka (maxPixInModule/maxiter in the kernel)
      blocks = phase2PixelTopology::numberOfModules;
#ifdef GPU_DEBUG
      std::cout << "CUDA findClus kernel launch with " << blocks << " blocks of " << threadsPerBlock << " threads\n";
#endif

      findClus<false><<<blocks, threadsPerBlock, 0, stream>>>(digis_d.view().moduleInd(),
                                                              digis_d.view().xx(),
                                                              digis_d.view().yy(),
                                                              clusters_d.moduleStart(),
                                                              clusters_d.clusInModule(),
                                                              clusters_d.moduleId(),
                                                              digis_d.view().clus(),
                                                              wordCounter);
      cudaCheck(cudaGetLastError());
#ifdef GPU_DEBUG
      cudaDeviceSynchronize();
      cudaCheck(cudaGetLastError());
#endif

      // apply charge cut
      clusterChargeCut<false><<<blocks, threadsPerBlock, 0, stream>>>(clusterThresholds,
                                                                      digis_d.view().moduleInd(),
                                                                      digis_d.view().adc(),
                                                                      clusters_d.moduleStart(),
                                                                      clusters_d.clusInModule(),
                                                                      clusters_d.moduleId(),
                                                                      digis_d.view().clus(),
                                                                      wordCounter);

      cudaCheck(cudaGetLastError());

      // count the module start indices already here (instead of
      // rechits) so that the number of clusters/hits can be made
      // available in the rechit producer without additional points of
      // synchronization/ExternalWork
      auto nModules_Clusters_d = cms::cuda::make_device_unique<uint32_t[]>(3, stream);
      // MUST be ONE block
      fillHitsModuleStart<false><<<1, 1024, 0, stream>>>(
          clusters_d.clusInModule(), clusters_d.clusModuleStart(), clusters_d.moduleStart(), nModules_Clusters_d.get());

      // copy to host
      nModules_Clusters_h = cms::cuda::make_host_unique<uint32_t[]>(3, stream);
      cudaCheck(cudaMemcpyAsync(
          nModules_Clusters_h.get(), nModules_Clusters_d.get(), 3 * sizeof(uint32_t), cudaMemcpyDefault, stream));

#ifdef GPU_DEBUG
      cudaDeviceSynchronize();
      cudaCheck(cudaGetLastError());
#endif

    }  // end clusterizer scope
  }

  void SiPixelRawToClusterGPUKernel::makePhase2ClustersAsync(const SiPixelClusterThresholds clusterThresholds,
                                                             const uint16_t *moduleIds,
                                                             const uint16_t *xDigis,
                                                             const uint16_t *yDigis,
                                                             const uint16_t *adcDigis,
                                                             const uint32_t *packedData,
                                                             const uint32_t *rawIds,
                                                             const uint32_t numDigis,
                                                             cudaStream_t stream) {
    using namespace gpuClustering;
    nDigis = numDigis;
    digis_d = SiPixelDigisCUDA(numDigis, stream);

    cudaCheck(
        cudaMemcpyAsync(digis_d.view().moduleInd(), moduleIds, sizeof(uint16_t) * numDigis, cudaMemcpyDefault, stream));
    cudaCheck(cudaMemcpyAsync(digis_d.view().xx(), xDigis, sizeof(uint16_t) * numDigis, cudaMemcpyDefault, stream));
    cudaCheck(cudaMemcpyAsync(digis_d.view().yy(), yDigis, sizeof(uint16_t) * numDigis, cudaMemcpyDefault, stream));
    cudaCheck(cudaMemcpyAsync(digis_d.view().adc(), adcDigis, sizeof(uint16_t) * numDigis, cudaMemcpyDefault, stream));
    cudaCheck(
        cudaMemcpyAsync(digis_d.view().pdigi(), packedData, sizeof(uint32_t) * numDigis, cudaMemcpyDefault, stream));
    cudaCheck(
        cudaMemcpyAsync(digis_d.view().rawIdArr(), rawIds, sizeof(uint32_t) * numDigis, cudaMemcpyDefault, stream));

    clusters_d = SiPixelClustersCUDA(phase2PixelTopology::numberOfModules, stream);

    nModules_Clusters_h = cms::cuda::make_host_unique<uint32_t[]>(2, stream);

    int threadsPerBlock = 512;
    int blocks = (int(numDigis) + threadsPerBlock - 1) / threadsPerBlock;

    gpuCalibPixel::calibDigisPhase2<<<blocks, threadsPerBlock, 0, stream>>>(digis_d.view().moduleInd(),
                                                                            digis_d.view().adc(),
                                                                            numDigis,
                                                                            clusters_d.moduleStart(),
                                                                            clusters_d.clusInModule(),
                                                                            clusters_d.clusModuleStart());

    cudaCheck(cudaGetLastError());

#ifdef GPU_DEBUG
    cudaDeviceSynchronize();
    cudaCheck(cudaGetLastError());
#endif

#ifdef GPU_DEBUG
    std::cout << "CUDA countModules kernel launch with " << blocks << " blocks of " << threadsPerBlock << " threads\n";
#endif

    countModules<true><<<blocks, threadsPerBlock, 0, stream>>>(
        digis_d.view().moduleInd(), clusters_d.moduleStart(), digis_d.view().clus(), numDigis);
    cudaCheck(cudaGetLastError());

    // read the number of modules into a data member, used by getProduct())
    cudaCheck(cudaMemcpyAsync(
        &(nModules_Clusters_h[0]), clusters_d.moduleStart(), sizeof(uint32_t), cudaMemcpyDefault, stream));

    threadsPerBlock = 256;
    blocks = phase2PixelTopology::numberOfModules;

    findClus<true><<<blocks, threadsPerBlock, 0, stream>>>(digis_d.view().moduleInd(),
                                                           digis_d.view().xx(),
                                                           digis_d.view().yy(),
                                                           clusters_d.moduleStart(),
                                                           clusters_d.clusInModule(),
                                                           clusters_d.moduleId(),
                                                           digis_d.view().clus(),
                                                           numDigis);

    cudaCheck(cudaGetLastError());
#ifdef GPU_DEBUG
    cudaDeviceSynchronize();
    cudaCheck(cudaGetLastError());
#endif

    // apply charge cut
    clusterChargeCut<true><<<blocks, threadsPerBlock, 0, stream>>>(clusterThresholds,
                                                                   digis_d.view().moduleInd(),
                                                                   digis_d.view().adc(),
                                                                   clusters_d.moduleStart(),
                                                                   clusters_d.clusInModule(),
                                                                   clusters_d.moduleId(),
                                                                   digis_d.view().clus(),
                                                                   numDigis);
    cudaCheck(cudaGetLastError());

    auto nModules_Clusters_d = cms::cuda::make_device_unique<uint32_t[]>(3, stream);
    // MUST be ONE block
    fillHitsModuleStart<true><<<1, 1024, 0, stream>>>(
        clusters_d.clusInModule(), clusters_d.clusModuleStart(), clusters_d.moduleStart(), nModules_Clusters_d.get());

    nModules_Clusters_h = cms::cuda::make_host_unique<uint32_t[]>(3, stream);
    cudaCheck(cudaMemcpyAsync(
        nModules_Clusters_h.get(), nModules_Clusters_d.get(), 3 * sizeof(uint32_t), cudaMemcpyDefault, stream));

#ifdef GPU_DEBUG
    cudaDeviceSynchronize();
    cudaCheck(cudaGetLastError());
#endif
  }  //
}  // namespace pixelgpudetails
