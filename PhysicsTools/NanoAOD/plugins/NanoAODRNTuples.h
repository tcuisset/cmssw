#ifndef PhysicsTools_NanoAOD_NanoAODRNTuples_h
#define PhysicsTools_NanoAOD_NanoAODRNTuples_h

#include "DataFormats/Provenance/interface/BranchType.h"
#include "DataFormats/Provenance/interface/ParameterSetBlob.h"
#include "DataFormats/Provenance/interface/ParameterSetID.h"
#include "FWCore/Framework/interface/LuminosityBlockForOutput.h"
#include "FWCore/Framework/interface/RunForOutput.h"
#include "FWCore/ParameterSet/interface/Registry.h"

#include "TFile.h"
#include <ROOT/RNTuple.hxx>
using ROOT::Experimental::RNTupleWriter;

#include "RNTupleFieldPtr.h"

class LumiNTuple {
public:
  LumiNTuple() = default;
  void fill(const edm::LuminosityBlockID& id, TFile& file);
  void finalizeWrite();
private:
  void createFields(const edm::LuminosityBlockID& id, TFile& file);
  std::unique_ptr<RNTupleWriter> m_ntuple;
  RNTupleFieldPtr<UInt_t> m_run;
  RNTupleFieldPtr<UInt_t> m_luminosityBlock;
};

class RunNTuple {
public:
  RunNTuple() = default;
  void fill(const edm::RunForOutput& iRun, TFile& file);
  void finalizeWrite();
private:
  void createFields(const edm::RunForOutput& iRun, TFile& file);
  std::unique_ptr<RNTupleWriter> m_ntuple;
  RNTupleFieldPtr<UInt_t> m_run;
  // TODO SummaryTableOutput fields
};

class PSetNTuple {
public:
  PSetNTuple() = default;
  void fill(edm::pset::Registry* pset, TFile& file);
  void finalizeWrite();
private:
  using PSetType = std::pair<std::string, std::string>;
  //using PSetType = std::pair<edm::ParameterSetID, edm::ParameterSetBlob>;
  void createFields(TFile& file);
  std::unique_ptr<RNTupleWriter> m_ntuple;
  RNTupleFieldPtr<PSetType> m_pset;
};

#endif
