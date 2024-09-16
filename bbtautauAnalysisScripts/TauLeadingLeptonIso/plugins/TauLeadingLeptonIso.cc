// -*- C++ -*-
//
// Package:    bbtautauAnalysisScripts/TauLeadingLeptonIso
// Class:      TauLeadingLeptonIso
// 
/**\class TauLeadingLeptonIso TauLeadingLeptonIso.cc bbtautauAnalysisScripts/TauLeadingLeptonIso/plugins/TauLeadingLeptonIso.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrew Loeliger
//         Created:  Wed, 20 Jul 2022 18:32:01 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/InputTag.h"

//#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "CommonTools/Egamma/interface/EffectiveAreas.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"


struct leptonInfo
{
float pt = -1.0;
float eta = 0.0;
float phi = 0.0;
float m = 0.0;

float sumPFChargedHadronPt = 0.0;
float sumPFNeutralHadronPt = 0.0;
float sumPFPhotonPt = 0.0;

float tausumPFChargedHadronPt = 0.0;
float tausumPFNeutralHadronPt = 0.0;
float tausumPFPhotonPt= 0.0;

float sumPUPt = 0.0;

float rho = 0.0;
float ea = 0.0;

//Adding the addtional delR variable
float delR = -404.0;
//end

float correctedSumPFChargedHadronPt = 0.0;
float correctedSumPFNeutralHadronPt = 0.0;

float correctedIso = 0.0;
};

//
// class declaration
//

class TauLeadingLeptonIso : public edm::stream::EDProducer<> {
   public:
      explicit TauLeadingLeptonIso(const edm::ParameterSet&);
      ~TauLeadingLeptonIso();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
      virtual void calculateCorrectedMuonIsoInformation(pat::Tau theTau, leptonInfo &theMuonInfo);
      virtual void calculateCorrectedElectronIsoInformation(pat::Tau theTau, leptonInfo &theElectronInfo);

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
  edm::EDGetTokenT< std::vector<pat::Muon> > muonCollection;
  edm::EDGetTokenT< std::vector<pat::Electron> > electronCollection;
  edm::EDGetTokenT< std::vector<pat::Tau> > TauCollection;
  edm::EDGetTokenT<double> rhoSrc;
  EffectiveAreas theEffectiveAreas;
  bool verboseDebug;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
TauLeadingLeptonIso::TauLeadingLeptonIso(const edm::ParameterSet& iConfig):
  muonCollection(consumes< std::vector<pat::Muon> >(iConfig.getParameter< edm::InputTag >("muonCollection"))),
  electronCollection(consumes< std::vector<pat::Electron> >(iConfig.getParameter< edm::InputTag >("electronCollection"))),
  TauCollection(consumes< std::vector<pat::Tau> >(iConfig.getParameter< edm::InputTag >("TauCollection"))),
  rhoSrc(consumes<double > (iConfig.getParameter< edm::InputTag >("rhoSrc"))),
  theEffectiveAreas(iConfig.getParameter< edm::FileInPath>("EAConfigFile").fullPath())
{
  verboseDebug = iConfig.exists("verboseDebug") ? iConfig.getParameter<bool>("verboseDebug"): false; 
  //Muon isolation products
  //Stores pt eta phi m
  //and implied isolation of leading, sub-leading, and sub-subleading muons
  //if this tau is used

  //Adding the addtional DelR branches for muons
  produces<edm::ValueMap<float>>("LeadingMuondelR");
  produces<edm::ValueMap<float>>("SubLeadingMuondelR");
  produces<edm::ValueMap<float>>("SubSubLeadingMuondelR");
  //end

  produces<edm::ValueMap<float>>("LeadingMuonPt");
  produces<edm::ValueMap<float>>("LeadingMuonEta");
  produces<edm::ValueMap<float>>("LeadingMuonPhi");
  produces<edm::ValueMap<float>>("LeadingMuonM");
  produces<edm::ValueMap<float>>("LeadingMuonCorrIso");
  produces<edm::ValueMap<float>>("LeadingMuonsumPFChargedHadronPt");
  produces<edm::ValueMap<float>>("LeadingMuonsumPFNeutralHadronPt");
  produces<edm::ValueMap<float>>("LeadingMuonsumPFPhotonPt");
  produces<edm::ValueMap<float>>("LeadingMuonsumPUPt");
  produces<edm::ValueMap<float>>("LeadingMuontausumPFChargedHadronPt");
  produces<edm::ValueMap<float>>("LeadingMuontausumPFNeutralHadronPt");
  produces<edm::ValueMap<float>>("LeadingMuontausumPFPhotonPt");
  

  produces<edm::ValueMap<float>>("SubLeadingMuonPt");
  produces<edm::ValueMap<float>>("SubLeadingMuonEta");
  produces<edm::ValueMap<float>>("SubLeadingMuonPhi");
  produces<edm::ValueMap<float>>("SubLeadingMuonM");
  produces<edm::ValueMap<float>>("SubLeadingMuonCorrIso");
  produces<edm::ValueMap<float>>("SubLeadingMuonsumPFChargedHadronPt");
  produces<edm::ValueMap<float>>("SubLeadingMuonsumPFNeutralHadronPt");
  produces<edm::ValueMap<float>>("SubLeadingMuonsumPFPhotonPt");
  produces<edm::ValueMap<float>>("SubLeadingMuonsumPUPt");
  produces<edm::ValueMap<float>>("SubLeadingMuontausumPFChargedHadronPt");
  produces<edm::ValueMap<float>>("SubLeadingMuontausumPFNeutralHadronPt");
  produces<edm::ValueMap<float>>("SubLeadingMuontausumPFPhotonPt");

  produces<edm::ValueMap<float>>("SubSubLeadingMuonPt");
  produces<edm::ValueMap<float>>("SubSubLeadingMuonEta");
  produces<edm::ValueMap<float>>("SubSubLeadingMuonPhi");
  produces<edm::ValueMap<float>>("SubSubLeadingMuonM");
  produces<edm::ValueMap<float>>("SubSubLeadingMuonCorrIso");
  produces<edm::ValueMap<float>>("SubSubLeadingMuonsumPFChargedHadronPt");
  produces<edm::ValueMap<float>>("SubSubLeadingMuonsumPFNeutralHadronPt");
  produces<edm::ValueMap<float>>("SubSubLeadingMuonsumPFPhotonPt");
  produces<edm::ValueMap<float>>("SubSubLeadingMuonsumPUPt");
  produces<edm::ValueMap<float>>("SubSubLeadingMuontausumPFChargedHadronPt");
  produces<edm::ValueMap<float>>("SubSubLeadingMuontausumPFNeutralHadronPt");
  produces<edm::ValueMap<float>>("SubSubLeadingMuontausumPFPhotonPt");

  produces<edm::ValueMap<int>>("Mcounter");

  //Electron isolation products
  //Stores pt eta phi m
  //and implied isolation of leading, sub-leading, and sub-subleading electons
  //if this tau is used

  //Adding the addtional DelR branches for electrons
  produces<edm::ValueMap<float>>("LeadingElectrondelR"); 
  produces<edm::ValueMap<float>>("SubLeadingElectrondelR");
  produces<edm::ValueMap<float>>("SubSubLeadingElectrondelR");
  //end 


  produces<edm::ValueMap<float>>("LeadingElectronPt");
  produces<edm::ValueMap<float>>("LeadingElectronEta");
  produces<edm::ValueMap<float>>("LeadingElectronPhi");
  produces<edm::ValueMap<float>>("LeadingElectronM");
  produces<edm::ValueMap<float>>("LeadingElectronCorrIso");
  produces<edm::ValueMap<float>>("LeadingElectronsumPFChargedHadronPt");
  produces<edm::ValueMap<float>>("LeadingElectronsumPFNeutralHadronPt");
  produces<edm::ValueMap<float>>("LeadingElectronsumPFPhotonPt");
  produces<edm::ValueMap<float>>("LeadingElectronea");
  produces<edm::ValueMap<float>>("LeadingElectronrho");
  produces<edm::ValueMap<float>>("LeadingElectrontausumPFChargedHadronPt");
  produces<edm::ValueMap<float>>("LeadingElectrontausumPFNeutralHadronPt");
  produces<edm::ValueMap<float>>("LeadingElectrontausumPFPhotonPt");


  produces<edm::ValueMap<float>>("SubLeadingElectronPt");
  produces<edm::ValueMap<float>>("SubLeadingElectronEta");
  produces<edm::ValueMap<float>>("SubLeadingElectronPhi");
  produces<edm::ValueMap<float>>("SubLeadingElectronM");
  produces<edm::ValueMap<float>>("SubLeadingElectronCorrIso");
  produces<edm::ValueMap<float>>("SubLeadingElectronsumPFChargedHadronPt");
  produces<edm::ValueMap<float>>("SubLeadingElectronsumPFNeutralHadronPt");
  produces<edm::ValueMap<float>>("SubLeadingElectronsumPFPhotonPt");
  produces<edm::ValueMap<float>>("SubLeadingElectronea");
  produces<edm::ValueMap<float>>("SubLeadingElectronrho");
  produces<edm::ValueMap<float>>("SubLeadingElectrontausumPFChargedHadronPt");
  produces<edm::ValueMap<float>>("SubLeadingElectrontausumPFNeutralHadronPt");
  produces<edm::ValueMap<float>>("SubLeadingElectrontausumPFPhotonPt");

  produces<edm::ValueMap<float>>("SubSubLeadingElectronPt");
  produces<edm::ValueMap<float>>("SubSubLeadingElectronEta");
  produces<edm::ValueMap<float>>("SubSubLeadingElectronPhi");
  produces<edm::ValueMap<float>>("SubSubLeadingElectronM");
  produces<edm::ValueMap<float>>("SubSubLeadingElectronCorrIso");
  produces<edm::ValueMap<float>>("SubSubLeadingElectronsumPFChargedHadronPt");
  produces<edm::ValueMap<float>>("SubSubLeadingElectronsumPFNeutralHadronPt");
  produces<edm::ValueMap<float>>("SubSubLeadingElectronsumPFPhotonPt");
  produces<edm::ValueMap<float>>("SubSubLeadingElectronea");
  produces<edm::ValueMap<float>>("SubSubLeadingElectronrho");
  produces<edm::ValueMap<float>>("SubSubLeadingElectrontausumPFChargedHadronPt");
  produces<edm::ValueMap<float>>("SubSubLeadingElectrontausumPFNeutralHadronPt");
  produces<edm::ValueMap<float>>("SubSubLeadingElectrontausumPFPhotonPt");


  produces<edm::ValueMap<int>>("Ecounter");
}


TauLeadingLeptonIso::~TauLeadingLeptonIso()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
TauLeadingLeptonIso::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   edm::Handle< std::vector<pat::Muon> > muonHandle;
   iEvent.getByToken(muonCollection, muonHandle);

   edm::Handle< std::vector<pat::Electron> > electronHandle;
   iEvent.getByToken(electronCollection, electronHandle);
   
   edm::Handle< std::vector<pat::Tau> > TauHandle;
   iEvent.getByToken(TauCollection, TauHandle);

   edm::Handle<double> rho;
   iEvent.getByToken(rhoSrc, rho);
   
   int nMuons = muonHandle->size();
   int nElectrons = electronHandle->size();
   if (verboseDebug) std::cout<<"nMuons: "<<nMuons<<std::endl;
   if (verboseDebug) std::cout<<"nElectrons: "<<nElectrons<<std::endl;

   double dRcutoff_ele = 0.6; //0.3 reco cone for tau and 0.3 for electron
   double dRcutoff_muon = 0.7; //0.3 reco cone for tau and 0.4 for muon
   double deltaR = 0.0;

   int eleCounter = 0;
   int mCounter = 0;


   //These vectors will store the information about corrected leptons that we make later
   std::vector<float> leadingMuonVector_delR, leadingMuonVector_pt, leadingMuonVector_eta, leadingMuonVector_phi, leadingMuonVector_m, leadingMuonVector_corrIso, leadingMuonVector_sumPFChargedHadronPt, leadingMuonVector_sumPFNeutralHadronPt, leadingMuonVector_sumPFPhotonPt, leadingMuonVector_sumPUPt, leadingMuonVector_tausumPFChargedHadronPt, leadingMuonVector_tausumPFNeutralHadronPt, leadingMuonVector_tausumPFPhotonPt;
   std::vector<float> subleadingMuonVector_delR, subleadingMuonVector_pt, subleadingMuonVector_eta, subleadingMuonVector_phi, subleadingMuonVector_m, subleadingMuonVector_corrIso, subleadingMuonVector_sumPFChargedHadronPt, subleadingMuonVector_sumPFNeutralHadronPt, subleadingMuonVector_sumPFPhotonPt, subleadingMuonVector_sumPUPt, subleadingMuonVector_tausumPFChargedHadronPt, subleadingMuonVector_tausumPFNeutralHadronPt, subleadingMuonVector_tausumPFPhotonPt;
   std::vector<float> subsubleadingMuonVector_delR, subsubleadingMuonVector_pt, subsubleadingMuonVector_eta, subsubleadingMuonVector_phi, subsubleadingMuonVector_m, subsubleadingMuonVector_corrIso, subsubleadingMuonVector_sumPFChargedHadronPt, subsubleadingMuonVector_sumPFNeutralHadronPt, subsubleadingMuonVector_sumPFPhotonPt, subsubleadingMuonVector_sumPUPt, subsubleadingMuonVector_tausumPFChargedHadronPt, subsubleadingMuonVector_tausumPFNeutralHadronPt, subsubleadingMuonVector_tausumPFPhotonPt;
   std::vector<int> Mcounter;

   std::vector<float> leadingElectronVector_delR, leadingElectronVector_pt, leadingElectronVector_eta, leadingElectronVector_phi, leadingElectronVector_m, leadingElectronVector_corrIso, leadingElectronVector_sumPFChargedHadronPt, leadingElectronVector_sumPFNeutralHadronPt, leadingElectronVector_sumPFPhotonPt, leadingElectronVector_ea, leadingElectronVector_rho, leadingElectronVector_tausumPFChargedHadronPt, leadingElectronVector_tausumPFNeutralHadronPt, leadingElectronVector_tausumPFPhotonPt;
   std::vector<float> subleadingElectronVector_delR, subleadingElectronVector_pt, subleadingElectronVector_eta, subleadingElectronVector_phi, subleadingElectronVector_m, subleadingElectronVector_corrIso, subleadingElectronVector_sumPFChargedHadronPt, subleadingElectronVector_sumPFNeutralHadronPt, subleadingElectronVector_sumPFPhotonPt, subleadingElectronVector_ea, subleadingElectronVector_rho, subleadingElectronVector_tausumPFChargedHadronPt, subleadingElectronVector_tausumPFNeutralHadronPt, subleadingElectronVector_tausumPFPhotonPt;
   std::vector<float> subsubleadingElectronVector_delR, subsubleadingElectronVector_pt, subsubleadingElectronVector_eta, subsubleadingElectronVector_phi, subsubleadingElectronVector_m, subsubleadingElectronVector_corrIso, subsubleadingElectronVector_sumPFChargedHadronPt, subsubleadingElectronVector_sumPFNeutralHadronPt, subsubleadingElectronVector_sumPFPhotonPt, subsubleadingElectronVector_ea, subsubleadingElectronVector_rho, subsubleadingElectronVector_tausumPFChargedHadronPt, subsubleadingElectronVector_tausumPFNeutralHadronPt, subsubleadingElectronVector_tausumPFPhotonPt;
   std::vector<int> Ecounter;

   //Okay, the idea here is that for each tau we have,
   //we go through and check each lepton
   //we first make a structure that contains it's 4 vector info,
   //and then the corrected isolation values that would be implied from that tau
   //we store this for the leading, subleading, and sub-sub-leading electrons, and muons

   leptonInfo nullInfo;

   for(std::vector<pat::Tau>::const_iterator theTau = TauHandle->begin();
       theTau != TauHandle->end();
       ++theTau)
     {
       eleCounter = 0;
       mCounter = 0;

       std::vector< leptonInfo > electronInformation;
       std::vector< leptonInfo > muonInformation;

       //For each tau we now go through the list of muons and electrons, and select 
       // the three leading candidates of each, and store their information in our vectors of
       // lepton information

       //once we have their information, we ca go through and calculate the corrected isolation
       //of each of the supplied leptons, with respect to the given tau
       for(std::vector<pat::Muon>::const_iterator theMuon = muonHandle->begin();
      theMuon != muonHandle->end();
      ++theMuon)
    {
      leptonInfo currentMuonInfo;
      currentMuonInfo.pt = theMuon->pt();
      currentMuonInfo.eta = theMuon->eta();
      currentMuonInfo.phi = theMuon->phi();
      currentMuonInfo.m = theMuon->mass();
      currentMuonInfo.sumPFChargedHadronPt = theMuon->pfIsolationR04().sumChargedHadronPt;
      currentMuonInfo.sumPFNeutralHadronPt = theMuon->pfIsolationR04().sumNeutralHadronEt;
      currentMuonInfo.sumPFPhotonPt = theMuon->pfIsolationR04().sumPhotonEt;
      currentMuonInfo.sumPUPt = theMuon->pfIsolationR04().sumPUPt;

      //loop through our information collection
      //If we have higher pt than the current entry, we insert this lepton's information before

     deltaR = reco::deltaR(currentMuonInfo.eta, currentMuonInfo.phi, theTau->eta(), theTau->phi());

     //add the delR variable to lepton info
     currentMuonInfo.delR = deltaR;
     //end     

     //if (deltaR < dRmin_muon && deltaR > 0.05 && theMuon->passed(reco::Muon::CutBasedIdLoose))
     if ((deltaR < dRcutoff_muon) && (deltaR > 0.05) && theMuon->passed(reco::Muon::CutBasedIdLoose) && (theMuon->pt()>15))
     {
        bool insertAtEnd = true;
        mCounter++;
        for (std::vector< leptonInfo >::const_iterator muonInfoIt = muonInformation.begin();
      muonInfoIt != muonInformation.end();
      ++muonInfoIt)
        {
          if (currentMuonInfo.pt > (*muonInfoIt).pt) 
       {
         insertAtEnd = false;
         muonInformation.insert(muonInfoIt, currentMuonInfo);
         break;
       }
        }
      //If we are at the end, insert the information at the end
      if (insertAtEnd) muonInformation.insert(muonInformation.end(), currentMuonInfo);
      //then if we have more than 3 entries in the list of information, get rid of the 
      //last entry
      if (muonInformation.size() > 3) muonInformation.pop_back();

     }
      
    }
       //Now that we have all of the leading muons and their information, let's go through
       //and calculated a rectified muon isolation for each of them
       for (std::vector< leptonInfo >::iterator muonInfoIt = muonInformation.begin();
       muonInfoIt != muonInformation.end();
       ++muonInfoIt) this->calculateCorrectedMuonIsoInformation(*theTau, *muonInfoIt);
       //if any slots in the information vector are empty, let's create null information
       //to fill them
       int nullMuonEntriesNeeded = (int)(3-muonInformation.size());
       for (int i=0; i< nullMuonEntriesNeeded;++i)
    {
      //leptonInfo nullInfo;
      muonInformation.push_back(nullInfo);
    }

       //Now we do something similar for the electrons
       for(std::vector<pat::Electron>::const_iterator theElectron = electronHandle->begin();
      theElectron != electronHandle->end();
      ++theElectron)
    {
      leptonInfo currentElectronInfo;
      currentElectronInfo.pt = theElectron->pt();
      currentElectronInfo.eta = theElectron->eta();
      currentElectronInfo.phi = theElectron->phi();
      currentElectronInfo.m = theElectron->mass();
      currentElectronInfo.sumPFChargedHadronPt = theElectron->pfIsolationVariables().sumChargedHadronPt;
      currentElectronInfo.sumPFNeutralHadronPt = theElectron->pfIsolationVariables().sumNeutralHadronEt;
      currentElectronInfo.sumPFPhotonPt = theElectron->pfIsolationVariables().sumPhotonEt;
      
      currentElectronInfo.rho = *rho;
      currentElectronInfo.ea = theEffectiveAreas.getEffectiveArea(fabs(theElectron->superCluster()->eta()));

      //loop through our information collection
      //if we have a higher pt than the current entry, we insert this lepton's information before
     deltaR = reco::deltaR(currentElectronInfo.eta, currentElectronInfo.phi, theTau->eta(), theTau->phi());

     //add the deltaR variale to lepton info
     currentElectronInfo.delR = deltaR;
     //end

     //if (deltaR < dRmin && deltaR > 0.05 && theElectron->electronID("cutBasedElectronID-Fall17-94X-V2-loose"))
     //std::cout<< "Debug Electron ID = "<<((theElectron->userInt("cutBasedElectronID-Fall17-94X-V2-loose") & 0x37F ) == 0x37F)<<" fullID = "<<theElectron->electronID("cutBasedElectronID-Fall17-94X-V2-loose")<<"\n";
     //if (deltaR < dRmin_ele && deltaR > 0.05 && ((theElectron->userInt("cutBasedElectronID-Fall17-94X-V2-loose") & 0x37F ) == 0x37F))
     if ((deltaR < dRcutoff_ele) && (deltaR > 0.05) && ((theElectron->userInt("cutBasedElectronID-Fall17-94X-V2-loose") & 0x37F ) == 0x37F) && (theElectron->pt()>5))
     {
     eleCounter++;
      bool insertAtEnd = true;
      for(std::vector< leptonInfo >::const_iterator electronInfoIt = electronInformation.begin();
          electronInfoIt != electronInformation.end();
          ++electronInfoIt)
        {
          if(currentElectronInfo.pt > (*electronInfoIt).pt)
       {
         insertAtEnd = false;
         electronInformation.insert(electronInfoIt, currentElectronInfo);
         break;
       }
        }
      //if we are at the end, insert the information at the end
      if (insertAtEnd) electronInformation.insert(electronInformation.end(), currentElectronInfo);
      //now, if we have more than 3 entries, get rid of the last entry in the list
      if (electronInformation.size() > 3) electronInformation.pop_back();
     }
    }
       //Now that we have all of the leading electrons and their information, let's go through
       //and calculated a rectified electron isolation for each of them
       for (std::vector< leptonInfo >::iterator electronInfoIt = electronInformation.begin();
       electronInfoIt != electronInformation.end();
       ++electronInfoIt) this->calculateCorrectedMuonIsoInformation(*theTau, *electronInfoIt);
       //if any slots in the information vector are empty, let's create null information
       //to fill them
       int nullElectronEntriesNeeded = (int)(3-electronInformation.size());
       for (int i=0; i<nullElectronEntriesNeeded;++i)
    {
      leptonInfo nullInfo;
      electronInformation.push_back(nullInfo);
    }
       //Now that we have all the correct information for this tau, we can read out all the information
       //to a series of vectors that we will store later
      //Add delR info for muons
       leadingMuonVector_delR.push_back(muonInformation[0].delR);
       subleadingMuonVector_delR.push_back(muonInformation[1].delR);
       subsubleadingMuonVector_delR.push_back(muonInformation[2].delR);
       //end

       leadingMuonVector_pt.push_back(muonInformation[0].pt); 
       leadingMuonVector_eta.push_back(muonInformation[0].eta); 
       leadingMuonVector_phi.push_back(muonInformation[0].phi);
       leadingMuonVector_m.push_back(muonInformation[0].m);
       leadingMuonVector_corrIso.push_back(muonInformation[0].correctedIso);
       leadingMuonVector_sumPFChargedHadronPt.push_back(muonInformation[0].sumPFChargedHadronPt);
       leadingMuonVector_sumPFNeutralHadronPt.push_back(muonInformation[0].sumPFNeutralHadronPt);
       leadingMuonVector_sumPFPhotonPt.push_back(muonInformation[0].sumPFPhotonPt);
       leadingMuonVector_sumPUPt.push_back(muonInformation[0].sumPUPt);
       leadingMuonVector_tausumPFChargedHadronPt.push_back(muonInformation[0].tausumPFChargedHadronPt);
       leadingMuonVector_tausumPFNeutralHadronPt.push_back(muonInformation[0].tausumPFNeutralHadronPt);
       leadingMuonVector_tausumPFPhotonPt.push_back(muonInformation[0].tausumPFPhotonPt);

       subleadingMuonVector_pt.push_back(muonInformation[1].pt); 
       subleadingMuonVector_eta.push_back(muonInformation[1].eta); 
       subleadingMuonVector_phi.push_back(muonInformation[1].phi);
       subleadingMuonVector_m.push_back(muonInformation[1].m);
       subleadingMuonVector_corrIso.push_back(muonInformation[1].correctedIso);
       subleadingMuonVector_sumPFChargedHadronPt.push_back(muonInformation[1].sumPFChargedHadronPt);
       subleadingMuonVector_sumPFNeutralHadronPt.push_back(muonInformation[1].sumPFNeutralHadronPt);
       subleadingMuonVector_sumPFPhotonPt.push_back(muonInformation[1].sumPFPhotonPt);
       subleadingMuonVector_sumPUPt.push_back(muonInformation[1].sumPUPt);
       subleadingMuonVector_tausumPFChargedHadronPt.push_back(muonInformation[1].tausumPFChargedHadronPt);
       subleadingMuonVector_tausumPFNeutralHadronPt.push_back(muonInformation[1].tausumPFNeutralHadronPt);
       subleadingMuonVector_tausumPFPhotonPt.push_back(muonInformation[1].tausumPFPhotonPt);


       subsubleadingMuonVector_pt.push_back(muonInformation[2].pt); 
       subsubleadingMuonVector_eta.push_back(muonInformation[2].eta); 
       subsubleadingMuonVector_phi.push_back(muonInformation[2].phi);
       subsubleadingMuonVector_m.push_back(muonInformation[2].m);
       subsubleadingMuonVector_corrIso.push_back(muonInformation[2].correctedIso);
       subsubleadingMuonVector_sumPFChargedHadronPt.push_back(muonInformation[2].sumPFChargedHadronPt);
       subsubleadingMuonVector_sumPFNeutralHadronPt.push_back(muonInformation[2].sumPFNeutralHadronPt);
       subsubleadingMuonVector_sumPFPhotonPt.push_back(muonInformation[2].sumPFPhotonPt);
       subsubleadingMuonVector_sumPUPt.push_back(muonInformation[2].sumPUPt);
       subsubleadingMuonVector_tausumPFChargedHadronPt.push_back(muonInformation[2].tausumPFChargedHadronPt);
       subsubleadingMuonVector_tausumPFNeutralHadronPt.push_back(muonInformation[2].tausumPFNeutralHadronPt);
       subsubleadingMuonVector_tausumPFPhotonPt.push_back(muonInformation[2].tausumPFPhotonPt);       


       Mcounter.push_back(mCounter);

      //Add delR info for electrons
       leadingElectronVector_delR.push_back(electronInformation[0].delR);
       subleadingElectronVector_delR.push_back(electronInformation[1].delR);
       subsubleadingElectronVector_delR.push_back(electronInformation[2].delR);
       //end


       leadingElectronVector_pt.push_back(electronInformation[0].pt); 
       leadingElectronVector_eta.push_back(electronInformation[0].eta); 
       leadingElectronVector_phi.push_back(electronInformation[0].phi);
       leadingElectronVector_m.push_back(electronInformation[0].m);
       leadingElectronVector_corrIso.push_back(electronInformation[0].correctedIso);
       leadingElectronVector_sumPFChargedHadronPt.push_back(electronInformation[0].sumPFChargedHadronPt);
       leadingElectronVector_sumPFNeutralHadronPt.push_back(electronInformation[0].sumPFNeutralHadronPt);       
       leadingElectronVector_sumPFPhotonPt.push_back(electronInformation[0].sumPFPhotonPt);
       leadingElectronVector_tausumPFChargedHadronPt.push_back(electronInformation[0].tausumPFChargedHadronPt);       
       leadingElectronVector_tausumPFNeutralHadronPt.push_back(electronInformation[0].tausumPFNeutralHadronPt);
       leadingElectronVector_tausumPFPhotonPt.push_back(electronInformation[0].tausumPFPhotonPt);
       leadingElectronVector_rho.push_back(electronInformation[0].rho);
       leadingElectronVector_ea.push_back(electronInformation[0].ea);       

       subleadingElectronVector_pt.push_back(electronInformation[1].pt); 
       subleadingElectronVector_eta.push_back(electronInformation[1].eta); 
       subleadingElectronVector_phi.push_back(electronInformation[1].phi);
       subleadingElectronVector_m.push_back(electronInformation[1].m);
       subleadingElectronVector_corrIso.push_back(electronInformation[1].correctedIso);
       subleadingElectronVector_sumPFChargedHadronPt.push_back(electronInformation[1].sumPFChargedHadronPt);
       subleadingElectronVector_sumPFNeutralHadronPt.push_back(electronInformation[1].sumPFNeutralHadronPt);       
       subleadingElectronVector_sumPFPhotonPt.push_back(electronInformation[1].sumPFPhotonPt);
       subleadingElectronVector_tausumPFChargedHadronPt.push_back(electronInformation[1].tausumPFChargedHadronPt);       
       subleadingElectronVector_tausumPFNeutralHadronPt.push_back(electronInformation[1].tausumPFNeutralHadronPt);
       subleadingElectronVector_tausumPFPhotonPt.push_back(electronInformation[1].tausumPFPhotonPt);
       subleadingElectronVector_rho.push_back(electronInformation[1].rho);
       subleadingElectronVector_ea.push_back(electronInformation[1].ea);        


       subsubleadingElectronVector_pt.push_back(electronInformation[2].pt); 
       subsubleadingElectronVector_eta.push_back(electronInformation[2].eta); 
       subsubleadingElectronVector_phi.push_back(electronInformation[2].phi);
       subsubleadingElectronVector_m.push_back(electronInformation[2].m);
       subsubleadingElectronVector_corrIso.push_back(electronInformation[2].correctedIso);
       subsubleadingElectronVector_sumPFChargedHadronPt.push_back(electronInformation[2].sumPFChargedHadronPt);
       subsubleadingElectronVector_sumPFNeutralHadronPt.push_back(electronInformation[2].sumPFNeutralHadronPt);       
       subsubleadingElectronVector_sumPFPhotonPt.push_back(electronInformation[2].sumPFPhotonPt);
       subsubleadingElectronVector_tausumPFChargedHadronPt.push_back(electronInformation[2].tausumPFChargedHadronPt);       
       subsubleadingElectronVector_tausumPFNeutralHadronPt.push_back(electronInformation[2].tausumPFNeutralHadronPt);
       subsubleadingElectronVector_tausumPFPhotonPt.push_back(electronInformation[2].tausumPFPhotonPt);
       subsubleadingElectronVector_rho.push_back(electronInformation[2].rho);
       subsubleadingElectronVector_ea.push_back(electronInformation[2].ea); 
       
       Ecounter.push_back(eleCounter);
       
     }

   //we have all of the information for the taus in this event. We read this out to the 
   //edm format, and we're done.
   //New valueMaps added

   //Adding MUONS value Maps - Stuff that was not originally inlcuded by Andrew

   std::unique_ptr< edm::ValueMap < float > > leadingMuonVector_sumPFChargedHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingMuonVector_sumPFChargedHadronPt_valueMap(*leadingMuonVector_sumPFChargedHadronPt_valueMap);
   filler_leadingMuonVector_sumPFChargedHadronPt_valueMap.insert(TauHandle, leadingMuonVector_sumPFChargedHadronPt.begin(), leadingMuonVector_sumPFChargedHadronPt.end());
   filler_leadingMuonVector_sumPFChargedHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingMuonVector_sumPFNeutralHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingMuonVector_sumPFNeutralHadronPt_valueMap(*leadingMuonVector_sumPFNeutralHadronPt_valueMap);
   filler_leadingMuonVector_sumPFNeutralHadronPt_valueMap.insert(TauHandle, leadingMuonVector_sumPFNeutralHadronPt.begin(), leadingMuonVector_sumPFNeutralHadronPt.end());
   filler_leadingMuonVector_sumPFNeutralHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingMuonVector_sumPFPhotonPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingMuonVector_sumPFPhotonPt_valueMap(*leadingMuonVector_sumPFPhotonPt_valueMap);
   filler_leadingMuonVector_sumPFPhotonPt_valueMap.insert(TauHandle, leadingMuonVector_sumPFPhotonPt.begin(), leadingMuonVector_sumPFPhotonPt.end());
   filler_leadingMuonVector_sumPFPhotonPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingMuonVector_sumPUPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingMuonVector_sumPUPt_valueMap(*leadingMuonVector_sumPUPt_valueMap);
   filler_leadingMuonVector_sumPUPt_valueMap.insert(TauHandle, leadingMuonVector_sumPUPt.begin(), leadingMuonVector_sumPUPt.end());
   filler_leadingMuonVector_sumPUPt_valueMap.fill(); 

   std::unique_ptr< edm::ValueMap < float > > leadingMuonVector_tausumPFChargedHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingMuonVector_tausumPFChargedHadronPt_valueMap(*leadingMuonVector_tausumPFChargedHadronPt_valueMap);
   filler_leadingMuonVector_tausumPFChargedHadronPt_valueMap.insert(TauHandle, leadingMuonVector_tausumPFChargedHadronPt.begin(), leadingMuonVector_tausumPFChargedHadronPt.end());
   filler_leadingMuonVector_tausumPFChargedHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingMuonVector_tausumPFNeutralHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingMuonVector_tausumPFNeutralHadronPt_valueMap(*leadingMuonVector_tausumPFNeutralHadronPt_valueMap);
   filler_leadingMuonVector_tausumPFNeutralHadronPt_valueMap.insert(TauHandle, leadingMuonVector_tausumPFNeutralHadronPt.begin(), leadingMuonVector_tausumPFNeutralHadronPt.end());
   filler_leadingMuonVector_tausumPFNeutralHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingMuonVector_tausumPFPhotonPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingMuonVector_tausumPFPhotonPt_valueMap(*leadingMuonVector_tausumPFPhotonPt_valueMap);
   filler_leadingMuonVector_tausumPFPhotonPt_valueMap.insert(TauHandle, leadingMuonVector_tausumPFPhotonPt.begin(), leadingMuonVector_tausumPFPhotonPt.end());
   filler_leadingMuonVector_tausumPFPhotonPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingMuonVector_sumPFChargedHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingMuonVector_sumPFChargedHadronPt_valueMap(*subleadingMuonVector_sumPFChargedHadronPt_valueMap);
   filler_subleadingMuonVector_sumPFChargedHadronPt_valueMap.insert(TauHandle, subleadingMuonVector_sumPFChargedHadronPt.begin(), subleadingMuonVector_sumPFChargedHadronPt.end());
   filler_subleadingMuonVector_sumPFChargedHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingMuonVector_sumPFNeutralHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingMuonVector_sumPFNeutralHadronPt_valueMap(*subleadingMuonVector_sumPFNeutralHadronPt_valueMap);
   filler_subleadingMuonVector_sumPFNeutralHadronPt_valueMap.insert(TauHandle, subleadingMuonVector_sumPFNeutralHadronPt.begin(), subleadingMuonVector_sumPFNeutralHadronPt.end());
   filler_subleadingMuonVector_sumPFNeutralHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingMuonVector_sumPFPhotonPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingMuonVector_sumPFPhotonPt_valueMap(*subleadingMuonVector_sumPFPhotonPt_valueMap);
   filler_subleadingMuonVector_sumPFPhotonPt_valueMap.insert(TauHandle, subleadingMuonVector_sumPFPhotonPt.begin(), subleadingMuonVector_sumPFPhotonPt.end());
   filler_subleadingMuonVector_sumPFPhotonPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingMuonVector_sumPUPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingMuonVector_sumPUPt_valueMap(*subleadingMuonVector_sumPUPt_valueMap);
   filler_subleadingMuonVector_sumPUPt_valueMap.insert(TauHandle, subleadingMuonVector_sumPUPt.begin(), subleadingMuonVector_sumPUPt.end());
   filler_subleadingMuonVector_sumPUPt_valueMap.fill(); 

   std::unique_ptr< edm::ValueMap < float > > subleadingMuonVector_tausumPFChargedHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingMuonVector_tausumPFChargedHadronPt_valueMap(*subleadingMuonVector_tausumPFChargedHadronPt_valueMap);
   filler_subleadingMuonVector_tausumPFChargedHadronPt_valueMap.insert(TauHandle, subleadingMuonVector_tausumPFChargedHadronPt.begin(), subleadingMuonVector_tausumPFChargedHadronPt.end());
   filler_subleadingMuonVector_tausumPFChargedHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingMuonVector_tausumPFNeutralHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingMuonVector_tausumPFNeutralHadronPt_valueMap(*subleadingMuonVector_tausumPFNeutralHadronPt_valueMap);
   filler_subleadingMuonVector_tausumPFNeutralHadronPt_valueMap.insert(TauHandle, subleadingMuonVector_tausumPFNeutralHadronPt.begin(), subleadingMuonVector_tausumPFNeutralHadronPt.end());
   filler_subleadingMuonVector_tausumPFNeutralHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingMuonVector_tausumPFPhotonPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingMuonVector_tausumPFPhotonPt_valueMap(*subleadingMuonVector_tausumPFPhotonPt_valueMap);
   filler_subleadingMuonVector_tausumPFPhotonPt_valueMap.insert(TauHandle, subleadingMuonVector_tausumPFPhotonPt.begin(), subleadingMuonVector_tausumPFPhotonPt.end());
   filler_subleadingMuonVector_tausumPFPhotonPt_valueMap.fill();     


   std::unique_ptr< edm::ValueMap < float > > subsubleadingMuonVector_sumPFChargedHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingMuonVector_sumPFChargedHadronPt_valueMap(*subsubleadingMuonVector_sumPFChargedHadronPt_valueMap);
   filler_subsubleadingMuonVector_sumPFChargedHadronPt_valueMap.insert(TauHandle, subsubleadingMuonVector_sumPFChargedHadronPt.begin(), subsubleadingMuonVector_sumPFChargedHadronPt.end());
   filler_subsubleadingMuonVector_sumPFChargedHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingMuonVector_sumPFNeutralHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingMuonVector_sumPFNeutralHadronPt_valueMap(*subsubleadingMuonVector_sumPFNeutralHadronPt_valueMap);
   filler_subsubleadingMuonVector_sumPFNeutralHadronPt_valueMap.insert(TauHandle, subsubleadingMuonVector_sumPFNeutralHadronPt.begin(), subsubleadingMuonVector_sumPFNeutralHadronPt.end());
   filler_subsubleadingMuonVector_sumPFNeutralHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingMuonVector_sumPFPhotonPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingMuonVector_sumPFPhotonPt_valueMap(*subsubleadingMuonVector_sumPFPhotonPt_valueMap);
   filler_subsubleadingMuonVector_sumPFPhotonPt_valueMap.insert(TauHandle, subsubleadingMuonVector_sumPFPhotonPt.begin(), subsubleadingMuonVector_sumPFPhotonPt.end());
   filler_subsubleadingMuonVector_sumPFPhotonPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingMuonVector_sumPUPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingMuonVector_sumPUPt_valueMap(*subsubleadingMuonVector_sumPUPt_valueMap);
   filler_subsubleadingMuonVector_sumPUPt_valueMap.insert(TauHandle, subsubleadingMuonVector_sumPUPt.begin(), subsubleadingMuonVector_sumPUPt.end());
   filler_subsubleadingMuonVector_sumPUPt_valueMap.fill(); 

   std::unique_ptr< edm::ValueMap < float > > subsubleadingMuonVector_tausumPFChargedHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingMuonVector_tausumPFChargedHadronPt_valueMap(*subsubleadingMuonVector_tausumPFChargedHadronPt_valueMap);
   filler_subsubleadingMuonVector_tausumPFChargedHadronPt_valueMap.insert(TauHandle, subsubleadingMuonVector_tausumPFChargedHadronPt.begin(), subsubleadingMuonVector_tausumPFChargedHadronPt.end());
   filler_subsubleadingMuonVector_tausumPFChargedHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingMuonVector_tausumPFNeutralHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingMuonVector_tausumPFNeutralHadronPt_valueMap(*subsubleadingMuonVector_tausumPFNeutralHadronPt_valueMap);
   filler_subsubleadingMuonVector_tausumPFNeutralHadronPt_valueMap.insert(TauHandle, subsubleadingMuonVector_tausumPFNeutralHadronPt.begin(), subsubleadingMuonVector_tausumPFNeutralHadronPt.end());
   filler_subsubleadingMuonVector_tausumPFNeutralHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingMuonVector_tausumPFPhotonPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingMuonVector_tausumPFPhotonPt_valueMap(*subsubleadingMuonVector_tausumPFPhotonPt_valueMap);
   filler_subsubleadingMuonVector_tausumPFPhotonPt_valueMap.insert(TauHandle, subsubleadingMuonVector_tausumPFPhotonPt.begin(), subsubleadingMuonVector_tausumPFPhotonPt.end());
   filler_subsubleadingMuonVector_tausumPFPhotonPt_valueMap.fill();   



   //End of MUONS value Maps - Stuff that was not originally inlcuded by Andrew



  //Adding ELECTRONS value Maps - Stuff that was not originally inlcuded by Andrew

   std::unique_ptr< edm::ValueMap < float > > leadingElectronVector_sumPFChargedHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingElectronVector_sumPFChargedHadronPt_valueMap(*leadingElectronVector_sumPFChargedHadronPt_valueMap);
   filler_leadingElectronVector_sumPFChargedHadronPt_valueMap.insert(TauHandle, leadingElectronVector_sumPFChargedHadronPt.begin(), leadingElectronVector_sumPFChargedHadronPt.end());
   filler_leadingElectronVector_sumPFChargedHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingElectronVector_sumPFNeutralHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingElectronVector_sumPFNeutralHadronPt_valueMap(*leadingElectronVector_sumPFNeutralHadronPt_valueMap);
   filler_leadingElectronVector_sumPFNeutralHadronPt_valueMap.insert(TauHandle, leadingElectronVector_sumPFNeutralHadronPt.begin(), leadingElectronVector_sumPFNeutralHadronPt.end());
   filler_leadingElectronVector_sumPFNeutralHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingElectronVector_sumPFPhotonPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingElectronVector_sumPFPhotonPt_valueMap(*leadingElectronVector_sumPFPhotonPt_valueMap);
   filler_leadingElectronVector_sumPFPhotonPt_valueMap.insert(TauHandle, leadingElectronVector_sumPFPhotonPt.begin(), leadingElectronVector_sumPFPhotonPt.end());
   filler_leadingElectronVector_sumPFPhotonPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingElectronVector_rho_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingElectronVector_rho_valueMap(*leadingElectronVector_rho_valueMap);
   filler_leadingElectronVector_rho_valueMap.insert(TauHandle, leadingElectronVector_rho.begin(), leadingElectronVector_rho.end());
   filler_leadingElectronVector_rho_valueMap.fill(); 

   std::unique_ptr< edm::ValueMap < float > > leadingElectronVector_ea_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingElectronVector_ea_valueMap(*leadingElectronVector_ea_valueMap);
   filler_leadingElectronVector_ea_valueMap.insert(TauHandle, leadingElectronVector_ea.begin(), leadingElectronVector_ea.end());
   filler_leadingElectronVector_ea_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingElectronVector_tausumPFChargedHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingElectronVector_tausumPFChargedHadronPt_valueMap(*leadingElectronVector_tausumPFChargedHadronPt_valueMap);
   filler_leadingElectronVector_tausumPFChargedHadronPt_valueMap.insert(TauHandle, leadingElectronVector_tausumPFChargedHadronPt.begin(), leadingElectronVector_tausumPFChargedHadronPt.end());
   filler_leadingElectronVector_tausumPFChargedHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingElectronVector_tausumPFNeutralHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingElectronVector_tausumPFNeutralHadronPt_valueMap(*leadingElectronVector_tausumPFNeutralHadronPt_valueMap);
   filler_leadingElectronVector_tausumPFNeutralHadronPt_valueMap.insert(TauHandle, leadingElectronVector_tausumPFNeutralHadronPt.begin(), leadingElectronVector_tausumPFNeutralHadronPt.end());
   filler_leadingElectronVector_tausumPFNeutralHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingElectronVector_tausumPFPhotonPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingElectronVector_tausumPFPhotonPt_valueMap(*leadingElectronVector_tausumPFPhotonPt_valueMap);
   filler_leadingElectronVector_tausumPFPhotonPt_valueMap.insert(TauHandle, leadingElectronVector_tausumPFPhotonPt.begin(), leadingElectronVector_tausumPFPhotonPt.end());
   filler_leadingElectronVector_tausumPFPhotonPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingElectronVector_sumPFChargedHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingElectronVector_sumPFChargedHadronPt_valueMap(*subleadingElectronVector_sumPFChargedHadronPt_valueMap);
   filler_subleadingElectronVector_sumPFChargedHadronPt_valueMap.insert(TauHandle, subleadingElectronVector_sumPFChargedHadronPt.begin(), subleadingElectronVector_sumPFChargedHadronPt.end());
   filler_subleadingElectronVector_sumPFChargedHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingElectronVector_sumPFNeutralHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingElectronVector_sumPFNeutralHadronPt_valueMap(*subleadingElectronVector_sumPFNeutralHadronPt_valueMap);
   filler_subleadingElectronVector_sumPFNeutralHadronPt_valueMap.insert(TauHandle, subleadingElectronVector_sumPFNeutralHadronPt.begin(), subleadingElectronVector_sumPFNeutralHadronPt.end());
   filler_subleadingElectronVector_sumPFNeutralHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingElectronVector_sumPFPhotonPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingElectronVector_sumPFPhotonPt_valueMap(*subleadingElectronVector_sumPFPhotonPt_valueMap);
   filler_subleadingElectronVector_sumPFPhotonPt_valueMap.insert(TauHandle, subleadingElectronVector_sumPFPhotonPt.begin(), subleadingElectronVector_sumPFPhotonPt.end());
   filler_subleadingElectronVector_sumPFPhotonPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingElectronVector_rho_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingElectronVector_rho_valueMap(*subleadingElectronVector_rho_valueMap);
   filler_subleadingElectronVector_rho_valueMap.insert(TauHandle, subleadingElectronVector_rho.begin(), subleadingElectronVector_rho.end());
   filler_subleadingElectronVector_rho_valueMap.fill(); 

   std::unique_ptr< edm::ValueMap < float > > subleadingElectronVector_ea_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingElectronVector_ea_valueMap(*subleadingElectronVector_ea_valueMap);
   filler_subleadingElectronVector_ea_valueMap.insert(TauHandle, subleadingElectronVector_ea.begin(), subleadingElectronVector_ea.end());
   filler_subleadingElectronVector_ea_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingElectronVector_tausumPFChargedHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingElectronVector_tausumPFChargedHadronPt_valueMap(*subleadingElectronVector_tausumPFChargedHadronPt_valueMap);
   filler_subleadingElectronVector_tausumPFChargedHadronPt_valueMap.insert(TauHandle, subleadingElectronVector_tausumPFChargedHadronPt.begin(), subleadingElectronVector_tausumPFChargedHadronPt.end());
   filler_subleadingElectronVector_tausumPFChargedHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingElectronVector_tausumPFNeutralHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingElectronVector_tausumPFNeutralHadronPt_valueMap(*subleadingElectronVector_tausumPFNeutralHadronPt_valueMap);
   filler_subleadingElectronVector_tausumPFNeutralHadronPt_valueMap.insert(TauHandle, subleadingElectronVector_tausumPFNeutralHadronPt.begin(), subleadingElectronVector_tausumPFNeutralHadronPt.end());
   filler_subleadingElectronVector_tausumPFNeutralHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingElectronVector_tausumPFPhotonPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingElectronVector_tausumPFPhotonPt_valueMap(*subleadingElectronVector_tausumPFPhotonPt_valueMap);
   filler_subleadingElectronVector_tausumPFPhotonPt_valueMap.insert(TauHandle, subleadingElectronVector_tausumPFPhotonPt.begin(), subleadingElectronVector_tausumPFPhotonPt.end());
   filler_subleadingElectronVector_tausumPFPhotonPt_valueMap.fill();     


   std::unique_ptr< edm::ValueMap < float > > subsubleadingElectronVector_sumPFChargedHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingElectronVector_sumPFChargedHadronPt_valueMap(*subsubleadingElectronVector_sumPFChargedHadronPt_valueMap);
   filler_subsubleadingElectronVector_sumPFChargedHadronPt_valueMap.insert(TauHandle, subsubleadingElectronVector_sumPFChargedHadronPt.begin(), subsubleadingElectronVector_sumPFChargedHadronPt.end());
   filler_subsubleadingElectronVector_sumPFChargedHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingElectronVector_sumPFNeutralHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingElectronVector_sumPFNeutralHadronPt_valueMap(*subsubleadingElectronVector_sumPFNeutralHadronPt_valueMap);
   filler_subsubleadingElectronVector_sumPFNeutralHadronPt_valueMap.insert(TauHandle, subsubleadingElectronVector_sumPFNeutralHadronPt.begin(), subsubleadingElectronVector_sumPFNeutralHadronPt.end());
   filler_subsubleadingElectronVector_sumPFNeutralHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingElectronVector_sumPFPhotonPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingElectronVector_sumPFPhotonPt_valueMap(*subsubleadingElectronVector_sumPFPhotonPt_valueMap);
   filler_subsubleadingElectronVector_sumPFPhotonPt_valueMap.insert(TauHandle, subsubleadingElectronVector_sumPFPhotonPt.begin(), subsubleadingElectronVector_sumPFPhotonPt.end());
   filler_subsubleadingElectronVector_sumPFPhotonPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingElectronVector_rho_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingElectronVector_rho_valueMap(*subsubleadingElectronVector_rho_valueMap);
   filler_subsubleadingElectronVector_rho_valueMap.insert(TauHandle, subsubleadingElectronVector_rho.begin(), subsubleadingElectronVector_rho.end());
   filler_subsubleadingElectronVector_rho_valueMap.fill(); 

   std::unique_ptr< edm::ValueMap < float > > subsubleadingElectronVector_ea_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingElectronVector_ea_valueMap(*subsubleadingElectronVector_ea_valueMap);
   filler_subsubleadingElectronVector_ea_valueMap.insert(TauHandle, subsubleadingElectronVector_ea.begin(), subsubleadingElectronVector_ea.end());
   filler_subsubleadingElectronVector_ea_valueMap.fill(); 

   std::unique_ptr< edm::ValueMap < float > > subsubleadingElectronVector_tausumPFChargedHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingElectronVector_tausumPFChargedHadronPt_valueMap(*subsubleadingElectronVector_tausumPFChargedHadronPt_valueMap);
   filler_subsubleadingElectronVector_tausumPFChargedHadronPt_valueMap.insert(TauHandle, subsubleadingElectronVector_tausumPFChargedHadronPt.begin(), subsubleadingElectronVector_tausumPFChargedHadronPt.end());
   filler_subsubleadingElectronVector_tausumPFChargedHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingElectronVector_tausumPFNeutralHadronPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingElectronVector_tausumPFNeutralHadronPt_valueMap(*subsubleadingElectronVector_tausumPFNeutralHadronPt_valueMap);
   filler_subsubleadingElectronVector_tausumPFNeutralHadronPt_valueMap.insert(TauHandle, subsubleadingElectronVector_tausumPFNeutralHadronPt.begin(), subsubleadingElectronVector_tausumPFNeutralHadronPt.end());
   filler_subsubleadingElectronVector_tausumPFNeutralHadronPt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingElectronVector_tausumPFPhotonPt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingElectronVector_tausumPFPhotonPt_valueMap(*subsubleadingElectronVector_tausumPFPhotonPt_valueMap);
   filler_subsubleadingElectronVector_tausumPFPhotonPt_valueMap.insert(TauHandle, subsubleadingElectronVector_tausumPFPhotonPt.begin(), subsubleadingElectronVector_tausumPFPhotonPt.end());
   filler_subsubleadingElectronVector_tausumPFPhotonPt_valueMap.fill(); 

     //End of ELECTRONS value Maps - Stuff that was not originally inlcuded by Andrew

   ///////////

   //MUON value maps - Originally included by Andrew- total correction, kinematic variables and new delRs

     //delR for muon valuemaps (for all three cases leading, subleading, subsubleading)


   std::unique_ptr< edm::ValueMap < float > > leadingMuonVector_delR_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingMuonVector_delR_valueMap(*leadingMuonVector_delR_valueMap);
   filler_leadingMuonVector_delR_valueMap.insert(TauHandle, leadingMuonVector_delR.begin(), leadingMuonVector_delR.end());
   filler_leadingMuonVector_delR_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingMuonVector_delR_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingMuonVector_delR_valueMap(*subleadingMuonVector_delR_valueMap);
   filler_subleadingMuonVector_delR_valueMap.insert(TauHandle, subleadingMuonVector_delR.begin(), subleadingMuonVector_delR.end());
   filler_subleadingMuonVector_delR_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingMuonVector_delR_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingMuonVector_delR_valueMap(*subsubleadingMuonVector_delR_valueMap);
   filler_subsubleadingMuonVector_delR_valueMap.insert(TauHandle, subsubleadingMuonVector_delR.begin(), subsubleadingMuonVector_delR.end());
   filler_subsubleadingMuonVector_delR_valueMap.fill();
     //end

   std::unique_ptr< edm::ValueMap < int > > Mcounter_valueMap(new edm::ValueMap < int >());
   edm::ValueMap< int >::Filler filler_Mcounter_valueMap(*Mcounter_valueMap);
   filler_Mcounter_valueMap.insert(TauHandle, Mcounter.begin(), Mcounter.end());
   filler_Mcounter_valueMap.fill();   

   std::unique_ptr< edm::ValueMap < int > > Ecounter_valueMap(new edm::ValueMap < int >());
   edm::ValueMap< int >::Filler filler_Ecounter_valueMap(*Ecounter_valueMap);
   filler_Ecounter_valueMap.insert(TauHandle, Ecounter.begin(), Ecounter.end());
   filler_Ecounter_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingMuonVector_pt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingMuonVector_pt_valueMap(*leadingMuonVector_pt_valueMap);
   filler_leadingMuonVector_pt_valueMap.insert(TauHandle, leadingMuonVector_pt.begin(), leadingMuonVector_pt.end());
   filler_leadingMuonVector_pt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingMuonVector_eta_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingMuonVector_eta_valueMap(*leadingMuonVector_eta_valueMap);
   filler_leadingMuonVector_eta_valueMap.insert(TauHandle, leadingMuonVector_eta.begin(), leadingMuonVector_eta.end());
   filler_leadingMuonVector_eta_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingMuonVector_phi_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingMuonVector_phi_valueMap(*leadingMuonVector_phi_valueMap);
   filler_leadingMuonVector_phi_valueMap.insert(TauHandle, leadingMuonVector_phi.begin(), leadingMuonVector_phi.end());
   filler_leadingMuonVector_phi_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingMuonVector_m_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingMuonVector_m_valueMap(*leadingMuonVector_m_valueMap);
   filler_leadingMuonVector_m_valueMap.insert(TauHandle, leadingMuonVector_m.begin(), leadingMuonVector_m.end());
   filler_leadingMuonVector_m_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingMuonVector_corrIso_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingMuonVector_corrIso_valueMap(*leadingMuonVector_corrIso_valueMap);
   filler_leadingMuonVector_corrIso_valueMap.insert(TauHandle, leadingMuonVector_corrIso.begin(), leadingMuonVector_corrIso.end());
   filler_leadingMuonVector_corrIso_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingMuonVector_pt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingMuonVector_pt_valueMap(*subleadingMuonVector_pt_valueMap);
   filler_subleadingMuonVector_pt_valueMap.insert(TauHandle, subleadingMuonVector_pt.begin(), subleadingMuonVector_pt.end());
   filler_subleadingMuonVector_pt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingMuonVector_eta_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingMuonVector_eta_valueMap(*subleadingMuonVector_eta_valueMap);
   filler_subleadingMuonVector_eta_valueMap.insert(TauHandle, subleadingMuonVector_eta.begin(), subleadingMuonVector_eta.end());
   filler_subleadingMuonVector_eta_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingMuonVector_phi_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingMuonVector_phi_valueMap(*subleadingMuonVector_phi_valueMap);
   filler_subleadingMuonVector_phi_valueMap.insert(TauHandle, subleadingMuonVector_phi.begin(), subleadingMuonVector_phi.end());
   filler_subleadingMuonVector_phi_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingMuonVector_m_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingMuonVector_m_valueMap(*subleadingMuonVector_m_valueMap);
   filler_subleadingMuonVector_m_valueMap.insert(TauHandle, subleadingMuonVector_m.begin(), subleadingMuonVector_m.end());
   filler_subleadingMuonVector_m_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingMuonVector_corrIso_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingMuonVector_corrIso_valueMap(*subleadingMuonVector_corrIso_valueMap);
   filler_subleadingMuonVector_corrIso_valueMap.insert(TauHandle, subleadingMuonVector_corrIso.begin(), subleadingMuonVector_corrIso.end());
   filler_subleadingMuonVector_corrIso_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingMuonVector_pt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingMuonVector_pt_valueMap(*subsubleadingMuonVector_pt_valueMap);
   filler_subsubleadingMuonVector_pt_valueMap.insert(TauHandle, subsubleadingMuonVector_pt.begin(), subsubleadingMuonVector_pt.end());
   filler_subsubleadingMuonVector_pt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingMuonVector_eta_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingMuonVector_eta_valueMap(*subsubleadingMuonVector_eta_valueMap);
   filler_subsubleadingMuonVector_eta_valueMap.insert(TauHandle, subsubleadingMuonVector_eta.begin(), subsubleadingMuonVector_eta.end());
   filler_subsubleadingMuonVector_eta_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingMuonVector_phi_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingMuonVector_phi_valueMap(*subsubleadingMuonVector_phi_valueMap);
   filler_subsubleadingMuonVector_phi_valueMap.insert(TauHandle, subsubleadingMuonVector_phi.begin(), subsubleadingMuonVector_phi.end());
   filler_subsubleadingMuonVector_phi_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingMuonVector_m_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingMuonVector_m_valueMap(*subsubleadingMuonVector_m_valueMap);
   filler_subsubleadingMuonVector_m_valueMap.insert(TauHandle, subsubleadingMuonVector_m.begin(), subsubleadingMuonVector_m.end());
   filler_subsubleadingMuonVector_m_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingMuonVector_corrIso_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingMuonVector_corrIso_valueMap(*subsubleadingMuonVector_corrIso_valueMap);
   filler_subsubleadingMuonVector_corrIso_valueMap.insert(TauHandle, subsubleadingMuonVector_corrIso.begin(), subsubleadingMuonVector_corrIso.end());
   filler_subsubleadingMuonVector_corrIso_valueMap.fill();


     //End of MUON value maps - Originally included by Andrew- total correction, kinematic variables and new delRs

   ////////////////

   //ELECTRONS value maps - Originally included by Andrew- total correction, kinematic variables and new delRs
  
  //delR for electrons valuemaps (for all three cases leading, subleading and subsubleading)
   std::unique_ptr< edm::ValueMap < float > > leadingElectronVector_delR_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingElectronVector_delR_valueMap(*leadingElectronVector_delR_valueMap);
   filler_leadingElectronVector_delR_valueMap.insert(TauHandle, leadingElectronVector_delR.begin(), leadingElectronVector_delR.end());
   filler_leadingElectronVector_delR_valueMap.fill();  

   std::unique_ptr< edm::ValueMap < float > > subleadingElectronVector_delR_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingElectronVector_delR_valueMap(*subleadingElectronVector_delR_valueMap);
   filler_subleadingElectronVector_delR_valueMap.insert(TauHandle, subleadingElectronVector_delR.begin(), subleadingElectronVector_delR.end());
   filler_subleadingElectronVector_delR_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingElectronVector_delR_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingElectronVector_delR_valueMap(*subsubleadingElectronVector_delR_valueMap);
   filler_subsubleadingElectronVector_delR_valueMap.insert(TauHandle, subsubleadingElectronVector_delR.begin(), subsubleadingElectronVector_delR.end());
   filler_subsubleadingElectronVector_delR_valueMap.fill();
   //end

   std::unique_ptr< edm::ValueMap < float > > leadingElectronVector_pt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingElectronVector_pt_valueMap(*leadingElectronVector_pt_valueMap);
   filler_leadingElectronVector_pt_valueMap.insert(TauHandle, leadingElectronVector_pt.begin(), leadingElectronVector_pt.end());
   filler_leadingElectronVector_pt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingElectronVector_eta_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingElectronVector_eta_valueMap(*leadingElectronVector_eta_valueMap);
   filler_leadingElectronVector_eta_valueMap.insert(TauHandle, leadingElectronVector_eta.begin(), leadingElectronVector_eta.end());
   filler_leadingElectronVector_eta_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingElectronVector_phi_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingElectronVector_phi_valueMap(*leadingElectronVector_phi_valueMap);
   filler_leadingElectronVector_phi_valueMap.insert(TauHandle, leadingElectronVector_phi.begin(), leadingElectronVector_phi.end());
   filler_leadingElectronVector_phi_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingElectronVector_m_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingElectronVector_m_valueMap(*leadingElectronVector_m_valueMap);
   filler_leadingElectronVector_m_valueMap.insert(TauHandle, leadingElectronVector_m.begin(), leadingElectronVector_m.end());
   filler_leadingElectronVector_m_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > leadingElectronVector_corrIso_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_leadingElectronVector_corrIso_valueMap(*leadingElectronVector_corrIso_valueMap);
   filler_leadingElectronVector_corrIso_valueMap.insert(TauHandle, leadingElectronVector_corrIso.begin(), leadingElectronVector_corrIso.end());
   filler_leadingElectronVector_corrIso_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingElectronVector_pt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingElectronVector_pt_valueMap(*subleadingElectronVector_pt_valueMap);
   filler_subleadingElectronVector_pt_valueMap.insert(TauHandle, subleadingElectronVector_pt.begin(), subleadingElectronVector_pt.end());
   filler_subleadingElectronVector_pt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingElectronVector_eta_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingElectronVector_eta_valueMap(*subleadingElectronVector_eta_valueMap);
   filler_subleadingElectronVector_eta_valueMap.insert(TauHandle, subleadingElectronVector_eta.begin(), subleadingElectronVector_eta.end());
   filler_subleadingElectronVector_eta_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingElectronVector_phi_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingElectronVector_phi_valueMap(*subleadingElectronVector_phi_valueMap);
   filler_subleadingElectronVector_phi_valueMap.insert(TauHandle, subleadingElectronVector_phi.begin(), subleadingElectronVector_phi.end());
   filler_subleadingElectronVector_phi_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingElectronVector_m_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingElectronVector_m_valueMap(*subleadingElectronVector_m_valueMap);
   filler_subleadingElectronVector_m_valueMap.insert(TauHandle, subleadingElectronVector_m.begin(), subleadingElectronVector_m.end());
   filler_subleadingElectronVector_m_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subleadingElectronVector_corrIso_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subleadingElectronVector_corrIso_valueMap(*subleadingElectronVector_corrIso_valueMap);
   filler_subleadingElectronVector_corrIso_valueMap.insert(TauHandle, subleadingElectronVector_corrIso.begin(), subleadingElectronVector_corrIso.end());
   filler_subleadingElectronVector_corrIso_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingElectronVector_pt_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingElectronVector_pt_valueMap(*subsubleadingElectronVector_pt_valueMap);
   filler_subsubleadingElectronVector_pt_valueMap.insert(TauHandle, subsubleadingElectronVector_pt.begin(), subsubleadingElectronVector_pt.end());
   filler_subsubleadingElectronVector_pt_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingElectronVector_eta_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingElectronVector_eta_valueMap(*subsubleadingElectronVector_eta_valueMap);
   filler_subsubleadingElectronVector_eta_valueMap.insert(TauHandle, subsubleadingElectronVector_eta.begin(), subsubleadingElectronVector_eta.end());
   filler_subsubleadingElectronVector_eta_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingElectronVector_phi_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingElectronVector_phi_valueMap(*subsubleadingElectronVector_phi_valueMap);
   filler_subsubleadingElectronVector_phi_valueMap.insert(TauHandle, subsubleadingElectronVector_phi.begin(), subsubleadingElectronVector_phi.end());
   filler_subsubleadingElectronVector_phi_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingElectronVector_m_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingElectronVector_m_valueMap(*subsubleadingElectronVector_m_valueMap);
   filler_subsubleadingElectronVector_m_valueMap.insert(TauHandle, subsubleadingElectronVector_m.begin(), subsubleadingElectronVector_m.end());
   filler_subsubleadingElectronVector_m_valueMap.fill();

   std::unique_ptr< edm::ValueMap < float > > subsubleadingElectronVector_corrIso_valueMap(new edm::ValueMap < float >());
   edm::ValueMap< float >::Filler filler_subsubleadingElectronVector_corrIso_valueMap(*subsubleadingElectronVector_corrIso_valueMap);
   filler_subsubleadingElectronVector_corrIso_valueMap.insert(TauHandle, subsubleadingElectronVector_corrIso.begin(), subsubleadingElectronVector_corrIso.end());
   filler_subsubleadingElectronVector_corrIso_valueMap.fill();

   //End of ELECTRONS value maps - Originally included by Andrew- total correction, kinematic variables and new delRs

   iEvent.put(std::move(Mcounter_valueMap), "Mcounter");
   iEvent.put(std::move(Ecounter_valueMap), "Ecounter");

   //put Muon delR in the event for all three cases
   iEvent.put(std::move(leadingMuonVector_delR_valueMap), "LeadingMuondelR");
   iEvent.put(std::move(subleadingMuonVector_delR_valueMap), "SubLeadingMuondelR");
   iEvent.put(std::move(subsubleadingMuonVector_delR_valueMap), "SubSubLeadingMuondelR");
   //end

   iEvent.put(std::move(leadingMuonVector_pt_valueMap), "LeadingMuonPt");
   iEvent.put(std::move(leadingMuonVector_eta_valueMap), "LeadingMuonEta");
   iEvent.put(std::move(leadingMuonVector_phi_valueMap), "LeadingMuonPhi");
   iEvent.put(std::move(leadingMuonVector_m_valueMap), "LeadingMuonM");
   iEvent.put(std::move(leadingMuonVector_corrIso_valueMap), "LeadingMuonCorrIso");
   iEvent.put(std::move(leadingMuonVector_sumPFChargedHadronPt_valueMap), "LeadingMuonsumPFChargedHadronPt");
   iEvent.put(std::move(leadingMuonVector_sumPFNeutralHadronPt_valueMap), "LeadingMuonsumPFNeutralHadronPt");
   iEvent.put(std::move(leadingMuonVector_sumPFPhotonPt_valueMap), "LeadingMuonsumPFPhotonPt");
   iEvent.put(std::move(leadingMuonVector_sumPUPt_valueMap), "LeadingMuonsumPUPt");
   iEvent.put(std::move(leadingMuonVector_tausumPFChargedHadronPt_valueMap), "LeadingMuontausumPFChargedHadronPt");
   iEvent.put(std::move(leadingMuonVector_tausumPFNeutralHadronPt_valueMap), "LeadingMuontausumPFNeutralHadronPt");
   iEvent.put(std::move(leadingMuonVector_tausumPFPhotonPt_valueMap), "LeadingMuontausumPFPhotonPt");



   iEvent.put(std::move(subleadingMuonVector_pt_valueMap), "SubLeadingMuonPt");
   iEvent.put(std::move(subleadingMuonVector_eta_valueMap), "SubLeadingMuonEta");
   iEvent.put(std::move(subleadingMuonVector_phi_valueMap), "SubLeadingMuonPhi");
   iEvent.put(std::move(subleadingMuonVector_m_valueMap), "SubLeadingMuonM");
   iEvent.put(std::move(subleadingMuonVector_corrIso_valueMap), "SubLeadingMuonCorrIso");
   iEvent.put(std::move(subleadingMuonVector_sumPFChargedHadronPt_valueMap), "SubLeadingMuonsumPFChargedHadronPt");
   iEvent.put(std::move(subleadingMuonVector_sumPFNeutralHadronPt_valueMap), "SubLeadingMuonsumPFNeutralHadronPt");
   iEvent.put(std::move(subleadingMuonVector_sumPFPhotonPt_valueMap), "SubLeadingMuonsumPFPhotonPt");
   iEvent.put(std::move(subleadingMuonVector_sumPUPt_valueMap), "SubLeadingMuonsumPUPt");
   iEvent.put(std::move(subleadingMuonVector_tausumPFChargedHadronPt_valueMap), "SubLeadingMuontausumPFChargedHadronPt");
   iEvent.put(std::move(subleadingMuonVector_tausumPFNeutralHadronPt_valueMap), "SubLeadingMuontausumPFNeutralHadronPt");
   iEvent.put(std::move(subleadingMuonVector_tausumPFPhotonPt_valueMap), "SubLeadingMuontausumPFPhotonPt");   


   iEvent.put(std::move(subsubleadingMuonVector_pt_valueMap), "SubSubLeadingMuonPt");
   iEvent.put(std::move(subsubleadingMuonVector_eta_valueMap), "SubSubLeadingMuonEta");
   iEvent.put(std::move(subsubleadingMuonVector_phi_valueMap), "SubSubLeadingMuonPhi");
   iEvent.put(std::move(subsubleadingMuonVector_m_valueMap), "SubSubLeadingMuonM");
   iEvent.put(std::move(subsubleadingMuonVector_corrIso_valueMap), "SubSubLeadingMuonCorrIso");
   iEvent.put(std::move(subsubleadingMuonVector_sumPFChargedHadronPt_valueMap), "SubSubLeadingMuonsumPFChargedHadronPt");
   iEvent.put(std::move(subsubleadingMuonVector_sumPFNeutralHadronPt_valueMap), "SubSubLeadingMuonsumPFNeutralHadronPt");
   iEvent.put(std::move(subsubleadingMuonVector_sumPFPhotonPt_valueMap), "SubSubLeadingMuonsumPFPhotonPt");
   iEvent.put(std::move(subsubleadingMuonVector_sumPUPt_valueMap), "SubSubLeadingMuonsumPUPt");
   iEvent.put(std::move(subsubleadingMuonVector_tausumPFChargedHadronPt_valueMap), "SubSubLeadingMuontausumPFChargedHadronPt");
   iEvent.put(std::move(subsubleadingMuonVector_tausumPFNeutralHadronPt_valueMap), "SubSubLeadingMuontausumPFNeutralHadronPt");
   iEvent.put(std::move(subsubleadingMuonVector_tausumPFPhotonPt_valueMap), "SubSubLeadingMuontausumPFPhotonPt");


     //put Electron delR in the event for all three cases
   
   iEvent.put(std::move(leadingElectronVector_delR_valueMap), "LeadingElectrondelR");
   iEvent.put(std::move(subleadingElectronVector_delR_valueMap), "SubLeadingElectrondelR");
   iEvent.put(std::move(subsubleadingElectronVector_delR_valueMap), "SubSubLeadingElectrondelR");
    //END
    

   iEvent.put(std::move(leadingElectronVector_pt_valueMap), "LeadingElectronPt");
   iEvent.put(std::move(leadingElectronVector_eta_valueMap), "LeadingElectronEta");
   iEvent.put(std::move(leadingElectronVector_phi_valueMap), "LeadingElectronPhi");
   iEvent.put(std::move(leadingElectronVector_m_valueMap), "LeadingElectronM");
   iEvent.put(std::move(leadingElectronVector_corrIso_valueMap), "LeadingElectronCorrIso");
   iEvent.put(std::move(leadingElectronVector_sumPFChargedHadronPt_valueMap), "LeadingElectronsumPFChargedHadronPt");
   iEvent.put(std::move(leadingElectronVector_sumPFNeutralHadronPt_valueMap), "LeadingElectronsumPFNeutralHadronPt");
   iEvent.put(std::move(leadingElectronVector_sumPFPhotonPt_valueMap), "LeadingElectronsumPFPhotonPt");
   iEvent.put(std::move(leadingElectronVector_tausumPFChargedHadronPt_valueMap), "LeadingElectrontausumPFChargedHadronPt");
   iEvent.put(std::move(leadingElectronVector_tausumPFNeutralHadronPt_valueMap), "LeadingElectrontausumPFNeutralHadronPt");
   iEvent.put(std::move(leadingElectronVector_tausumPFPhotonPt_valueMap), "LeadingElectrontausumPFPhotonPt");
   iEvent.put(std::move(leadingElectronVector_rho_valueMap), "LeadingElectronrho");
   iEvent.put(std::move(leadingElectronVector_ea_valueMap), "LeadingElectronea");

   iEvent.put(std::move(subleadingElectronVector_pt_valueMap), "SubLeadingElectronPt");
   iEvent.put(std::move(subleadingElectronVector_eta_valueMap), "SubLeadingElectronEta");
   iEvent.put(std::move(subleadingElectronVector_phi_valueMap), "SubLeadingElectronPhi");
   iEvent.put(std::move(subleadingElectronVector_m_valueMap), "SubLeadingElectronM");
   iEvent.put(std::move(subleadingElectronVector_corrIso_valueMap), "SubLeadingElectronCorrIso");
   iEvent.put(std::move(subleadingElectronVector_sumPFChargedHadronPt_valueMap), "SubLeadingElectronsumPFChargedHadronPt");
   iEvent.put(std::move(subleadingElectronVector_sumPFNeutralHadronPt_valueMap), "SubLeadingElectronsumPFNeutralHadronPt");
   iEvent.put(std::move(subleadingElectronVector_sumPFPhotonPt_valueMap), "SubLeadingElectronsumPFPhotonPt");
   iEvent.put(std::move(subleadingElectronVector_tausumPFChargedHadronPt_valueMap), "SubLeadingElectrontausumPFChargedHadronPt");
   iEvent.put(std::move(subleadingElectronVector_tausumPFNeutralHadronPt_valueMap), "SubLeadingElectrontausumPFNeutralHadronPt");
   iEvent.put(std::move(subleadingElectronVector_tausumPFPhotonPt_valueMap), "SubLeadingElectrontausumPFPhotonPt");
   iEvent.put(std::move(subleadingElectronVector_rho_valueMap), "SubLeadingElectronrho");
   iEvent.put(std::move(subleadingElectronVector_ea_valueMap), "SubLeadingElectronea");   


   iEvent.put(std::move(subsubleadingElectronVector_pt_valueMap), "SubSubLeadingElectronPt");
   iEvent.put(std::move(subsubleadingElectronVector_eta_valueMap), "SubSubLeadingElectronEta");
   iEvent.put(std::move(subsubleadingElectronVector_phi_valueMap), "SubSubLeadingElectronPhi");
   iEvent.put(std::move(subsubleadingElectronVector_m_valueMap), "SubSubLeadingElectronM");
   iEvent.put(std::move(subsubleadingElectronVector_corrIso_valueMap), "SubSubLeadingElectronCorrIso");
   iEvent.put(std::move(subsubleadingElectronVector_sumPFChargedHadronPt_valueMap), "SubSubLeadingElectronsumPFChargedHadronPt");
   iEvent.put(std::move(subsubleadingElectronVector_sumPFNeutralHadronPt_valueMap), "SubSubLeadingElectronsumPFNeutralHadronPt");
   iEvent.put(std::move(subsubleadingElectronVector_sumPFPhotonPt_valueMap), "SubSubLeadingElectronsumPFPhotonPt");
   iEvent.put(std::move(subsubleadingElectronVector_tausumPFChargedHadronPt_valueMap), "SubSubLeadingElectrontausumPFChargedHadronPt");
   iEvent.put(std::move(subsubleadingElectronVector_tausumPFNeutralHadronPt_valueMap), "SubSubLeadingElectrontausumPFNeutralHadronPt");
   iEvent.put(std::move(subsubleadingElectronVector_tausumPFPhotonPt_valueMap), "SubSubLeadingElectrontausumPFPhotonPt");
   iEvent.put(std::move(subsubleadingElectronVector_rho_valueMap), "SubSubLeadingElectronrho");
   iEvent.put(std::move(subsubleadingElectronVector_ea_valueMap), "SubSubLeadingElectronea");   
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
TauLeadingLeptonIso::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
TauLeadingLeptonIso::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
TauLeadingLeptonIso::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
TauLeadingLeptonIso::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TauLeadingLeptonIso::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TauLeadingLeptonIso::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauLeadingLeptonIso::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void TauLeadingLeptonIso::calculateCorrectedMuonIsoInformation(pat::Tau theTau, leptonInfo &theMuonInfo)
{
  double tauSumChargedHadronPt = 0.0;
  double tauSumNeutralHadronEt = 0.0;
  double tauSumPhotonEt        = 0.0;
  //Okay, let's go through and check for any charged hadrons of the tau in the cone of the muon
  //If we find any, we add it to the sum and correct the isolation for it
  
  for(size_t hadrCandInd = 0;
      hadrCandInd < theTau.signalChargedHadrCands().size();
      ++hadrCandInd)
    {
      double dRConst = reco::deltaR(theMuonInfo.eta, theMuonInfo.phi, theTau.signalChargedHadrCands()[hadrCandInd]->eta(), theTau.signalChargedHadrCands()[hadrCandInd]->phi());
      if (dRConst < 0.4) tauSumChargedHadronPt += theTau.signalChargedHadrCands()[hadrCandInd]->pt();
    }

  for(size_t neutCandInd = 0;
      neutCandInd < theTau.signalNeutrHadrCands().size();
      ++neutCandInd)
    {
      double dRConst = reco::deltaR(theMuonInfo.eta, theMuonInfo.phi, theTau.signalNeutrHadrCands()[neutCandInd]->eta(), theTau.signalNeutrHadrCands()[neutCandInd]->phi());
      if (dRConst < 0.4) tauSumNeutralHadronEt += theTau.signalNeutrHadrCands()[neutCandInd]->pt();
    }

  for(size_t photonCandInd = 0;
      photonCandInd < theTau.signalGammaCands().size();
      ++photonCandInd)
    {
      double dRConst = reco::deltaR(theMuonInfo.eta, theMuonInfo.phi, theTau.signalGammaCands()[photonCandInd]->eta(), theTau.signalGammaCands()[photonCandInd]->phi());
      if (dRConst < 0.4) tauSumPhotonEt += theTau.signalGammaCands()[photonCandInd]->pt();
    }
  
  theMuonInfo.tausumPFChargedHadronPt = tauSumChargedHadronPt;
  theMuonInfo.tausumPFNeutralHadronPt = tauSumNeutralHadronEt;
  theMuonInfo.tausumPFPhotonPt = tauSumPhotonEt;


  theMuonInfo.correctedSumPFChargedHadronPt = std::max(0.0, theMuonInfo.sumPFChargedHadronPt - tauSumChargedHadronPt);
  theMuonInfo.correctedSumPFNeutralHadronPt = std::max(0.0, theMuonInfo.sumPFNeutralHadronPt - tauSumNeutralHadronEt + theMuonInfo.sumPFPhotonPt - tauSumPhotonEt);

  theMuonInfo.correctedIso = theMuonInfo.correctedSumPFChargedHadronPt + std::max(theMuonInfo.correctedSumPFNeutralHadronPt - theMuonInfo.sumPUPt/2.0, 0.0);
}

void TauLeadingLeptonIso::calculateCorrectedElectronIsoInformation(pat::Tau theTau, leptonInfo &theElectronInfo)
{
  double tauSumChargedHadronPt = 0.0;
  double tauSumNeutralHadronEt = 0.0;
  double tauSumPhotonEt        = 0.0;
  //Okay, let's go through and check for any charged hadrons of the tau in the cone of the muon
  //If we find any, we add it to the sum and correct the isolation for it
  
  for(size_t hadrCandInd = 0;
      hadrCandInd < theTau.signalChargedHadrCands().size();
      ++hadrCandInd)
    {
      double dRConst = reco::deltaR(theElectronInfo.eta, theElectronInfo.phi, theTau.signalChargedHadrCands()[hadrCandInd]->eta(), theTau.signalChargedHadrCands()[hadrCandInd]->phi());
      if (dRConst < 0.3) tauSumChargedHadronPt += theTau.signalChargedHadrCands()[hadrCandInd]->pt();
    }

  for(size_t neutCandInd = 0;
      neutCandInd < theTau.signalNeutrHadrCands().size();
      ++neutCandInd)
    {
      double dRConst = reco::deltaR(theElectronInfo.eta, theElectronInfo.phi, theTau.signalNeutrHadrCands()[neutCandInd]->eta(), theTau.signalNeutrHadrCands()[neutCandInd]->phi());
      if (dRConst < 0.3) tauSumNeutralHadronEt += theTau.signalNeutrHadrCands()[neutCandInd]->pt();
    }

  for(size_t photonCandInd = 0;
      photonCandInd < theTau.signalGammaCands().size();
      ++photonCandInd)
    {
      double dRConst = reco::deltaR(theElectronInfo.eta, theElectronInfo.phi, theTau.signalGammaCands()[photonCandInd]->eta(), theTau.signalGammaCands()[photonCandInd]->phi());
      if (dRConst < 0.3) tauSumPhotonEt += theTau.signalGammaCands()[photonCandInd]->pt();
    }
  
  theElectronInfo.tausumPFChargedHadronPt = tauSumChargedHadronPt;
  theElectronInfo.tausumPFNeutralHadronPt = tauSumNeutralHadronEt;
  theElectronInfo.tausumPFPhotonPt = tauSumPhotonEt;
  
  
  theElectronInfo.correctedSumPFChargedHadronPt = std::max(0.0, theElectronInfo.sumPFChargedHadronPt - tauSumChargedHadronPt);
  theElectronInfo.correctedSumPFNeutralHadronPt = std::max(0.0, theElectronInfo.sumPFNeutralHadronPt - tauSumNeutralHadronEt + theElectronInfo.sumPFPhotonPt - tauSumPhotonEt);
  
  theElectronInfo.correctedIso = theElectronInfo.correctedSumPFChargedHadronPt +std::max((double)(theElectronInfo.correctedSumPFNeutralHadronPt - theElectronInfo.rho * theElectronInfo.ea), 0.0);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauLeadingLeptonIso);
