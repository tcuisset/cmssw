#ifndef SimGeneral_CaloAnalysis_SimClusterMergerByFastJet_h
#define SimGeneral_CaloAnalysis_SimClusterMergerByFastJet_h

#include <vector>
#include <ranges>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "SimGeneral/MixingModule/interface/DecayGraph.h"

#include <fastjet/ClusterSequence.hh>

#include "SimGeneral/CaloAnalysis/interface/SimClusterBuilder.h"

namespace {
template <typename ClusterParentIndexRecorderT>
class SimClusterMergerByFastJet : public SubClusterMergerBase {
public:
SimClusterMergerByFastJet(SimClusterCollection &clusters,
                    ClusterParentIndexRecorderT parentIndexRecorder,
                    fastjet::JetDefinition const& jetDefinition)
    : clusters_(clusters), indexRecorder(parentIndexRecorder),jetDefinition_(jetDefinition) {}

void start_subcluster(DecayChain::edge_descriptor e, const DecayChain &g, std::size_t currentSubClusterIndex) {
    auto edge_property = get(edge_weight, g, e);
    SimTrack const& simtrack = *edge_property.simTrack;

    /* Build the particle 4-vector such that the energy is the energy of SimTrack at boundary,
    the momentum 3-vector points to the boundary position, and the mass is zero */
    auto energyAtBoundary = simtrack.getMomentumAtBoundary().E();
    auto momentum3D = energyAtBoundary * simtrack.getPositionAtBoundary().Vect().Unit();
    auto& jet = fjInputs.emplace_back(momentum3D.X(), momentum3D.Y(), momentum3D.Z(), energyAtBoundary);
    jet.set_user_index(currentSubClusterIndex); // Store the index in SimCLuster collection for later
}

void end_parentCluster(std::span<SimCluster> subClustersBuilt, std::size_t currentCaloParticleIndex) {
    // Clustering
    fastjet::ClusterSequence sequence(fjInputs, jetDefinition_);
    auto jets = sequence.inclusive_jets();

    // Merging
    for (fastjet::PseudoJet const& jet : jets) {
        auto constituents = fastjet::sorted_by_E(jet.constituents());
        assert(constituents.size() >= 1);
        assert(constituents[0].user_index()<static_cast<int>(subClustersBuilt.size()));
        SimCluster& mergedSimCluster = clusters_.emplace_back(SimCluster::mergeHitsFromCollection(constituents | std::views::transform([&](fastjet::PseudoJet const& pseudoJet) { return subClustersBuilt[static_cast<std::size_t>(pseudoJet.user_index())];})));
        
        // TODO set pdgId to something

        indexRecorder.recordParentClusterIndex(currentCaloParticleIndex);
    }

    fjInputs.clear();
}

private:
    SimClusterCollection &clusters_;

    std::vector<fastjet::PseudoJet> fjInputs;

    ClusterParentIndexRecorderT indexRecorder;

    fastjet::JetDefinition const& jetDefinition_;
};


class SimClusterMergerByFastJetConfig : public SimClusterConfig {
public:
SimClusterMergerByFastJetConfig(edm::ProducesCollector& c, std::string tag, const edm::ParameterSet& ps) 
: SimClusterConfig(c, tag), jetClusteringRadius_(ps.getParameter<double>("jetClusteringRadius")),
jetDefinition_(fastjet::antikt_algorithm, jetClusteringRadius_) {}

void fillPSetDescription(edm::ParameterSetDescription& desc) {
    desc.add<double>("jetClusteringRadius", 0.05)->setComment("Distance parameter for clustering algorithm");
}

template <typename ClusterParentIndexRecorderT>
auto getMerger(ClusterParentIndexRecorderT parentIndexRecorder) {
    return SimClusterMergerByFastJet(outputClusters, parentIndexRecorder, jetDefinition_);
}

private:
 double jetClusteringRadius_;
 fastjet::JetDefinition jetDefinition_;
};
};

#endif