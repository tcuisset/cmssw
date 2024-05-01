import FWCore.ParameterSet.Config as cms

# This modifier is for running TICL v5.

ticl_v5 =  cms.Modifier()

# Modifier for using Mustache for superclustering in TICLv5 instead of the new superclustering DNN
ticl_v5_mustache = cms.Modifier()
