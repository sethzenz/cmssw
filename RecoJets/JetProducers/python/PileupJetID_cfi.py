import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.PileupJetIDParams_cfi import * 

stdalgos_4x = cms.VPSet(full,   cutbased,PhilV1)
stdalgos_5x = cms.VPSet(full_53x,cutbased,PhilV1,full_53x_met)

chsalgos_4x = cms.VPSet(full,   cutbased)
chsalgos_5x = cms.VPSet(full_53x_chs,cutbased)
#chsalgos = chsalgos_53x

import os
try:
    cmssw_version = os.environ["CMSSW_VERSION"].replace("CMSSW_","")
except:
    cmssw_version = "5_X"

if cmssw_version.startswith("4"):
    stdalgos    = stdalgos_4x
    chsalgos    = chsalgos_4x
else:
    stdalgos    = stdalgos_5x
    chsalgos    = chsalgos_5x

pileupJetIdProducer = cms.EDProducer('PileupJetIdProducer',
                         produceJetIds = cms.bool(True),
                         jetids = cms.InputTag(""),
                         runMvas = cms.bool(True),
                         jets = cms.InputTag("selectedPatJetsPFlow"),
                         vertexes = cms.InputTag("offlinePrimaryVertices"),
                         algos = cms.VPSet(stdalgos),
                                     
                         rho     = cms.InputTag("fixedGridRhoFastjetAll"),
                         jec     = cms.string("AK4PF"),
                         applyJec = cms.bool(False),
                         inputIsCorrected = cms.bool(True),                                     
                         residualsFromTxt = cms.bool(False),
                         residualsTxt     = cms.FileInPath("RecoJets/JetProducers/data/dummy.txt"),
)

pileupJetIdProducerChs = cms.EDProducer('PileupJetIdProducer',
                         produceJetIds = cms.bool(True),
                         jetids = cms.InputTag(""),
                         runMvas = cms.bool(True),
                         jets = cms.InputTag("selectedPatJetsPFlow"),
                         vertexes = cms.InputTag("offlinePrimaryVertices"),
                         algos = cms.VPSet(chsalgos),
                                        
                         rho     = cms.InputTag("fixedGridRhoFastjetAll"),
                         jec     = cms.string("AK4PFchs"),
                         applyJec = cms.bool(False),
                         inputIsCorrected = cms.bool(True),
                         residualsFromTxt = cms.bool(False),
                         residualsTxt     = cms.FileInPath("RecoJets/JetProducers/data/dummy.txt"),

)

