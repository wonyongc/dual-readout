from Gaudi.Configuration import *
from Configurables import ApplicationMgr

from Configurables import k4DataSvc
dataservice = k4DataSvc("EventDataSvc", input="digi.root")

from Configurables import PodioInput
podioinput = PodioInput("PodioInput",
  collections = [
    "DigiCalorimeterHits",
    "DigiWaveforms",
    "RawTimeStructs",
    "RawCalorimeterHits",
    "SimCalorimeterHits",
    "RawWavlenStructs",
    "GenParticles",
    "Leakages"
  ],
  OutputLevel = DEBUG
)

from Configurables import GeoSvc
geoservice = GeoSvc(
    "GeoSvc",
    detectors = [
        'file:/afs/cern.ch/user/w/wochung/private/dual-readout/install/share/compact/DRcalo.xml'
    ]
)

from Configurables import DRcalib2D
calib2d = DRcalib2D("DRcalib2D", OutputLevel=DEBUG, calibPath="calib.csv")

from Configurables import PodioOutput
podiooutput = PodioOutput("PodioOutput", filename = "reco2D.root", OutputLevel = DEBUG)
podiooutput.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg = [
        podioinput,
        calib2d,
        podiooutput
    ],
    EvtSel = 'NONE',
    EvtMax = 10,
    ExtSvc = [dataservice, geoservice]
)
