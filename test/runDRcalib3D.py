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
    "Sim3dCalorimeterHits",
    "RawWavlenStructs",
    "GenParticles",
    "Leakages"
  ]
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

from Configurables import DRcalib3D
calib3d = DRcalib3D("DRcalib3D", OutputLevel=DEBUG, veloFile="velo.root")

from Configurables import PodioOutput
podiooutput = PodioOutput("PodioOutput", filename = "reco3D.root")
podiooutput.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg = [
        podioinput,
        calib2d,
        calib3d,
        podiooutput
    ],
    EvtSel = 'NONE',
    EvtMax = 10,
    ExtSvc = [dataservice, geoservice]
)
