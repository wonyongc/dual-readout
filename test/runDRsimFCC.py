### | read events from a HepMC file | convert `HepMC::GenEvent` to EDM | geometry taken from GDML file (no sensitive detectors!) | FTFP_BERT physics list | empty action initialisation list | write the EDM output to ROOT file using PODIO |


from Gaudi.Configuration import *
from GaudiKernel import SystemOfUnits as units
from GaudiKernel import PhysicalConstants as constants


from Configurables import ApplicationMgr
ApplicationMgr().EvtSel = 'NONE'
ApplicationMgr().EvtMax = 100
ApplicationMgr().OutputLevel = INFO
ApplicationMgr().ExtSvc = ['RndmGenSvc']


from Configurables import FCCDataSvc
## Data service
podioevent = FCCDataSvc("EventDataSvc")
podioevent.input = "/home/wonyongc/src/hep/dual-readout/test/hepmc/events_185677613.root"
ApplicationMgr().ExtSvc += [podioevent]



from Configurables import MomentumRangeParticleGun
guntool = MomentumRangeParticleGun()
guntool.ThetaMin = 0
guntool.ThetaMax = 2 * constants.pi
guntool.PdgCodes = [11]


from Configurables import GenAlg
gen = GenAlg()
gen.SignalProvider=guntool
gen.hepmc.Path = "hepmc"
ApplicationMgr().TopAlg += [gen]


from Configurables import HepMCToEDMConverter
## reads an HepMC::GenEvent from the data service and writes a collection of EDM Particles
hepmc_converter = HepMCToEDMConverter("Converter")
hepmc_converter.hepmc.Path="hepmc"
hepmc_converter.GenParticles.Path="allGenParticles"
hepmc_converter.OutputLevel = DEBUG
ApplicationMgr().TopAlg += [hepmc_converter]


from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")

geoservice.detectors = detectors = [
    # 'file:/afs/cern.ch/user/w/wochung/private/dual-readout/install/test/SCEPCAL.xml'
    'file:/home/wonyongc/src/hep/dual-readout/test/SCEPCAL.xml'
]
geoservice.OutputLevel = INFO
ApplicationMgr().ExtSvc += [geoservice]


from Configurables import SimG4Svc, SimG4FastSimPhysicsList, SimG4FastSimOpFiberRegion, SimG4OpticalPhysicsList
regionTool = SimG4FastSimOpFiberRegion("fastfiber")
opticalPhysicsTool = SimG4OpticalPhysicsList("opticalPhysics", fullphysics="SimG4FtfpBert")
physicslistTool = SimG4FastSimPhysicsList("Physics", fullphysics=opticalPhysicsTool)

from Configurables import SimG4DRcaloActions
actionTool = SimG4DRcaloActions("SimG4DRcaloActions")

geantservice = SimG4Svc("SimG4Svc",
                        physicslist = physicslistTool,
                        regions = ["SimG4FastSimOpFiberRegion/fastfiber"],
                        actions = actionTool
                        )
ApplicationMgr().ExtSvc += [geantservice]

from Configurables import SimG4SaveDRcaloHits, SimG4SaveDRcaloMCTruth
saveDRcaloTool = SimG4SaveDRcaloHits("saveDRcaloTool", readoutNames = ["DRcaloSiPMreadout"])
saveMCTruthTool = SimG4SaveDRcaloMCTruth("saveMCTruthTool") # need SimG4DRcaloActions


from Configurables import SimG4Alg
geantsim = SimG4Alg("SimG4Alg")

geantsim.outputs = [
    "SimG4SaveDRcaloHits/saveDRcaloTool",
    "SimG4SaveDRcaloMCTruth/saveMCTruthTool"
]

from Configurables import SimG4PrimariesFromEdmTool
geantsim.eventProvider = SimG4PrimariesFromEdmTool("EdmConverter")
geantsim.eventProvider.GenParticles.Path = "allGenParticles"
ApplicationMgr().TopAlg += [geantsim]

# PODIO algorithm
from Configurables import PodioOutput
podiooutput = PodioOutput()
podiooutput.outputCommands = ["keep *"]
podiooutput.filename = "out_geant_fullsim.root"
ApplicationMgr().TopAlg += [podiooutput]

