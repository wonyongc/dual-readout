#include "SimG4SaveDRcaloMCTruth.h"

#include "G4RunManager.hh"

DECLARE_COMPONENT(SimG4SaveDRcaloMCTruth)

SimG4SaveDRcaloMCTruth::SimG4SaveDRcaloMCTruth(const std::string& aType, const std::string& aName, const IInterface* aParent)
: GaudiTool(aType, aName, aParent), m_geantSvc("SimG4Svc", aName) {
  declareInterface<ISimG4SaveOutputTool>(this);
}

SimG4SaveDRcaloMCTruth::~SimG4SaveDRcaloMCTruth() {}

StatusCode SimG4SaveDRcaloMCTruth::initialize() {
  if (GaudiTool::initialize().isFailure())
    return StatusCode::FAILURE;

  if (!m_geantSvc) {
    error() << "Unable to locate Geant Simulation Service" << endmsg;
    return StatusCode::FAILURE;
  }

  auto* runManager = G4RunManager::GetRunManager();
  const auto* eventAction = dynamic_cast<const drc::SimG4DRcaloEventAction*>(runManager->GetUserEventAction());

  if (!eventAction) {
    error() << "Unable to cast to SimG4DRcaloEventAction from G4UserEventAction" << endmsg;
    return StatusCode::FAILURE;
  }

  m_eventAction = const_cast<drc::SimG4DRcaloEventAction*>(eventAction); // HACK!!!

  return StatusCode::SUCCESS;
}

StatusCode SimG4SaveDRcaloMCTruth::finalize() { return GaudiTool::finalize(); }

StatusCode SimG4SaveDRcaloMCTruth::saveOutput(const G4Event&) {

  auto* edepsHCAL = m_eventAction->getEdepsCollectionHCAL();
  auto* edeps3dHCAL = m_eventAction->getEdeps3dCollectionHCAL();

  auto* edepsECALF = m_eventAction->getEdepsCollectionECALF();
  auto* edeps3dECALF = m_eventAction->getEdeps3dCollectionECALF();

  auto* edepsECALR = m_eventAction->getEdepsCollectionECALR();
  auto* edeps3dECALR = m_eventAction->getEdeps3dCollectionECALR();

  auto* edepsECALFCher = m_eventAction->getEdepsCollectionECALFCher();
  auto* edepsECALRCher = m_eventAction->getEdepsCollectionECALRCher();

  auto* leakages = m_eventAction->getLeakagesCollection();

  m_EdepsHCAL.put(edepsHCAL);
  m_Edeps3dHCAL.put(edeps3dHCAL);

  m_EdepsECALF.put(edepsECALF);
  m_Edeps3dECALF.put(edeps3dECALF);

  m_EdepsECALR.put(edepsECALR);
  m_Edeps3dECALR.put(edeps3dECALR);

  m_EdepsECALFCher.put(edepsECALFCher);
  m_EdepsECALRCher.put(edepsECALRCher);

  m_Leakages.put(leakages);



  return StatusCode::SUCCESS;
}
