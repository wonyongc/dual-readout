#include "SimG4DRcaloEventAction.h"

#include "G4RunManager.hh"

namespace drc {
SimG4DRcaloEventAction::SimG4DRcaloEventAction(): G4UserEventAction() {}

SimG4DRcaloEventAction::~SimG4DRcaloEventAction() {}

void SimG4DRcaloEventAction::BeginOfEventAction(const G4Event*) {

  m_EdepsHCAL = new edm4hep::SimCalorimeterHitCollection();
  m_Edeps3dHCAL = new edm4hep::SimCalorimeterHitCollection();

  m_EdepsECALF = new edm4hep::SimCalorimeterHitCollection();
  m_Edeps3dECALF = new edm4hep::SimCalorimeterHitCollection();

  m_EdepsECALR = new edm4hep::SimCalorimeterHitCollection();
  m_Edeps3dECALR = new edm4hep::SimCalorimeterHitCollection();

  m_EdepsECALFCher = new edm4hep::SimCalorimeterHitCollection();
  m_EdepsECALRCher = new edm4hep::SimCalorimeterHitCollection();

  m_Leakages = new edm4hep::MCParticleCollection();


  pSteppingAction->setEdepsCollectionHCAL(m_EdepsHCAL);
  pSteppingAction->setEdeps3dCollectionHCAL(m_Edeps3dHCAL);

  pSteppingAction->setEdepsCollectionECALF(m_EdepsECALF);
  pSteppingAction->setEdeps3dCollectionECALF(m_Edeps3dECALF);

  pSteppingAction->setEdepsCollectionECALR(m_EdepsECALR);
  pSteppingAction->setEdeps3dCollectionECALR(m_Edeps3dECALR);

  pSteppingAction->setEdepsCollectionECALFCher(m_EdepsECALFCher);
  pSteppingAction->setEdepsCollectionECALRCher(m_EdepsECALRCher);

  pSteppingAction->setLeakagesCollection(m_Leakages);

  return;
}

void SimG4DRcaloEventAction::EndOfEventAction(const G4Event*) {}
} // namespace drc
