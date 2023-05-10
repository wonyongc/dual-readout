#ifndef SimG4DRcaloEventAction_h
#define SimG4DRcaloEventAction_h 1

#include "G4UserEventAction.hh"
#include "G4Event.hh"

#include "SimG4DRcaloSteppingAction.h"

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

namespace drc {
class SimG4DRcaloEventAction : public G4UserEventAction {
    
public:
  SimG4DRcaloEventAction();
  virtual ~SimG4DRcaloEventAction();

  virtual void BeginOfEventAction(const G4Event*) final;
  virtual void EndOfEventAction(const G4Event*) final;

  void setSteppingAction(SimG4DRcaloSteppingAction* steppingAction) { pSteppingAction = steppingAction; }

  edm4hep::SimCalorimeterHitCollection* getEdepsCollectionHCAL() { return m_EdepsHCAL; }
  edm4hep::SimCalorimeterHitCollection* getEdeps3dCollectionHCAL() { return m_Edeps3dHCAL; }

  edm4hep::SimCalorimeterHitCollection* getEdepsCollectionECALF() { return m_EdepsECALF; }
  edm4hep::SimCalorimeterHitCollection* getEdeps3dCollectionECALF() { return m_Edeps3dECALF; }

  edm4hep::SimCalorimeterHitCollection* getEdepsCollectionECALR() { return m_EdepsECALR; }
  edm4hep::SimCalorimeterHitCollection* getEdeps3dCollectionECALR() { return m_Edeps3dECALR; }

  edm4hep::SimCalorimeterHitCollection* getEdepsCollectionECALFCher() { return m_EdepsECALFCher; }
  edm4hep::SimCalorimeterHitCollection* getEdepsCollectionECALRCher() { return m_EdepsECALRCher; }

  edm4hep::MCParticleCollection* getLeakagesCollection() { return m_Leakages; }

private:
  SimG4DRcaloSteppingAction* pSteppingAction;

  // ownership of collections transferred to DataWrapper<T>
  edm4hep::SimCalorimeterHitCollection* m_EdepsHCAL;
  edm4hep::SimCalorimeterHitCollection* m_Edeps3dHCAL;

  edm4hep::SimCalorimeterHitCollection* m_EdepsECALF;
  edm4hep::SimCalorimeterHitCollection* m_Edeps3dECALF;

  edm4hep::SimCalorimeterHitCollection* m_EdepsECALR;
  edm4hep::SimCalorimeterHitCollection* m_Edeps3dECALR;

  edm4hep::SimCalorimeterHitCollection* m_EdepsECALFCher;
  edm4hep::SimCalorimeterHitCollection* m_EdepsECALRCher;

  edm4hep::MCParticleCollection* m_Leakages;

};

}

#endif
