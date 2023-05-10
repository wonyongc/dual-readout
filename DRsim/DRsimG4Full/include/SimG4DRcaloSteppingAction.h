#ifndef SimG4DRcaloSteppingAction_h
#define SimG4DRcaloSteppingAction_h 1

#include "GridDRcalo.h"
#include "SCEPCALGridDRcalo.h"

#include "G4UserSteppingAction.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"

// Data model
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

namespace drc {
class SimG4DRcaloSteppingAction : public G4UserSteppingAction {
public:
  SimG4DRcaloSteppingAction();
  virtual ~SimG4DRcaloSteppingAction();

  virtual void UserSteppingAction(const G4Step*);

  void setLeakagesCollection(edm4hep::MCParticleCollection* data) { m_Leakages = data; }

  void setSegmentationHCAL(dd4hep::DDSegmentation::GridDRcalo* seg) { pSegHCAL = seg; }

  void setEdepsCollectionHCAL(edm4hep::SimCalorimeterHitCollection* data) { m_EdepsHCAL = data; }
  void setEdeps3dCollectionHCAL(edm4hep::SimCalorimeterHitCollection* data) { m_Edeps3dHCAL = data; }

  void setThresholdHCAL(const double thres) { m_thresHCAL = thres; }


  void setSegmentationECAL(dd4hep::DDSegmentation::SCEPCALGridDRcalo* seg) { pSegECAL = seg; }

  void setEdepsCollectionECALF(edm4hep::SimCalorimeterHitCollection* data) { m_EdepsECALF = data; }
  void setEdeps3dCollectionECALF(edm4hep::SimCalorimeterHitCollection* data) { m_Edeps3dECALF = data; }

  void setEdepsCollectionECALR(edm4hep::SimCalorimeterHitCollection* data) { m_EdepsECALR = data; }
  void setEdeps3dCollectionECALR(edm4hep::SimCalorimeterHitCollection* data) { m_Edeps3dECALR = data; }

  void setEdepsCollectionECALFCher(edm4hep::SimCalorimeterHitCollection* data) { m_EdepsECALFCher = data; }
  void setEdepsCollectionECALRCher(edm4hep::SimCalorimeterHitCollection* data) { m_EdepsECALRCher = data; }


  void setThresholdECAL(const double thres) { m_thresECAL = thres; }

private:
  void ECALSteppingAction(const G4Step* step, G4StepPoint* presteppoint, G4ParticleDefinition* particle);
  void HCALSteppingAction(const G4Step* step, G4StepPoint* presteppoint, G4ParticleDefinition* particle);

  void accumulateHCAL(unsigned int &prev, dd4hep::DDSegmentation::CellID& id64, float edep);

  void accumulateECALF(unsigned int &prev, dd4hep::DDSegmentation::CellID& id64, float edep);
  void accumulateECALR(unsigned int &prev, dd4hep::DDSegmentation::CellID& id64, float edep);

  void accumulateECALFCher(unsigned int &prev, dd4hep::DDSegmentation::CellID& id64);
  void accumulateECALRCher(unsigned int &prev, dd4hep::DDSegmentation::CellID& id64);

  bool checkId(edm4hep::SimCalorimeterHit edep, dd4hep::DDSegmentation::CellID& id64);

  void saveLeakage(G4Track* track, G4StepPoint* pre);

  unsigned int fPrevTowerECALF;
  unsigned int fPrevTowerECALR;

  unsigned int fPrevTowerECALFCher;
  unsigned int fPrevTowerECALRCher;

  unsigned int fPrevTowerHCAL;
  unsigned int fPrevFiberHCAL;

  dd4hep::DDSegmentation::GridDRcalo* pSegHCAL;
  dd4hep::DDSegmentation::SCEPCALGridDRcalo* pSegECAL;

  // collections owned by SimG4DRcaloEventAction
  edm4hep::MCParticleCollection* m_Leakages;

  edm4hep::SimCalorimeterHitCollection* m_EdepsHCAL;
  edm4hep::SimCalorimeterHitCollection* m_Edeps3dHCAL;

  edm4hep::SimCalorimeterHitCollection* m_EdepsECALF;
  edm4hep::SimCalorimeterHitCollection* m_Edeps3dECALF;

  edm4hep::SimCalorimeterHitCollection* m_EdepsECALR;
  edm4hep::SimCalorimeterHitCollection* m_Edeps3dECALR;

  edm4hep::SimCalorimeterHitCollection* m_EdepsECALFCher;
  edm4hep::SimCalorimeterHitCollection* m_EdepsECALRCher;

  double m_thresHCAL;
  double m_thresECAL;

};
}

#endif
