#include "SimG4DRcaloActionInitialization.h"
#include "SimG4DRcaloSteppingAction.h"
#include "SimG4DRcaloEventAction.h"
#include "CLHEP/Units/SystemOfUnits.h"

namespace drc {
SimG4DRcaloActionInitialization::SimG4DRcaloActionInitialization(): G4VUserActionInitialization() {}

SimG4DRcaloActionInitialization::~SimG4DRcaloActionInitialization() {}

void SimG4DRcaloActionInitialization::setBirksConstantHCAL(const std::string scintName, const double birks) {
  m_scintNameHCAL = scintName;
  m_birksHCAL = birks;
}

void SimG4DRcaloActionInitialization::setBirksConstantECAL(const std::string scintName, const double birks) {
  m_scintNameECAL = scintName;
  m_birksECAL = birks;
}

void SimG4DRcaloActionInitialization::Build() const {

  SimG4DRcaloSteppingAction* steppingAction = new SimG4DRcaloSteppingAction(); // deleted by G4

  steppingAction->setSegmentationHCAL(pSegHCAL);
  steppingAction->setThresholdHCAL(m_thresHCAL);

  steppingAction->setSegmentationECAL(pSegECAL);
  steppingAction->setThresholdECAL(m_thresECAL);

  SetUserAction(steppingAction);

  SimG4DRcaloEventAction* eventAction = new SimG4DRcaloEventAction(); // deleted by G4

  eventAction->setSteppingAction(steppingAction);
  SetUserAction(eventAction);

  G4Material::GetMaterial(m_scintNameHCAL)->GetIonisation()->SetBirksConstant(m_birksHCAL*CLHEP::millimeter/CLHEP::MeV); // makeshift for DD4hep
  G4Material::GetMaterial(m_scintNameECAL)->GetIonisation()->SetBirksConstant(m_birksECAL*CLHEP::millimeter/CLHEP::MeV); // makeshift for DD4hep
}
}
