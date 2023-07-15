#include "SimG4DRcaloSteppingAction.h"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VProcess.hh"

#include "CLHEP/Units/SystemOfUnits.h"
#include "DD4hep/DD4hepUnits.h"
#include <bitset>

namespace drc {

SimG4DRcaloSteppingAction::SimG4DRcaloSteppingAction()
: G4UserSteppingAction(), fPrevTowerHCAL(0), fPrevFiberHCAL(0), fPrevTowerECALF(0), fPrevTowerECALR(0),
  fPrevTowerECALFCher(0), fPrevTowerECALRCher(0) {}

SimG4DRcaloSteppingAction::~SimG4DRcaloSteppingAction() {}


void SimG4DRcaloSteppingAction::UserSteppingAction(const G4Step* step) {

  G4StepPoint* presteppoint = step->GetPreStepPoint();
  G4StepPoint* poststeppoint = step->GetPostStepPoint();

  G4Track* track = step->GetTrack();
  G4ParticleDefinition* particle = track->GetDefinition();

  G4VPhysicalVolume* PreStepVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  std::cout<<"Pre-step Volume Name: " << PreStepVolume->GetName() << std::endl;
  std::cout<<"Particle Name, Type, SubType: " << particle->GetParticleName()<< " "<<particle->GetParticleType() << " "<<particle->GetParticleSubType()<< std::endl;


  G4TouchableHandle theTouchable = presteppoint->GetTouchableHandle();

  if ( theTouchable->GetHistoryDepth()<2 ) return; // skip particles in the world or assembly volume

  //historyDepth    world
  //historyDepth-1  experimentalhall (detector)
  //historyDepth-2  towerAssemblyVol, assemblyEnvelopVol

  int volID = theTouchable->GetCopyNumber();

  std::bitset<32> volIdbits(volID);

  int system = (volID) & (32-1);
  int eta = (volID >> 5) & (1024-1);
  int phi = (volID >> 15) & (1024-1);
  int depth = (volID >> 25) & (8-1);

  std::cout << "Stepping Action start: history depth: " << theTouchable->GetHistoryDepth() <<std::endl;
  std::cout << "Stepping Action start: volID bits: " << volIdbits <<std::endl;
  std::cout << "Stepping Action start: system: " << system << " eta: "<<eta<<" phi: "<<phi<<" depth: "<< depth <<std::endl;


  if (system==1) ECALSteppingAction(step, presteppoint, particle);
  else if (system==2) HCALSteppingAction(step, presteppoint, particle);

  // leakage particles
  if (poststeppoint->GetStepStatus() == fWorldBoundary) {
    saveLeakage(track,presteppoint);
    return;
  }
}

void SimG4DRcaloSteppingAction::ECALSteppingAction(const G4Step* step, G4StepPoint* presteppoint, G4ParticleDefinition* particle) {
  G4Track* track = step->GetTrack();
  G4TouchableHandle theTouchable = presteppoint->GetTouchableHandle();

  float edep=step->GetTotalEnergyDeposit()*CLHEP::MeV/CLHEP::GeV;
//  int towerNum32=theTouchable->GetCopyNumber(theTouchable->GetHistoryDepth()-2);
  int towerNum32=theTouchable->GetCopyNumber();

  std::bitset<32> towerNum32bits(towerNum32);

  int system = (towerNum32) & (32-1);
  int eta = (towerNum32 >> 5) & (1024-1);
  int phi = (towerNum32 >> 15) & (1024-1);
  int depth = (towerNum32 >> 25) & (8-1);

  std::cout << "  ECALSteppingAction start: history depth: " << theTouchable->GetHistoryDepth() <<std::endl;
  std::cout << "  ECALSteppingAction start: towerNum32 bits: " << towerNum32bits <<std::endl;
  std::cout<<"  ECAL crystal eta: "<<eta<<" phi: "<<phi<<" depth: "<< depth <<std::endl;
  std::cout<<"  ECAL crystal edep: "<<edep<<std::endl;

  auto towerNum64=pSegECAL->convertFirst32to64(towerNum32);

  if (depth==1) {
    if (edep>m_thresECAL) {
      auto simEdep3dECALF=m_Edeps3dECALF->create();

      simEdep3dECALF.setCellID(static_cast<unsigned long long>(towerNum64));
      simEdep3dECALF.setEnergy(edep);

      auto &pos=presteppoint->GetPosition();

      std::cout<<"    ECAL step position x: "<<pos.x()*CLHEP::millimeter<<" y: "<<pos.y()*CLHEP::millimeter<<" z: "<< pos.z()*CLHEP::millimeter <<std::endl;

      simEdep3dECALF.setPosition({static_cast<float>(pos.x()*CLHEP::millimeter),
                             static_cast<float>(pos.y()*CLHEP::millimeter),
                             static_cast<float>(pos.z()*CLHEP::millimeter)});
    }
    accumulateECALF(fPrevTowerECALF, towerNum64, edep);
  } else if (depth==2) {
    if (edep>m_thresECAL) {
      auto simEdep3dECALR=m_Edeps3dECALR->create();

      simEdep3dECALR.setCellID(static_cast<unsigned long long>(towerNum64));
      simEdep3dECALR.setEnergy(edep);

      auto &pos=presteppoint->GetPosition();
      std::cout<<"    ECAL step position x: "<<pos.x()*CLHEP::millimeter<<" y: "<<pos.y()*CLHEP::millimeter<<" z: "<< pos.z()*CLHEP::millimeter <<std::endl;

      simEdep3dECALR.setPosition({static_cast<float>(pos.x()*CLHEP::millimeter),
                             static_cast<float>(pos.y()*CLHEP::millimeter),
                             static_cast<float>(pos.z()*CLHEP::millimeter)});
    }
    accumulateECALR(fPrevTowerECALR, towerNum64, edep);
  }

  if ((track->GetCurrentStepNumber()==1)&&particle==G4OpticalPhoton::OpticalPhotonDefinition()) {

    G4String processName=track->GetCreatorProcess()->GetProcessName();

    if (processName=="Cerenkov") {
      //kill very long or short wavelengths
      float photWL=1239.84187/(track->GetTotalEnergy()*CLHEP::eV/CLHEP::GeV);

      if (photWL>1000||photWL<300) {
        std::cout<<"Cerenkov photon killed with WL: "<<photWL<<std::endl;
        track->SetTrackStatus(fKillTrackAndSecondaries);
      } else {
        if (depth==1) {
          accumulateECALFCher(fPrevTowerECALFCher, towerNum64);
        } else if (depth==2) {
          accumulateECALRCher(fPrevTowerECALRCher, towerNum64);
        }

        //do not propagate the photon
//        theTrack->SetTrackStatus(fKillTrackAndSecondaries);
      }

    };


    return;
  }
}

void SimG4DRcaloSteppingAction::HCALSteppingAction(const G4Step* step, G4StepPoint* presteppoint, G4ParticleDefinition* particle) {

  if ( particle == G4OpticalPhoton::OpticalPhotonDefinition() ) return;
  G4TouchableHandle theTouchable = presteppoint->GetTouchableHandle();

  float edep = step->GetTotalEnergyDeposit()*CLHEP::MeV/CLHEP::GeV;
  int towerNum32 = theTouchable->GetCopyNumber( theTouchable->GetHistoryDepth()-2 );

  std::bitset<32> towerNum32bits(towerNum32);

  int system = (towerNum32) & (32-1);
  int eta = (towerNum32 >> 5) & (1024-1);
  int phi = (towerNum32 >> 15) & (1024-1);
  int depth = (towerNum32 >> 25) & (8-1);

  std::cout << "  HCALSteppingAction start: history depth: " << theTouchable->GetHistoryDepth() <<std::endl;
  std::cout << "  HCALSteppingAction start: towerNum32 bits: " << towerNum32bits <<std::endl;
  std::cout << "  HCAL tower eta: " << eta << " phi: " << phi << " depth: " << depth << std::endl;
  std::cout << "  HCAL crystal edep: "<<edep<<std::endl;

  auto towerNum64 = pSegHCAL->convertFirst32to64( towerNum32 );

  if (edep > m_thresHCAL) {
    auto simEdep3dHCAL = m_Edeps3dHCAL->create();

    simEdep3dHCAL.setCellID( static_cast<unsigned long long>(towerNum64) );
    simEdep3dHCAL.setEnergy(edep);

    auto& pos = presteppoint->GetPosition();
    simEdep3dHCAL.setPosition( { static_cast<float>(pos.x()*CLHEP::millimeter),
                             static_cast<float>(pos.y()*CLHEP::millimeter),
                             static_cast<float>(pos.z()*CLHEP::millimeter) } );
  }
  accumulateHCAL(fPrevTowerHCAL,towerNum64,edep);
  return;
}

void SimG4DRcaloSteppingAction::accumulateECALF(unsigned int &prev, dd4hep::DDSegmentation::CellID& id64, float edep) {
  std::cout<<"      Accumulate ECAL F: "<<std::endl;
  std::cout<<"      Accumulate ECAL F edep: "<<edep<<std::endl;

  // search for the element
  bool found = false;
  edm4hep::SimCalorimeterHit* thePtr = nullptr;

  if ( m_EdepsECALF->size() > prev ) { // check previous element
    auto element = m_EdepsECALF->at(prev);
    if ( checkId(element, id64) ) {
      thePtr = &element;
      found = true;
    }
  }

  if (!found) { // fall back to loop
    std::cout<<"      Accumulate ECAL F loop: "<<std::endl;

    for (unsigned int iElement = 0; iElement<m_EdepsECALF->size(); iElement++) {
      auto element = m_EdepsECALF->at(iElement);
      if ( checkId(element, id64) ) {
        found = true;
        prev = iElement;
        thePtr = &element;

        break;
      }
    }
  }

  if (!found) { // create
    std::cout<<"      Accumulate ECAL F create: "<<std::endl;

    auto simEdepECALF = m_EdepsECALF->create();
    simEdepECALF.setCellID( static_cast<unsigned long long>(id64) );
    simEdepECALF.setEnergy(0.); // added later

    auto pos = pSegECAL->position(id64);
    simEdepECALF.setPosition( { static_cast<float>(pos.x()*CLHEP::millimeter/dd4hep::millimeter),
                               static_cast<float>(pos.y()*CLHEP::millimeter/dd4hep::millimeter),
                               static_cast<float>(pos.z()*CLHEP::millimeter/dd4hep::millimeter) } );

    prev = m_EdepsECALF->size();
    thePtr = &simEdepECALF;
  }

  auto edepPrev = thePtr->getEnergy();
  thePtr->setEnergy( edepPrev + edep );
}
void SimG4DRcaloSteppingAction::accumulateECALR(unsigned int &prev, dd4hep::DDSegmentation::CellID& id64, float edep) {
  std::cout<<"      Accumulate ECAL R: "<<std::endl;
  std::cout<<"      Accumulate ECAL R edep: "<<edep<<std::endl;

  // search for the element
  bool found = false;
  edm4hep::SimCalorimeterHit* thePtr = nullptr;

  if ( m_EdepsECALR->size() > prev ) { // check previous element
    auto element = m_EdepsECALR->at(prev);
    if ( checkId(element, id64) ) {
      thePtr = &element;
      found = true;
    }
  }

  if (!found) { // fall back to loop
    std::cout<<"      Accumulate ECAL R loop: "<<std::endl;

    for (unsigned int iElement = 0; iElement<m_EdepsECALR->size(); iElement++) {
      auto element = m_EdepsECALR->at(iElement);
      if ( checkId(element, id64) ) {
        found = true;
        prev = iElement;
        thePtr = &element;

        break;
      }
    }
  }

  if (!found) { // create
    std::cout<<"      Accumulate ECAL R create: "<<std::endl;

    auto simEdepECALR = m_EdepsECALR->create();
    simEdepECALR.setCellID( static_cast<unsigned long long>(id64) );
    simEdepECALR.setEnergy(0.); // added later

    auto pos = pSegECAL->position(id64);
    simEdepECALR.setPosition( { static_cast<float>(pos.x()*CLHEP::millimeter/dd4hep::millimeter),
                               static_cast<float>(pos.y()*CLHEP::millimeter/dd4hep::millimeter),
                               static_cast<float>(pos.z()*CLHEP::millimeter/dd4hep::millimeter) } );

    prev = m_EdepsECALR->size();
    thePtr = &simEdepECALR;
  }

  auto edepPrev = thePtr->getEnergy();
  thePtr->setEnergy( edepPrev + edep );
}
void SimG4DRcaloSteppingAction::accumulateECALFCher(unsigned int &prev, dd4hep::DDSegmentation::CellID& id64) {
  std::cout<<"      Accumulate ECAL F Cher: "<<std::endl;

  // search for the element
  bool found = false;
  edm4hep::SimCalorimeterHit* thePtr = nullptr;

  if ( m_EdepsECALFCher->size() > prev ) { // check previous element
    auto element = m_EdepsECALFCher->at(prev);
    if ( checkId(element, id64) ) {
      thePtr = &element;
      found = true;
    }
  }

  if (!found) { // fall back to loop
    std::cout<<"      Accumulate ECAL F Cher loop: "<<std::endl;

    for (unsigned int iElement = 0; iElement<m_EdepsECALFCher->size(); iElement++) {
      auto element = m_EdepsECALFCher->at(iElement);
      if ( checkId(element, id64) ) {
        found = true;
        prev = iElement;
        thePtr = &element;

        break;
      }
    }
  }

  if (!found) { // create
    std::cout<<"      Accumulate ECAL F Cher create: "<<std::endl;

    auto simEdepECALFCher = m_EdepsECALFCher->create();
    simEdepECALFCher.setCellID( static_cast<unsigned long long>(id64) );
    simEdepECALFCher.setEnergy(0.); // added later

    auto pos = pSegECAL->position(id64);
    simEdepECALFCher.setPosition( { static_cast<float>(pos.x()*CLHEP::millimeter/dd4hep::millimeter),
                                static_cast<float>(pos.y()*CLHEP::millimeter/dd4hep::millimeter),
                                static_cast<float>(pos.z()*CLHEP::millimeter/dd4hep::millimeter) } );

    prev = m_EdepsECALFCher->size();
    thePtr = &simEdepECALFCher;
  }

  auto edepPrev = thePtr->getEnergy();
  thePtr->setEnergy( edepPrev + 1 );
}
void SimG4DRcaloSteppingAction::accumulateECALRCher(unsigned int &prev, dd4hep::DDSegmentation::CellID& id64) {
  std::cout<<"      Accumulate ECAL R Cher: "<<std::endl;

  // search for the element
  bool found = false;
  edm4hep::SimCalorimeterHit* thePtr = nullptr;

  if ( m_EdepsECALRCher->size() > prev ) { // check previous element
    auto element = m_EdepsECALRCher->at(prev);
    if ( checkId(element, id64) ) {
      thePtr = &element;
      found = true;
    }
  }

  if (!found) { // fall back to loop
    std::cout<<"      Accumulate ECAL R Cher loop: "<<std::endl;

    for (unsigned int iElement = 0; iElement<m_EdepsECALRCher->size(); iElement++) {
      auto element = m_EdepsECALRCher->at(iElement);
      if ( checkId(element, id64) ) {
        found = true;
        prev = iElement;
        thePtr = &element;

        break;
      }
    }
  }

  if (!found) { // create
    std::cout<<"      Accumulate ECAL R Cher create: "<<std::endl;

    auto simEdepECALRCher = m_EdepsECALRCher->create();
    simEdepECALRCher.setCellID( static_cast<unsigned long long>(id64) );
    simEdepECALRCher.setEnergy(0.); // added later

    auto pos = pSegECAL->position(id64);
    simEdepECALRCher.setPosition( { static_cast<float>(pos.x()*CLHEP::millimeter/dd4hep::millimeter),
                                    static_cast<float>(pos.y()*CLHEP::millimeter/dd4hep::millimeter),
                                    static_cast<float>(pos.z()*CLHEP::millimeter/dd4hep::millimeter) } );

    prev = m_EdepsECALRCher->size();
    thePtr = &simEdepECALRCher;
  }

  auto edepPrev = thePtr->getEnergy();
  thePtr->setEnergy( edepPrev + 1 );
}

void SimG4DRcaloSteppingAction::accumulateHCAL(unsigned int &prev, dd4hep::DDSegmentation::CellID& id64, float edep) {
  std::cout<<"      Accumulate HCAL: "<<std::endl;
  std::cout<<"      Accumulate HCAL edep: "<<edep<<std::endl;

  // search for the element
  bool found = false;
  edm4hep::SimCalorimeterHit* thePtr = nullptr;

  if ( m_EdepsHCAL->size() > prev ) { // check previous element
    auto element = m_EdepsHCAL->at(prev);
    if ( checkId(element, id64) ) {
      thePtr = &element;
      found = true;
    }
  }

  if (!found) { // fall back to loop
    std::cout<<"      Accumulate HCAL loop: "<<std::endl;

    for (unsigned int iElement = 0; iElement<m_EdepsHCAL->size(); iElement++) {
      auto element = m_EdepsHCAL->at(iElement);
      if ( checkId(element, id64) ) {
        found = true;
        prev = iElement;
        thePtr = &element;

        break;
      }
    }
  }

  if (!found) { // create
    std::cout<<"      Accumulate HCAL create: "<<std::endl;

    auto simEdepHCAL = m_EdepsHCAL->create();
    simEdepHCAL.setCellID( static_cast<unsigned long long>(id64) );
    simEdepHCAL.setEnergy(0.); // added later

    auto pos = pSegHCAL->position(id64);
    simEdepHCAL.setPosition( { static_cast<float>(pos.x()*CLHEP::millimeter/dd4hep::millimeter),
                               static_cast<float>(pos.y()*CLHEP::millimeter/dd4hep::millimeter),
                               static_cast<float>(pos.z()*CLHEP::millimeter/dd4hep::millimeter) } );

    prev = m_EdepsHCAL->size();
    thePtr = &simEdepHCAL;
  }

  auto edepPrev = thePtr->getEnergy();
  thePtr->setEnergy( edepPrev + edep );
}

bool SimG4DRcaloSteppingAction::checkId(edm4hep::SimCalorimeterHit edep, dd4hep::DDSegmentation::CellID& id64) {
  return ( edep.getCellID()==static_cast<unsigned long long>(id64) );
}

void SimG4DRcaloSteppingAction::saveLeakage(G4Track* track, G4StepPoint* presteppoint) {
  std::cout<<"Leakage particle: "<<std::endl;
  std::cout<<"Leakage position x: "<<presteppoint->GetPosition().x()*CLHEP::millimeter<<" y: "<<presteppoint->GetPosition().y()*CLHEP::millimeter<<" z: "<< presteppoint->GetPosition().z()*CLHEP::millimeter <<std::endl;

  auto leakage = m_Leakages->create();
  leakage.setPDG( track->GetDefinition()->GetPDGEncoding() );
  leakage.setGeneratorStatus(1); // leakages naturally belong to final states
  leakage.setCharge( track->GetDefinition()->GetPDGCharge() );
  leakage.setMass( track->GetDefinition()->GetPDGMass()*CLHEP::MeV/CLHEP::GeV );
  leakage.setMomentum( { static_cast<float>(track->GetMomentum().x()*CLHEP::MeV/CLHEP::GeV),
                         static_cast<float>(track->GetMomentum().y()*CLHEP::MeV/CLHEP::GeV),
                         static_cast<float>(track->GetMomentum().z()*CLHEP::MeV/CLHEP::GeV) } );
  leakage.setVertex( { static_cast<float>(presteppoint->GetPosition().x()*CLHEP::millimeter),
                       static_cast<float>(presteppoint->GetPosition().y()*CLHEP::millimeter),
                       static_cast<float>(presteppoint->GetPosition().z()*CLHEP::millimeter) } );
}

} // namespace drc
