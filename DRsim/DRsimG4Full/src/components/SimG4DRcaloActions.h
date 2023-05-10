#ifndef SimG4DRcaloActions_h
#define SimG4DRcaloActions_h 1

#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/MCParticleCollection.h"

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SystemOfUnits.h"

#include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/ISimG4ActionTool.h"

#include "SimG4DRcaloActionInitialization.h"

#include "GridDRcalo.h"
#include "SCEPCALGridDRcalo.h"

class SimG4DRcaloActions : public AlgTool, virtual public ISimG4ActionTool {
public:
  explicit SimG4DRcaloActions(const std::string& type, const std::string& name, const IInterface* parent);
  virtual ~SimG4DRcaloActions();

  virtual StatusCode initialize() final;
  virtual StatusCode finalize() final;

  virtual G4VUserActionInitialization* userActionInitialization() final;

private:
  ServiceHandle<IGeoSvc> m_geoSvc;

  dd4hep::DDSegmentation::GridDRcalo* pSegHCAL;
  dd4hep::DDSegmentation::SCEPCALGridDRcalo* pSegECAL;

  Gaudi::Property<std::string> m_readoutNameHCAL{this, "readoutNameHCAL", "DRcaloSiPMreadout", "readout name of DRcalo"};
  Gaudi::Property<std::string> m_readoutNameECAL{this, "readoutNameECAL", "SCEPCALDRcaloSiPMreadout", "readout name of SCEPCALDRcalo"};

  Gaudi::Property<std::string> m_scintNameHCAL{this, "scintNameHCAL", "DR_Polystyrene", "Name of the scintillators"};
  Gaudi::Property<std::string> m_scintNameECAL{this, "scintNameECAL", "PbWO", "Name of the scintillators"};

  Gaudi::Property<double> m_birksHCAL{this, "birksHCAL", 0.126, "Birk's constant for the scintillators in mm/MeV"};
  Gaudi::Property<double> m_birksECAL{this, "birksECAL", 0.126, "Birk's constant for the scintillators in mm/MeV"};

  Gaudi::Property<double> m_thresHCAL{this, "thresHCAL", 0.0, "Energy threshold to store 3d SimCalorimeterHits in GeV"};
  Gaudi::Property<double> m_thresECAL{this, "thresECAL", 0.0, "Energy threshold to store 3d SimCalorimeterHits in GeV"};
};

#endif
