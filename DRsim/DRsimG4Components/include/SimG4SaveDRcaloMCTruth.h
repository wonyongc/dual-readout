#ifndef SimG4SaveDRcaloMCTruth_h
#define SimG4SaveDRcaloMCTruth_h 1

// Data model
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#include "GaudiAlg/GaudiTool.h"
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ISimG4Svc.h"
#include "k4Interface/ISimG4SaveOutputTool.h"

#include "SimG4DRcaloEventAction.h"

class IGeoSvc;

class SimG4SaveDRcaloMCTruth : public GaudiTool, virtual public ISimG4SaveOutputTool {
public:
  explicit SimG4SaveDRcaloMCTruth(const std::string& aType, const std::string& aName, const IInterface* aParent);
  virtual ~SimG4SaveDRcaloMCTruth();

  virtual StatusCode initialize();
  virtual StatusCode finalize();

  virtual StatusCode saveOutput(const G4Event& aEvent) final;

private:
  DataHandle<edm4hep::SimCalorimeterHitCollection> m_EdepsHCAL{"SimCalorimeterHitsHCAL", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::SimCalorimeterHitCollection> m_Edeps3dHCAL{"Sim3dCalorimeterHitsHCAL", Gaudi::DataHandle::Writer, this};

  DataHandle<edm4hep::SimCalorimeterHitCollection> m_EdepsECALF{"SimCalorimeterHitsECALF", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::SimCalorimeterHitCollection> m_Edeps3dECALF{"Sim3dCalorimeterHitsECALF", Gaudi::DataHandle::Writer, this};

  DataHandle<edm4hep::SimCalorimeterHitCollection> m_EdepsECALR{"SimCalorimeterHitsECALR", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::SimCalorimeterHitCollection> m_Edeps3dECALR{"Sim3dCalorimeterHitsECALR", Gaudi::DataHandle::Writer, this};

  DataHandle<edm4hep::SimCalorimeterHitCollection> m_EdepsECALFCher{"SimCalorimeterHitsECALFCher", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::SimCalorimeterHitCollection> m_EdepsECALRCher{"SimCalorimeterHitsECALRCher", Gaudi::DataHandle::Writer, this};

  DataHandle<edm4hep::MCParticleCollection> m_Leakages{"Leakages", Gaudi::DataHandle::Writer, this};

  drc::SimG4DRcaloEventAction* m_eventAction;

  ServiceHandle<ISimG4Svc> m_geantSvc;
};

#endif
