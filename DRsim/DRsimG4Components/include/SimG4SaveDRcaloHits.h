#ifndef SimG4SaveDRcaloHits_h
#define SimG4SaveDRcaloHits_h 1

#include "DRcaloSiPMHit.h"

// Data model
#include "edm4hep/RawCalorimeterHitCollection.h"
#include "edm4hep/SparseVectorCollection.h"

#include "GaudiAlg/GaudiTool.h"
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ISimG4SaveOutputTool.h"
#include "k4Interface/IGeoSvc.h"

class IGeoSvc;

class SimG4SaveDRcaloHits : public GaudiTool, virtual public ISimG4SaveOutputTool {
public:
  explicit SimG4SaveDRcaloHits(const std::string& aType, const std::string& aName, const IInterface* aParent);
  virtual ~SimG4SaveDRcaloHits();

  virtual StatusCode initialize();
  virtual StatusCode finalize();

  virtual StatusCode saveOutput(const G4Event& aEvent) final;

private:
  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;

  Gaudi::Property<std::vector<std::string>> m_readoutNames{this, "readoutNames", {"DRcaloSiPMreadout", "SCEPCALDRcaloSiPMreadout"},
                                                           "Name of the readouts (hits collections) to save"};

  DataHandle<edm4hep::RawCalorimeterHitCollection> mRawCaloHitsHCAL{"RawCalorimeterHitsHCAL", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::SparseVectorCollection> mTimeStructHCAL{"RawTimeStructsHCAL", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::SparseVectorCollection> mWavlenStructHCAL{"RawWavlenStructsHCAL", Gaudi::DataHandle::Writer, this};

  DataHandle<edm4hep::RawCalorimeterHitCollection> mRawCaloHitsECAL{"RawCalorimeterHitsECAL", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::SparseVectorCollection> mTimeStructECAL{"RawTimeStructsECAL", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::SparseVectorCollection> mWavlenStructECAL{"RawWavlenStructsECAL", Gaudi::DataHandle::Writer, this};

};

#endif
