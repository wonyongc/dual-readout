#ifndef SimG4DRcaloActionInitialization_h
#define SimG4DRcaloActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

#include "GridDRcalo.h"
#include "SCEPCALGridDRcalo.h"

namespace drc {
class SimG4DRcaloActionInitialization : public G4VUserActionInitialization {
public:
  SimG4DRcaloActionInitialization();
  virtual ~SimG4DRcaloActionInitialization();

  virtual void Build() const final;

  void setSegmentationHCAL(dd4hep::DDSegmentation::GridDRcalo* seg) { pSegHCAL = seg; }
  void setThresholdHCAL(const double thres) { m_thresHCAL = thres; }
  void setBirksConstantHCAL(const std::string scintName, const double birks);

  void setSegmentationECAL(dd4hep::DDSegmentation::SCEPCALGridDRcalo* seg) { pSegECAL = seg; }
  void setThresholdECAL(const double thres) { m_thresECAL = thres; }
  void setBirksConstantECAL(const std::string scintName, const double birks);

private:
  dd4hep::DDSegmentation::GridDRcalo* pSegHCAL;
  dd4hep::DDSegmentation::SCEPCALGridDRcalo* pSegECAL;

  std::string m_scintNameHCAL;
  std::string m_scintNameECAL;

  double m_birksHCAL;
  double m_birksECAL;

  double m_thresHCAL;
  double m_thresECAL;
};
}

#endif
