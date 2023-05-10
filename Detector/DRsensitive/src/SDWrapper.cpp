#include "DRcaloSiPMSD.h"
#include "SCEPCALDRcaloSiPMSD.h"

#include "DD4hep/Detector.h"
#include "DDG4/Factories.h"

namespace dd4hep {
namespace sim {

  static G4VSensitiveDetector* create_DRcaloSiPM_sd(const std::string& aDetectorName, dd4hep::Detector& aLcdd) {
    std::string readoutName = aLcdd.sensitiveDetector(aDetectorName).readout().name();
    return new drc::DRcaloSiPMSD(aDetectorName,readoutName,aLcdd.sensitiveDetector(aDetectorName).readout().segmentation());
  }

  static G4VSensitiveDetector* create_SCEPCALDRcaloSiPM_sd(const std::string& aDetectorName, dd4hep::Detector& aLcdd) {
      std::string readoutName = aLcdd.sensitiveDetector(aDetectorName).readout().name();
      return new drc::SCEPCALDRcaloSiPMSD(aDetectorName,readoutName,aLcdd.sensitiveDetector(aDetectorName).readout().segmentation());
  }

}
}
DECLARE_EXTERNAL_GEANT4SENSITIVEDETECTOR(DRcaloSiPMSD,dd4hep::sim::create_DRcaloSiPM_sd)
DECLARE_EXTERNAL_GEANT4SENSITIVEDETECTOR(SCEPCALDRcaloSiPMSD,dd4hep::sim::create_SCEPCALDRcaloSiPM_sd)
