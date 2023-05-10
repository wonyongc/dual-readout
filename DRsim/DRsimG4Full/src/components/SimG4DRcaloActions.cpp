#include "SimG4DRcaloActions.h"

DECLARE_COMPONENT(SimG4DRcaloActions)

SimG4DRcaloActions::SimG4DRcaloActions(const std::string& type, const std::string& name, const IInterface* parent)
: AlgTool(type, name, parent), m_geoSvc("GeoSvc", name) {
  declareInterface<ISimG4ActionTool>(this);
}

SimG4DRcaloActions::~SimG4DRcaloActions() {}

StatusCode SimG4DRcaloActions::initialize() {
  if (AlgTool::initialize().isFailure())
    return StatusCode::FAILURE;

  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }

  pSegHCAL = dynamic_cast<dd4hep::DDSegmentation::GridDRcalo*>(m_geoSvc->lcdd()->readout(m_readoutNameHCAL).segmentation().segmentation());
  pSegECAL = dynamic_cast<dd4hep::DDSegmentation::SCEPCALGridDRcalo*>(m_geoSvc->lcdd()->readout(m_readoutNameECAL).segmentation().segmentation());

  return StatusCode::SUCCESS;
}

StatusCode SimG4DRcaloActions::finalize() { return AlgTool::finalize(); }

G4VUserActionInitialization* SimG4DRcaloActions::userActionInitialization() {
  auto* actions = new drc::SimG4DRcaloActionInitialization();

  actions->setSegmentationHCAL(pSegHCAL);
  actions->setBirksConstantHCAL(m_scintNameHCAL,m_birksHCAL);
  actions->setThresholdHCAL(m_thresHCAL);

  actions->setSegmentationECAL(pSegECAL);
  actions->setBirksConstantECAL(m_scintNameECAL,m_birksECAL);
  actions->setThresholdECAL(m_thresECAL);

  return actions;
}
