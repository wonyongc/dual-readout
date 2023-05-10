//=========================================================================
// Author: Wonyong Chung
//=========================================================================

#define VERBOSE_LEVEL 0

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "XML/Layering.h"
#include "TGeoTrd2.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"
#include "DDSegmentation/Segmentation.h"

using namespace std;

using dd4hep::BUILD_ENVELOPE;
using dd4hep::Box;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::Layering;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::Readout;
using dd4hep::Ref_t;
using dd4hep::RotationZYX;
using dd4hep::Segmentation;
using dd4hep::SensitiveDetector;
using dd4hep::Transform3D;
using dd4hep::Translation3D;
using dd4hep::Trapezoid;
using dd4hep::EightPointSolid;
using dd4hep::Volume;
using dd4hep::_toString;

using dd4hep::rec::LayeredCalorimeterData;

static Ref_t create_detector(Detector& theDetector, xml_h e, SensitiveDetector sens)  {

    xml_det_t     x_det     = e;
    int           det_id    = x_det.id();
    string        det_name  = x_det.nameStr();
    DetElement    sdet      (det_name, det_id);

    // --- create an envelope volume and position it into the world ---------------------

    Volume envelope = dd4hep::xml::createPlacedEnvelope(theDetector, e, sdet);
    dd4hep::xml::setDetectorTypeFlag( e, sdet ) ;

    if(theDetector.buildType() == BUILD_ENVELOPE) return sdet;

    //-----------------------------------------------------------------------------------
    dd4hep::xml::Component xmlParameter = x_det.child(_Unicode(parameter));

    const double Rin = xmlParameter.attr<double>(_Unicode(Rmin));
    const double EBz = xmlParameter.attr<double>(_Unicode(EBz));
    const double nomfw = xmlParameter.attr<double>(_Unicode(nominalfacewidth));
    const double Fdz = xmlParameter.attr<double>(_Unicode(Fdz));
    const double Rdz = xmlParameter.attr<double>(_Unicode(Rdz));

    const string region = xmlParameter.attr<string>(_Unicode(region));
    const string cal_limits = xmlParameter.attr<string>(_Unicode(limits));
    const string cFvis = xmlParameter.attr<string>(_Unicode(cFvis));
    const string cRvis = xmlParameter.attr<string>(_Unicode(cRvis));
    const string cCvis = xmlParameter.attr<string>(_Unicode(cCvis));

    Material PbWO4 = theDetector.material(xmlParameter.attr<string>(_Unicode(crystalMaterial)));
//    Material LYSO = theDetector.material(xmlParameter.attr<string>(_Unicode(crystalMaterial)));
//    Material Aluminum = theDetector.material(xmlParameter.attr<string>(_Unicode(instMaterial)));

    int nThetaBarrel = floor(EBz/nomfw);
    int nThetaEndcap = floor(Rin/nomfw);
    double dTheta = (M_PI/2)/(nThetaBarrel+nThetaEndcap);

    double thetaBarrel = dTheta*nThetaBarrel;
    double thetaEndcap = dTheta*nThetaEndcap;

    int nTheta = nThetaEndcap;
    double thetaRange = thetaEndcap;

    int nPhi = floor(2*M_PI*Rin/nomfw);
    double dPhi = 2*M_PI/nPhi;

    Readout readout = sens.readout();
    Segmentation seg = readout.segmentation();

    std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0);
    //Assume uniform cell sizes, provide dummy cellID
    double cell_sizeX      = cellSizeVector[0];
    double cell_sizeY      = cellSizeVector[1];

    Layering      layering (e);
    Material      air       = theDetector.air();
    xml_comp_t    x_staves  = x_det.staves();
    // Create extension objects for reconstruction -----------------

    // Create caloData object to extend driver with data required for reconstruction

    LayeredCalorimeterData* caloData = new LayeredCalorimeterData;
    caloData->layoutType = LayeredCalorimeterData::BarrelLayout;
    caloData->inner_symmetry = nPhi;
    caloData->outer_symmetry = nPhi;

    caloData->inner_phi0 = 0.;
    caloData->outer_phi0 = 0.;

    // One should ensure that these sensitivity gaps are correctly used
    caloData->gap0 = 0;  // the 4 gaps between the 5 towers, along z
    caloData->gap1 = 0; // gaps between stacks in a module, along z
    caloData->gap2 = 0; // gaps where the staves overlap

    // extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
    caloData->extent[0] = Rin;
    caloData->extent[1] = Rin+Fdz+Rdz; // or r_max ?
    caloData->extent[2] = 0.;      // NN: for barrel detectors this is 0
    caloData->extent[3] = 2*EBz;

    // Compute stave dimensions
    double rStave = EBz/cos(nTheta*dTheta);

    double verticesStave[] = {rStave*tan(thetaRange+dTheta/2), rStave*tan(dPhi/2),
                              rStave*tan(thetaRange+dTheta/2), -rStave*tan(dPhi/2),
                              -rStave*tan(thetaRange+dTheta/2), -rStave*tan(dPhi/2),
                              -rStave*tan(thetaRange+dTheta/2), rStave*tan(dPhi/2),
                              (rStave+Fdz+Rdz)*tan(thetaRange+dTheta/2), (rStave+Fdz+Rdz)*tan(dPhi/2),
                              (rStave+Fdz+Rdz)*tan(thetaRange+dTheta/2), -(rStave+Fdz+Rdz)*tan(dPhi/2),
                              -(rStave+Fdz+Rdz)*tan(thetaRange+dTheta/2), -(rStave+Fdz+Rdz)*tan(dPhi/2),
                              -(rStave+Fdz+Rdz)*tan(thetaRange+dTheta/2), (rStave+Fdz+Rdz)*tan(dPhi/2)};

    EightPointSolid stave((Fdz+Rdz)/2, verticesStave);

    Volume stave_vol("stave", stave, air);

    sens.setType("calorimeter");

    DetElement stave_det("stave1", det_id);

    LayeredCalorimeterData::Layer caloLayerCrystalsF;
    LayeredCalorimeterData::Layer caloLayerCrystalsR;

    caloLayerCrystalsF.distance = Rin;
    caloLayerCrystalsF.sensitive_thickness       = Fdz ;
    caloLayerCrystalsF.inner_nRadiationLengths   = 1;
    caloLayerCrystalsF.inner_nInteractionLengths = 1;
    caloLayerCrystalsF.inner_thickness           = Fdz;
    caloLayerCrystalsF.outer_nRadiationLengths   = 1;
    caloLayerCrystalsF.outer_nInteractionLengths = 1;
    caloLayerCrystalsF.outer_thickness           = Fdz;
    caloLayerCrystalsF.absorberThickness         = Fdz;
    caloLayerCrystalsF.cellSize0 = cell_sizeX;
    caloLayerCrystalsF.cellSize1 = cell_sizeY;

    caloLayerCrystalsR.distance = Rin;
    caloLayerCrystalsR.sensitive_thickness       = Rdz ;
    caloLayerCrystalsR.inner_nRadiationLengths   = 1;
    caloLayerCrystalsR.inner_nInteractionLengths = 1;
    caloLayerCrystalsR.inner_thickness           = Rdz;
    caloLayerCrystalsR.outer_nRadiationLengths   = 1;
    caloLayerCrystalsR.outer_nInteractionLengths = 1;
    caloLayerCrystalsR.outer_thickness           = Rdz;
    caloLayerCrystalsR.absorberThickness         = Rdz;
    caloLayerCrystalsR.cellSize0 = cell_sizeX;
    caloLayerCrystalsR.cellSize1 = cell_sizeY;

    caloData->layers.push_back( caloLayerCrystalsF ) ;
    caloData->layers.push_back( caloLayerCrystalsR ) ;

    for (int t = 0 ; t < nTheta; t++) {
        if (abs(t) < 2) continue;

        double thC = t*dTheta;

        double r0 = EBz/cos(thC);
        double y0 = r0*tan(dTheta/2);

        double r1 = r0+Fdz;
        double y1 = r1*tan(dTheta/2);

        double r2 = r1+Rdz;
        double y2 = r2*tan(dTheta/2);

        double verticesF[] = { (r0*cos(thC)+y0*sin(thC)) *tan(thC-dTheta/2) *tan(dPhi/2),  y0,
                               (r0*cos(thC)-y0*sin(thC)) *tan(thC+dTheta/2) *tan(dPhi/2), -y0,
                              -(r0*cos(thC)-y0*sin(thC)) *tan(thC+dTheta/2) *tan(dPhi/2), -y0,
                              -(r0*cos(thC)+y0*sin(thC)) *tan(thC-dTheta/2) *tan(dPhi/2),  y0,

                               (r1*cos(thC)+y1*sin(thC)) *tan(thC-dTheta/2) *tan(dPhi/2),  y1,
                               (r1*cos(thC)-y1*sin(thC)) *tan(thC+dTheta/2) *tan(dPhi/2), -y1,
                               -(r1*cos(thC)-y1*sin(thC)) *tan(thC+dTheta/2) *tan(dPhi/2), -y1,
                               -(r1*cos(thC)+y1*sin(thC)) *tan(thC-dTheta/2) *tan(dPhi/2),  y1};

        double verticesR[] = { (r1*cos(thC)+y1*sin(thC)) *tan(thC-dTheta/2) *tan(dPhi/2),  y1,
                               (r1*cos(thC)-y1*sin(thC)) *tan(thC+dTheta/2) *tan(dPhi/2), -y1,
                               -(r1*cos(thC)-y1*sin(thC)) *tan(thC+dTheta/2) *tan(dPhi/2), -y1,
                               -(r1*cos(thC)+y1*sin(thC)) *tan(thC-dTheta/2) *tan(dPhi/2),  y1,

                               (r2*cos(thC)+y2*sin(thC)) *tan(thC-dTheta/2) *tan(dPhi/2),  y2,
                               (r2*cos(thC)-y2*sin(thC)) *tan(thC+dTheta/2) *tan(dPhi/2), -y2,
                               -(r2*cos(thC)-y2*sin(thC)) *tan(thC+dTheta/2) *tan(dPhi/2), -y2,
                               -(r2*cos(thC)+y2*sin(thC)) *tan(thC-dTheta/2) *tan(dPhi/2),  y2};

        double verticesC[] = { (r0*cos(thC)+y0*sin(thC)) *tan(thC-dTheta/2) *tan(dPhi/2),  y0,
                               (r0*cos(thC)-y0*sin(thC)) *tan(thC+dTheta/2) *tan(dPhi/2), -y0,
                               -(r0*cos(thC)-y0*sin(thC)) *tan(thC+dTheta/2) *tan(dPhi/2), -y0,
                               -(r0*cos(thC)+y0*sin(thC)) *tan(thC-dTheta/2) *tan(dPhi/2),  y0,

                               (r2*cos(thC)+y2*sin(thC)) *tan(thC-dTheta/2) *tan(dPhi/2),  y2,
                               (r2*cos(thC)-y2*sin(thC)) *tan(thC+dTheta/2) *tan(dPhi/2), -y2,
                               -(r2*cos(thC)-y2*sin(thC)) *tan(thC+dTheta/2) *tan(dPhi/2), -y2,
                               -(r2*cos(thC)+y2*sin(thC)) *tan(thC-dTheta/2) *tan(dPhi/2),  y2};


        EightPointSolid crystalsegmentF(Fdz/2, verticesF);
        EightPointSolid crystalsegmentR(Rdz/2, verticesR);
        EightPointSolid crystalsegmentC((Fdz+Rdz)/2, verticesC);

        string cF_name = _toString(t, "crystalF%d");
        string cR_name = _toString(t, "crystalR%d");
        string cC_name = _toString(t, "crystalC%d");

        Volume crystalF(cF_name, crystalsegmentF, PbWO4);
        Volume crystalR(cR_name, crystalsegmentR, PbWO4);
        Volume crystalC(cC_name, crystalsegmentC, air);

        DetElement cF_det(stave_det, cF_name, det_id);
        DetElement cR_det(stave_det, cR_name, det_id);
        DetElement cC_det(stave_det, cC_name, det_id);

        crystalF.setSensitiveDetector(sens);
        crystalR.setSensitiveDetector(sens);

        crystalF.setAttributes(theDetector, region, cal_limits, cFvis);
        crystalR.setAttributes(theDetector, region, cal_limits, cRvis);
        crystalC.setAttributes(theDetector, region, cal_limits, cCvis);

        PlacedVolume cF = crystalC.placeVolume(crystalF, Position(0,0,r0+Fdz/2));
        PlacedVolume cR = crystalC.placeVolume(crystalR, Position(0,0,r1+Rdz/2));

        cF.addPhysVolID("layer", t);
        cF.addPhysVolID("layer", t);

        cF_det.setPlacement(cF);
        cR_det.setPlacement(cR);

        // Place the staves
        for (int i = 0; i < nPhi; i++) { // i is the stave number
            stave_vol.placeVolume(crystalC, RotationZYX(i*dPhi, thC, 0));
        }

    }

    stave_det.setVisAttributes(theDetector, x_staves.visStr(), stave_vol);

    DetElement stave_det2 = stave_det.clone("stave2");


    // Place the staves
    PlacedVolume pv1 = envelope.placeVolume(stave_vol, RotationZYX(0, 0, 0));
    pv1.addPhysVolID("side", 1);

    PlacedVolume pv2 = envelope.placeVolume(stave_vol, Transform3D(RotationZYX(nPhi%2==1? dPhi/2:0, 0, M_PI), Position(0,0,0)) );
    pv2.addPhysVolID("side", 2);


    stave_det.setPlacement(pv1);
    stave_det2.setPlacement(pv2);

    sdet.add(stave_det);
    sdet.add(stave_det2);

    // Set envelope volume attributes
    envelope.setAttributes(theDetector,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());

    sdet.addExtension< LayeredCalorimeterData >( caloData ) ;

    return sdet;
}

DECLARE_DETELEMENT(SCEPCALEndcap, create_detector)
