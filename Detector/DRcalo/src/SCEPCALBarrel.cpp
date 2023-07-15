//=========================================================================
// Author: Wonyong Chung
//=========================================================================

#define VERBOSE_LEVEL 0

#include "SCEPCALGridDRcaloHandle.h"

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DetectorTools.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Detector.h"

#include "TGeoTrd2.h"
#include <bitset>

using dd4hep::Transform3D;
using dd4hep::RotationZYX;
using dd4hep::Rotation3D;
using dd4hep::Position;
//using dd4hep::PlacedVolume::VolIDs;

namespace ddDRcalo {
    static dd4hep::Ref_t
    create_detector(dd4hep::Detector &theDetector, xml_h xmlElement, dd4hep::SensitiveDetector sensDet) {

      xml_det_t detectorXML=xmlElement;
      std::string name=detectorXML.nameStr();
      dd4hep::DetElement drDet(name, detectorXML.id());

      dd4hep::xml::Dimension sensDetType=detectorXML.child(_Unicode(sensitive));
      xml_comp_t dimXML=detectorXML.child(_Unicode(dim));
      xml_comp_t crystalFXML=detectorXML.child(_Unicode(crystalF));
      xml_comp_t crystalRXML=detectorXML.child(_Unicode(crystalR));
      xml_comp_t timingXML=detectorXML.child(_Unicode(timingLayer));
      xml_comp_t instXML=detectorXML.child(_Unicode(inst));
      xml_comp_t towerAssemblyXML=detectorXML.child(_Unicode(towerAssembly));

      //-----------------------------------------------------------------------------------

      sensDet.setType(sensDetType.typeStr());

      auto segmentation=dynamic_cast<dd4hep::DDSegmentation::SCEPCALGridDRcalo *>( sensDet.readout().segmentation().segmentation());

      dd4hep::Assembly experimentalHall("hall");

      const double Fdz=crystalFXML.attr<double>(_Unicode(length));
      const double Rdz=crystalRXML.attr<double>(_Unicode(length));
      const double nomfw =dimXML.attr<double>(_Unicode(towerFaceWidthNominal));
      const double nomfwT=dimXML.attr<double>(_Unicode(crystalFaceWidthNominal));
      const double EBz=dimXML.attr<double>(_Unicode(barrelHalfZ));
      const double Rin=dimXML.attr<double>(_Unicode(barrelInnerR));
      const double cube=towerAssemblyXML.attr<double>(_Unicode(cube));

      dd4hep::Material crystalFMat=theDetector.material(crystalFXML.materialStr());
      dd4hep::Material crystalRMat=theDetector.material(crystalRXML.materialStr());
      dd4hep::Material timingMat=theDetector.material(timingXML.materialStr());
      dd4hep::Material instMat=theDetector.material(instXML.materialStr());

      int nThetaBarrel=floor(EBz/nomfw);
      int nThetaEndcap=floor(Rin/nomfw);
      int nPhi=std::floor(2*M_PI*Rin/nomfw);

      double dTheta=(M_PI/2)/(nThetaBarrel+nThetaEndcap);
      double dPhi=2*M_PI/nPhi;

      for (int iTheta=0; iTheta<2*nThetaBarrel+1; iTheta++) {

        double thC=iTheta*dTheta+dTheta*nThetaEndcap;

        double r0=Rin/sin(thC);
        double r1=r0+Fdz;
        double r2=r1+Rdz;

        double y0=r0*tan(dTheta/2.);
        double y1=r1*tan(dTheta/2.);
        double y2=r2*tan(dTheta/2.);

        double x0y0 = (r0*cos(thC) +y0*sin(thC)) *tan(thC -dTheta/2.) *tan(dPhi/2.);
        double x1y0 = (r0*cos(thC) -y0*sin(thC)) *tan(thC +dTheta/2.) *tan(dPhi/2.);
        double x0y1 = (r1*cos(thC) +y1*sin(thC)) *tan(thC -dTheta/2.) *tan(dPhi/2.);
        double x1y1 = (r1*cos(thC) -y1*sin(thC)) *tan(thC +dTheta/2.) *tan(dPhi/2.);
        double x0y2 = (r2*cos(thC) +y2*sin(thC)) *tan(thC -dTheta/2.) *tan(dPhi/2.);
        double x1y2 = (r2*cos(thC) -y2*sin(thC)) *tan(thC +dTheta/2.) *tan(dPhi/2.);

        double verticesF[]={x0y0,y0,x1y0,-y0,-x1y0,-y0,-x0y0,y0,
                            x0y1,y1,x1y1,-y1,-x1y1,-y1,-x0y1,y1};

        double verticesR[]={x0y1,y1,x1y1,-y1,-x1y1,-y1,-x0y1,y1,
                            x0y2,y2,x1y2,-y2,-x1y2,-y2,-x0y2,y2};

        double verticesA[]={x0y0,y0,x1y0,-y0,-x1y0,-y0,-x0y0,y0,
                            x0y2,y2,x1y2,-y2,-x1y2,-y2,-x0y2,y2};

        dd4hep::EightPointSolid crystalFShape(Fdz/2, verticesF);
        dd4hep::EightPointSolid crystalRShape(Rdz/2, verticesR);
        dd4hep::EightPointSolid towerAssemblyShapeBarrel((Fdz+Rdz)/2, verticesA);
/*

        int nPhiT=
        int timingDir=iPhi%2-iTheta%2 //Checker pattern

        double x0y0step = 2*x0y0/nPhiT;


        for (int iPhiT=0; iPhiT<nPhiT; iPhiT++) {
          for (int iPhiT=0; iPhiT<nPhiT; iPhiT++) {

            double x0y0c = x0y0
            double x1y0c = x1y0
            double x0y1c = x0y1
            double x1y1c = x1y1
            double x0y2c = x0y2
            double x1y2c = x1y2

            double verticesF[]={x0y0,y0,x1y0,-y0,-x1y0,-y0,-x0y0,y0,
                                x0y1,y1,x1y1,-y1,-x1y1,-y1,-x0y1,y1};

            double verticesR[]={x0y1,y1,x1y1,-y1,-x1y1,-y1,-x0y1,y1,
                                x0y2,y2,x1y2,-y2,-x1y2,-y2,-x0y2,y2};

            double verticesA[]={x0y0,y0,x1y0,-y0,-x1y0,-y0,-x0y0,y0,
                                x0y2,y2,x1y2,-y2,-x1y2,-y2,-x0y2,y2};


          dd4hep::EightPointSolid crystalFShape(Fdz/2, verticesF);
          dd4hep::EightPointSolid crystalRShape(Rdz/2, verticesR);
*/

          dd4hep::Volume crystalFVol("BarrelCrystalR", crystalFShape, crystalFMat);
          dd4hep::Volume crystalRVol("BarrelCrystalF", crystalRShape, crystalRMat);

          crystalFVol.setVisAttributes(theDetector, crystalFXML.visStr());
          crystalRVol.setVisAttributes(theDetector, crystalRXML.visStr());
/*
          auto crystalFId64=segmentation->setVolumeID(1, nThetaBarrel-iTheta, iPhi, 1);
          auto crystalRId64=segmentation->setVolumeID(1, nThetaBarrel-iTheta, iPhi, 2);

          int crystalFId32=segmentation->getFirst32bits(crystalFId64);
          int crystalRId32=segmentation->getFirst32bits(crystalRId64);

          dd4hep::PlacedVolume crystalFp=towerAssemblyVol.placeVolume(crystalFVol,crystalFId32,Transform3D(RotationZYX(0,0,thC>M_PI/2?M_PI:0),Position(0,0,(thC>M_PI/2?-1:1)*(r0+Fdz/2))));
          dd4hep::PlacedVolume crystalRp=towerAssemblyVol.placeVolume(crystalRVol,crystalRId32,Transform3D(RotationZYX(0,0,thC>M_PI/2?M_PI:0),Position(0,0,(thC>M_PI/2?-1:1)*(r1+Rdz/2))));


          crystalFp.addPhysVolID("eta", nThetaBarrel-iTheta);
          crystalFp.addPhysVolID("phi", iPhi);
          crystalFp.addPhysVolID("depth", 1);

          crystalRp.addPhysVolID("eta", nThetaBarrel-iTheta);
          crystalRp.addPhysVolID("phi", iPhi);
          crystalRp.addPhysVolID("depth", 2);

          }
        }
*/

        dd4hep::Box towerAssemblyBox(cube,cube,cube);

        for (int iPhi=0; iPhi<nPhi; iPhi++) {

          auto crystalFId64=segmentation->setVolumeID(1, (nThetaBarrel-iTheta) *(iTheta>nThetaBarrel?-1:1), iPhi, 1);
          auto crystalRId64=segmentation->setVolumeID(1, (nThetaBarrel-iTheta) *(iTheta>nThetaBarrel?-1:1), iPhi, 2);

          int crystalFId32=segmentation->getFirst32bits(crystalFId64);
          int crystalRId32=segmentation->getFirst32bits(crystalRId64);


          dd4hep::Volume towerAssemblyVol("towerAssemblyVol", towerAssemblyShapeBarrel, theDetector.material("Vacuum"));
          towerAssemblyVol.setVisAttributes(theDetector, towerAssemblyXML.visStr());

          double rt=r0+(Fdz+Rdz)/2.;
          double phi=iPhi*dPhi;

          RotationZYX rot(M_PI/2, thC, 0);
          ROOT::Math::RotationZ rotZ = ROOT::Math::RotationZ(phi);
          rot = rotZ*rot;

          Position disp(rt*sin(thC)*cos(phi),
                        rt*sin(thC)*sin(phi),
                        rt*cos(thC));

          experimentalHall.placeVolume( towerAssemblyVol, Transform3D(rot,disp) );

          dd4hep::PlacedVolume crystalFp = towerAssemblyVol.placeVolume( crystalFVol, crystalFId32, Position(0,0,-Rdz/2) );
          dd4hep::PlacedVolume crystalRp = towerAssemblyVol.placeVolume( crystalRVol, crystalRId32, Position(0,0,Fdz/2) );

          crystalFp.addPhysVolID("eta", nThetaBarrel-iTheta *(iTheta>nThetaBarrel?-1:1));
          crystalFp.addPhysVolID("phi", iPhi);
          crystalFp.addPhysVolID("depth", 1);
          crystalFp.addPhysVolID("system", 1);

          crystalRp.addPhysVolID("eta", nThetaBarrel-iTheta *(iTheta>nThetaBarrel?-1:1));
          crystalRp.addPhysVolID("phi", iPhi);
          crystalRp.addPhysVolID("depth", 2);
          crystalRp.addPhysVolID("system", 1);

          std::bitset<10> _eta((nThetaBarrel-iTheta) *(iTheta>nThetaBarrel?-1:1));
          std::bitset<10> _phi(iPhi);
          std::bitset<3> depthF(1);
          std::bitset<3> depthR(2);
          std::bitset<32> id32F(crystalFId32);
          std::bitset<32> id32R(crystalRId32);

//          VolIDs crystalFpVID = static_cast<VolIDs> crystalFp.VolIDs();
//          VolIDs crystalRpVID = static_cast<VolIDs> crystalRp.VolIDs();

          std::cout << "B crystalF eta: " << ((nThetaBarrel-iTheta) *(iTheta>nThetaBarrel?-1:1)) << " phi: " << iPhi << " depth: " << 1 << std::endl;
          std::cout << "B crystalF eta: " << _eta << " phi: " << _phi << " depth: " << depthF << std::endl;
          std::cout << "B crystalFId32: " << id32F << std::endl;
//          std::cout << "B crystalF copyNum: " << crystalFp.copyNumber() << " VolIDs: " << dd4hep::detail::tools::toString( crystalFpVID ) << std::endl;

          std::cout << "B crystalR eta: " << ((nThetaBarrel-iTheta) *(iTheta>nThetaBarrel?-1:1)) << " phi: " << iPhi << " depth: " << 2 << std::endl;
          std::cout << "B crystalR eta: " << _eta << " phi: " << _phi << " depth: " << depthR << std::endl;
          std::cout << "B crystalRId32: " << id32R << std::endl;
//          std::cout << "B crystalR copyNum: " << crystalRp.copyNumber() << " VolIDs: " << dd4hep::detail::tools::toString( crystalRpVID ) << std::endl;
        }
      }

      for (int iTheta=0; iTheta<nThetaEndcap; iTheta++) {
        if (abs(iTheta)<2) continue;

        double thC=iTheta*dTheta;

        double r0=EBz/cos(thC);
        double y0=r0*tan(dTheta/2.);

        double r1=r0+Fdz;
        double y1=r1*tan(dTheta/2.);

        double r2=r1+Rdz;
        double y2=r2*tan(dTheta/2.);

        double verticesF[]={(r0*cos(thC)+y0*sin(thC))*tan(thC-dTheta/2.)*tan(dPhi/2.), y0,
                            (r0*cos(thC)-y0*sin(thC))*tan(thC+dTheta/2.)*tan(dPhi/2.),-y0,
                           -(r0*cos(thC)-y0*sin(thC))*tan(thC+dTheta/2.)*tan(dPhi/2.),-y0,
                           -(r0*cos(thC)+y0*sin(thC))*tan(thC-dTheta/2.)*tan(dPhi/2.), y0,

                            (r1*cos(thC)+y1*sin(thC))*tan(thC-dTheta/2.)*tan(dPhi/2.), y1,
                            (r1*cos(thC)-y1*sin(thC))*tan(thC+dTheta/2.)*tan(dPhi/2.),-y1,
                           -(r1*cos(thC)-y1*sin(thC))*tan(thC+dTheta/2.)*tan(dPhi/2.),-y1,
                           -(r1*cos(thC)+y1*sin(thC))*tan(thC-dTheta/2.)*tan(dPhi/2.), y1};

        double verticesR[]={(r1*cos(thC)+y1*sin(thC))*tan(thC-dTheta/2.)*tan(dPhi/2.), y1,
                            (r1*cos(thC)-y1*sin(thC))*tan(thC+dTheta/2.)*tan(dPhi/2.),-y1,
                           -(r1*cos(thC)-y1*sin(thC))*tan(thC+dTheta/2.)*tan(dPhi/2.),-y1,
                           -(r1*cos(thC)+y1*sin(thC))*tan(thC-dTheta/2.)*tan(dPhi/2.), y1,

                            (r2*cos(thC)+y2*sin(thC))*tan(thC-dTheta/2.)*tan(dPhi/2.), y2,
                            (r2*cos(thC)-y2*sin(thC))*tan(thC+dTheta/2.)*tan(dPhi/2.),-y2,
                           -(r2*cos(thC)-y2*sin(thC))*tan(thC+dTheta/2.)*tan(dPhi/2.),-y2,
                           -(r2*cos(thC)+y2*sin(thC))*tan(thC-dTheta/2.)*tan(dPhi/2.), y2};

        double verticesA[]={(r0*cos(thC)+y0*sin(thC))*tan(thC-dTheta/2.)*tan(dPhi/2.), y0,
                            (r0*cos(thC)-y0*sin(thC))*tan(thC+dTheta/2.)*tan(dPhi/2.),-y0,
                           -(r0*cos(thC)-y0*sin(thC))*tan(thC+dTheta/2.)*tan(dPhi/2.),-y0,
                           -(r0*cos(thC)+y0*sin(thC))*tan(thC-dTheta/2.)*tan(dPhi/2.), y0,

                            (r2*cos(thC)+y2*sin(thC))*tan(thC-dTheta/2.)*tan(dPhi/2.), y2,
                            (r2*cos(thC)-y2*sin(thC))*tan(thC+dTheta/2.)*tan(dPhi/2.),-y2,
                           -(r2*cos(thC)-y2*sin(thC))*tan(thC+dTheta/2.)*tan(dPhi/2.),-y2,
                           -(r2*cos(thC)+y2*sin(thC))*tan(thC-dTheta/2.)*tan(dPhi/2.), y2};


        dd4hep::EightPointSolid crystalFShape(Fdz/2, verticesF);
        dd4hep::EightPointSolid crystaslRShape(Rdz/2, verticesR);
        dd4hep::EightPointSolid towerAssemblyShapeEndcap((Fdz+Rdz)/2, verticesA);

        dd4hep::Volume crystalFVol("EndcapCrystalR", crystalFShape, crystalFMat);
        dd4hep::Volume crystalRVol("EndcapCrystalF", crystaslRShape, crystalRMat);

        crystalFVol.setVisAttributes(theDetector, crystalFXML.visStr());
        crystalRVol.setVisAttributes(theDetector, crystalRXML.visStr());

        dd4hep::Box towerAssemblyBoxEndcap(cube,cube,cube);

        for (int iPhi=0; iPhi<nPhi; iPhi++) {

          auto crystalFId64=segmentation->setVolumeID(1,(nThetaBarrel+nThetaEndcap-iTheta), iPhi, 1);
          auto crystalRId64=segmentation->setVolumeID(1,(nThetaBarrel+nThetaEndcap-iTheta), iPhi, 2);
          auto crystalFId641=segmentation->setVolumeID(1, -(nThetaBarrel+nThetaEndcap-iTheta), iPhi, 1);
          auto crystalRId641=segmentation->setVolumeID(1, -(nThetaBarrel+nThetaEndcap-iTheta), iPhi, 2);

          int crystalFId32=segmentation->getFirst32bits(crystalFId64);
          int crystalRId32=segmentation->getFirst32bits(crystalRId64);
          int crystalFId321=segmentation->getFirst32bits(crystalFId641);
          int crystalRId321=segmentation->getFirst32bits(crystalRId641);

          dd4hep::Volume towerAssemblyVolEndcap("towerAssemblyVolEndcap", towerAssemblyShapeEndcap, theDetector.material("Vacuum"));
          dd4hep::Volume towerAssemblyVolEndcap1("towerAssemblyVolEndcap1", towerAssemblyShapeEndcap, theDetector.material("Vacuum"));
          towerAssemblyVolEndcap.setVisAttributes(theDetector, towerAssemblyXML.visStr());
          towerAssemblyVolEndcap1.setVisAttributes(theDetector, towerAssemblyXML.visStr());

          double rt=r0+(Fdz+Rdz)/2.;
          double phi=iPhi*dPhi;

          ROOT::Math::RotationZ rotZ = ROOT::Math::RotationZ(phi);

          RotationZYX rot(M_PI/2, thC, 0);
          rot = rotZ*rot;
          Position disp(rt*sin(thC)*cos(phi),
                        rt*sin(thC)*sin(phi),
                        rt*cos(thC));

          RotationZYX rot1(M_PI/2, thC, M_PI);
          rot1 = rotZ*rot1;
          Position disp1(rt*sin(thC)*cos(phi),
                        rt*sin(thC)*sin(phi),
                        -rt*cos(thC));

          experimentalHall.placeVolume( towerAssemblyVolEndcap, Transform3D(rot,disp) );
          experimentalHall.placeVolume( towerAssemblyVolEndcap1, Transform3D(rot1,disp1) );

          dd4hep::PlacedVolume crystalFp = towerAssemblyVolEndcap.placeVolume( crystalFVol, crystalFId32, Position(0,0,-Rdz/2) );
          dd4hep::PlacedVolume crystalRp = towerAssemblyVolEndcap.placeVolume( crystalRVol, crystalRId32, Position(0,0,Fdz/2) );

          dd4hep::PlacedVolume crystalFp1 = towerAssemblyVolEndcap1.placeVolume( crystalFVol, crystalFId321, Position(0,0,-Rdz/2) );
          dd4hep::PlacedVolume crystalRp1 = towerAssemblyVolEndcap1.placeVolume( crystalRVol, crystalRId321, Position(0,0,Fdz/2) );

          crystalFp.addPhysVolID("eta", (nThetaBarrel+nThetaEndcap-iTheta));
          crystalFp.addPhysVolID("phi", iPhi);
          crystalFp.addPhysVolID("depth", 1);
          crystalFp.addPhysVolID("system", 1);

          crystalRp.addPhysVolID("eta", (nThetaBarrel+nThetaEndcap-iTheta));
          crystalRp.addPhysVolID("phi", iPhi);
          crystalRp.addPhysVolID("depth", 2);
          crystalRp.addPhysVolID("system", 1);


          crystalFp1.addPhysVolID("eta", -(nThetaBarrel+nThetaEndcap-iTheta));
          crystalFp1.addPhysVolID("phi", iPhi);
          crystalFp1.addPhysVolID("depth", 1);
          crystalFp1.addPhysVolID("system", 1);

          crystalRp1.addPhysVolID("eta", -(nThetaBarrel+nThetaEndcap-iTheta));
          crystalRp1.addPhysVolID("phi", iPhi);
          crystalRp1.addPhysVolID("depth", 2);
          crystalRp1.addPhysVolID("system", 1);

          std::bitset<10> _eta((nThetaBarrel+nThetaEndcap-iTheta));
          std::bitset<10> _eta1(-(nThetaBarrel+nThetaEndcap-iTheta));
          std::bitset<10> _phi(iPhi);
          std::bitset<3> depthF(1);
          std::bitset<3> depthR(2);
          std::bitset<32> id32F(crystalFId32);
          std::bitset<32> id32F1(crystalFId321);
          std::bitset<32> id32R(crystalRId32);
          std::bitset<32> id32R1(crystalRId321);



          std::cout << "E crystalF eta: " << (nThetaBarrel+nThetaEndcap-iTheta) << " phi: " << iPhi << " depth: " << 1 << std::endl;
          std::cout << "E crystalF eta: " << _eta << " phi: " << _phi << " depth: " << depthF << std::endl;
          std::cout << "E crystalFId32: " << id32F << std::endl;
//          std::cout << "E crystalF copyNum: " << crystalFp.copyNumber() << " VolIDs: " << dd4hep::detail::tools::toString( crystalFp.VolIDs()) << std::endl;

          std::cout << "E crystalR eta: " << (nThetaBarrel+nThetaEndcap-iTheta) << " phi: " << iPhi << " depth: " << 2 << std::endl;
          std::cout << "E crystalR eta: " << _eta << " phi: " << _phi << " depth: " << depthR << std::endl;
          std::cout << "E crystalRId32: " << id32R << std::endl;
//          std::cout << "E crystalR copyNum: " << crystalRp.copyNumber() << " VolIDs: " << dd4hep::detail::tools::toString( crystalRp.VolIDs()) << std::endl;

          std::cout << "E crystalF1 eta: " << -(nThetaBarrel+nThetaEndcap-iTheta) << " phi: " << iPhi << " depth: " << 1 << std::endl;
          std::cout << "E crystalF1 eta: " << _eta1 << " phi: " << _phi << " depth: " << depthF << std::endl;
          std::cout << "E crystalFId321: " << id32F1 << std::endl;
//          std::cout << "E crystalF1 copyNum: " << crystalFp1.copyNumber() << " VolIDs: " << dd4hep::detail::tools::toString(crystalFp1.VolIDs()) << std::endl;

          std::cout << "E crystalR1 eta: " << -(nThetaBarrel+nThetaEndcap-iTheta) << " phi: " << iPhi << " depth: " << 2 << std::endl;
          std::cout << "E crystalR1 eta: " << _eta1 << " phi: " << _phi << " depth: " << depthR << std::endl;
          std::cout << "E crystalRId321: " << id32R1 << std::endl;
//          std::cout << "E crystalR1 copyNum: " << crystalRp1.copyNumber() << " VolIDs: " << dd4hep::detail::tools::toString(crystalRp1.VolIDs()) << std::endl;
        }


      }

      dd4hep::PlacedVolume hallPlace=theDetector.pickMotherVolume(drDet).placeVolume(experimentalHall);
//      hallPlace.addPhysVolID("system", detectorXML.id());
      hallPlace.addPhysVolID("system", 0);

      drDet.setPlacement(hallPlace);

      return drDet;
    }
}


DECLARE_DETELEMENT(SCEPCALBarrel, ddDRcalo::create_detector)
