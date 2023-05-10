#include "SCEPCALGridDRcalo.h"

#include <climits>
#include <cmath>
#include <stdexcept>

namespace dd4hep {
namespace DDSegmentation {

/// default constructor using an encoding string
        SCEPCALGridDRcalo::SCEPCALGridDRcalo(const std::string& cellEncoding) : Segmentation(cellEncoding) {

            _type = "SCEPCALGridDRcalo";
            _description = "DRcalo segmentation based on the tower / (Cerenkov or Scintillation) fiber / SiPM hierarchy";

            registerIdentifier("identifier_system", "Cell ID identifier for numSystem", fSystemId, "system");
            registerIdentifier("identifier_eta", "Cell ID identifier for numEta", fEtaId, "eta");
            registerIdentifier("identifier_phi", "Cell ID identifier for numPhi", fPhiId, "phi");
            registerIdentifier("identifier_depth", "Cell ID identifier for numDepth", fDepthId, "depth");
            registerIdentifier("identifier_sipm", "Cell ID identifier for numSipm", fSipmId, "sipm");
            registerIdentifier("identifier_IsCerenkov", "Cell ID identifier for IsCerenkov", fIsCerenkovId, "c");
            registerIdentifier("identifier_module", "Cell ID identifier for module", fModule, "module");

        }

        SCEPCALGridDRcalo::SCEPCALGridDRcalo(const BitFieldCoder* decoder) : Segmentation(decoder) {

            _type = "SCEPCALGridDRcalo";
            _description = "DRcalo segmentation based on the tower / (Cerenkov or Scintillation) fiber / SiPM hierarchy";

            registerIdentifier("identifier_system", "Cell ID identifier for numSystem", fSystemId, "system");
            registerIdentifier("identifier_eta", "Cell ID identifier for Eta", fEtaId, "eta");
            registerIdentifier("identifier_phi", "Cell ID identifier for Phi", fPhiId, "phi");
            registerIdentifier("identifier_depth", "Cell ID identifier for Depth", fDepthId, "depth");
            registerIdentifier("identifier_sipm", "Cell ID identifier for Sipm", fSipmId, "sipm");
            registerIdentifier("identifier_IsCerenkov", "Cell ID identifier for IsCerenkov", fIsCerenkovId, "c");
            registerIdentifier("identifier_module", "Cell ID identifier for module", fModule, "module");

        }

        SCEPCALGridDRcalo::~SCEPCALGridDRcalo() {}

        Vector3D SCEPCALGridDRcalo::position(const CellID& cID) const {

          int nEta_in = Eta(cID);
          int nPhi_in = Phi(cID);
          int nDepth_in = Depth(cID);

          double EBz = 2.0;
          double Rin = 1.5;
          double nomfw = 0.1;
          double Fdz = 0.05;
          double Rdz = 0.15;

          int nThetaBarrel=floor(EBz/nomfw);
          int nThetaEndcap=floor(Rin/nomfw);
          int nPhi=std::floor(2*M_PI*Rin/nomfw);

          double dTheta=(M_PI/2)/(nThetaBarrel+nThetaEndcap);
          double dPhi=2*M_PI/nPhi;


          int nTheta = nEta_in>0? (nThetaBarrel+nThetaEndcap)-nEta_in : (nThetaBarrel+nThetaEndcap)-nEta_in;
          double thC = nTheta*dTheta;
          double phi = nPhi_in*dPhi;

          double r0=EBz/cos(thC);
          double r1=r0+Fdz;
          double r2=r1+Rdz;

          double R = nDepth_in==1? (r0+r1)/2:(r1+r2)/2;
          double x = R*sin(thC)*cos(phi);
          double y = R*sin(thC)*sin(phi);
          double z = R*cos(thC);

          return Vector3D(x,y,z);

        }

/*
        Vector3D SCEPCALGridDRcalo::localPosition(const CellID& cID) const {
            int numx = numX(cID);
            int numy = numY(cID);
            int numz = numZ(cID);
            int x_ = x(cID);
            int y_ = y(cID);
            int z_ = z(cID);

            return localPosition(numx,numy,numz,x_,y_,z_);
        }

        Vector3D SCEPCALGridDRcalo::localPosition(int numx, int numy, int numz, int x_, int y_, int z_) const {
            float ptX = -fGridSize*static_cast<float>(numx/2) + static_cast<float>(x_)*fGridSize + ( numx%2==0 ? fGridSize/2. : 0. );
            float ptY = -fGridSize*static_cast<float>(numy/2) + static_cast<float>(y_)*fGridSize + ( numy%2==0 ? fGridSize/2. : 0. );
            float ptZ = -fGridSize*static_cast<float>(numz/2) + static_cast<float>(z_)*fGridSize + ( numz%2==0 ? fGridSize/2. : 0. );

            return Vector3D(ptX,ptY,ptZ);
        }
*/

        CellID SCEPCALGridDRcalo::cellID(const Vector3D& /*localPosition*/, const Vector3D& /*globalPosition*/, const VolumeID& vID) const {

            return setCellID(System(vID), Eta(vID), Phi(vID), Depth(vID) );
        }

        VolumeID SCEPCALGridDRcalo::setVolumeID(int System, int Eta, int Phi, int Depth) const {
            VolumeID SystemId = static_cast<VolumeID>(System);
            VolumeID EtaId = static_cast<VolumeID>(Eta);
            VolumeID PhiId = static_cast<VolumeID>(Phi);
            VolumeID DepthId = static_cast<VolumeID>(Depth);
            VolumeID vID = 0;
            _decoder->set(vID, fSystemId, SystemId);
            _decoder->set(vID, fEtaId, EtaId);
            _decoder->set(vID, fPhiId, PhiId);
            _decoder->set(vID, fDepthId, DepthId);

            VolumeID module = 0; // Tower
            _decoder->set(vID, fModule, module);

            return vID;
        }

        CellID SCEPCALGridDRcalo::setCellID(int System, int Eta, int Phi, int Depth) const {
          VolumeID SystemId = static_cast<VolumeID>(System);
          VolumeID EtaId = static_cast<VolumeID>(Eta);
            VolumeID PhiId = static_cast<VolumeID>(Phi);
            VolumeID DepthId = static_cast<VolumeID>(Depth);
            VolumeID vID = 0;
            _decoder->set(vID, fSystemId, SystemId);
            _decoder->set(vID, fEtaId, EtaId);
            _decoder->set(vID, fPhiId, PhiId);
            _decoder->set(vID, fDepthId, DepthId);

            VolumeID module = 1; // Fiber, SiPM, etc.
            _decoder->set(vID, fModule, module);

            VolumeID isCeren = IsCerenkov(Eta,Phi) ? 1 : 0;
            _decoder->set(vID, fIsCerenkovId, isCeren);

            return vID;
        }

        int SCEPCALGridDRcalo::System(const CellID& aCellID) const {
          VolumeID System = static_cast<VolumeID>(_decoder->get(aCellID, fSystemId));
          return static_cast<int>(System);
        }

        int SCEPCALGridDRcalo::Eta(const CellID& aCellID) const {
            VolumeID Eta = static_cast<VolumeID>(_decoder->get(aCellID, fEtaId));
            return static_cast<int>(Eta);
        }

        int SCEPCALGridDRcalo::Phi(const CellID& aCellID) const {
            VolumeID Phi = static_cast<VolumeID>(_decoder->get(aCellID, fPhiId));
            return static_cast<int>(Phi);
        }

        int SCEPCALGridDRcalo::Depth(const CellID& aCellID) const {
            VolumeID Depth = static_cast<VolumeID>(_decoder->get(aCellID, fDepthId));
            return static_cast<int>(Depth);
        }

        int SCEPCALGridDRcalo::Sipm(const CellID& aCellID) const {
          VolumeID Sipm = static_cast<VolumeID>(_decoder->get(aCellID, fSipmId));
          return static_cast<int>(Sipm);
        }

        bool SCEPCALGridDRcalo::IsCerenkov(const CellID& aCellID) const {
            VolumeID isCeren = static_cast<VolumeID>(_decoder->get(aCellID, fIsCerenkovId));
            return static_cast<bool>(isCeren);
        }

        bool SCEPCALGridDRcalo::IsCerenkov(int eta, int phi) const {
          bool isCeren = false;
          if ( eta%2 == 1 ) { isCeren = !isCeren; }
          if ( phi%2 == 1 ) { isCeren = !isCeren; }
          return isCeren;
        }

        bool SCEPCALGridDRcalo::IsTower(const CellID& aCellID) const {
            VolumeID module = static_cast<VolumeID>(_decoder->get(aCellID, fModule));
            return module==0;
        }

        bool SCEPCALGridDRcalo::IsSiPM(const CellID& aCellID) const {
            VolumeID module = static_cast<VolumeID>(_decoder->get(aCellID, fModule));
            return module==1;
        }

        int SCEPCALGridDRcalo::getLast32bits(const CellID& aCellID) const {
            CellID aId64 = aCellID >> sizeof(int)*CHAR_BIT;
            int aId32 = (int)aId64;

            return aId32;
        }


        CellID SCEPCALGridDRcalo::convertLast32to64(const int aId32) const {
            CellID aId64 = (CellID)aId32;
            aId64 <<= sizeof(int)*CHAR_BIT;

            return aId64;
        }

/*        DRparamBase* SCEPCALGridDRcalo::setParamBase(int noEta) const {
            DRparamBase* paramBase = nullptr;

            if ( fParamEndcap->unsignedTowerNo(noEta) >= fParamBarrel->GetTotTowerNum() ) paramBase = static_cast<DRparamBase*>(fParamEndcap);
            else paramBase = static_cast<DRparamBase*>(fParamBarrel);

            if ( paramBase->GetCurrentTowerNum()==noEta ) return paramBase;

            // This should not be called while building detector geometry
            if (!paramBase->IsFinalized()) throw std::runtime_error("SCEPCALGridDRcalo::position should not be called while building detector geometry!");

            paramBase->SetDeltaThetaByTowerNo(noEta, fParamBarrel->GetTotTowerNum());
            paramBase->SetThetaOfCenterByTowerNo(noEta, fParamBarrel->GetTotTowerNum());
            paramBase->SetIsRHSByTowerNo(noEta);
            paramBase->SetCurrentTowerNum(noEta);
            paramBase->init();

            return paramBase;
        }*/


    }
}
