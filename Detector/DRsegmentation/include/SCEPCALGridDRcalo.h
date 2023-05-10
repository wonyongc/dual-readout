#ifndef SCEPCALGridDRcalo_h
#define SCEPCALGridDRcalo_h 1


#include "DDSegmentation/Segmentation.h"

#include "TVector3.h"
#include "DD4hep/DetFactoryHelper.h"

#include <vector>
#include <cmath>

namespace dd4hep {
namespace DDSegmentation {
class SCEPCALGridDRcalo : public Segmentation {
        public:
            SCEPCALGridDRcalo(const std::string& aCellEncoding);
            SCEPCALGridDRcalo(const BitFieldCoder* decoder);
            virtual ~SCEPCALGridDRcalo() override;

            virtual Vector3D position(const CellID& aCellID) const;

            virtual CellID cellID(const Vector3D& aLocalPosition, const Vector3D& aGlobalPosition,
                                  const VolumeID& aVolumeID) const;

            VolumeID setVolumeID(int System, int Eta, int Phi, int Depth) const;
            CellID setCellID(int System, int Eta, int Phi, int Depth) const;

            int System(const CellID& aCellID) const;
            int Eta(const CellID& aCellID) const;
            int Phi(const CellID& aCellID) const;
            int Depth(const CellID& aCellID) const;
            int Sipm(const CellID& aCellID) const;

            bool IsCerenkov(const CellID& aCellID) const;
            bool IsCerenkov(int eta, int phi) const;

            bool IsTower(const CellID& aCellID) const;
            bool IsSiPM(const CellID& aCellID) const;

            int getFirst32bits(const CellID& aCellID) const { return (int)aCellID; }
            int getLast32bits(const CellID& aCellID) const;

            CellID convertFirst32to64(const int aId32) const { return (CellID)aId32; }
            CellID convertLast32to64(const int aId32) const;

            int System(const int& aId32) const { return System( convertFirst32to64(aId32) ); }
            int Eta(const int& aId32) const { return Eta( convertFirst32to64(aId32) ); }
            int Phi(const int& aId32) const { return Phi( convertFirst32to64(aId32) ); }
            int Depth(const int& aId32) const { return Depth( convertFirst32to64(aId32) ); }
            int Sipm(const int& aId32) const { return Sipm( convertFirst32to64(aId32) ); }

            bool IsCerenkov(const int& aId32) const { return IsCerenkov( convertFirst32to64(aId32) ); }

            bool IsTower(const int& aId32) const { return IsTower( convertFirst32to64(aId32) ); }
            bool IsSiPM(const int& aId32) const { return IsSiPM( convertFirst32to64(aId32) ); }

            inline const std::string& fieldNameEta() const { return fEtaId; }
            inline const std::string& fieldNamePhi() const { return fPhiId; }
            inline const std::string& fieldNameDepth() const { return fDepthId; }
            inline const std::string& fieldNameSipm() const { return fSipmId; }
            inline const std::string& fieldNameIsCerenkov() const { return fIsCerenkovId; }
            inline const std::string& fieldNameModule() const { return fModule; }

            inline void setFieldNameEta(const std::string& fieldName) { fEtaId = fieldName; }
            inline void setFieldNamePhi(const std::string& fieldName) { fPhiId = fieldName; }
            inline void setFieldNameDepth(const std::string& fieldName) { fDepthId = fieldName; }
            inline void setFieldNameSipm(const std::string& fieldName) { fSipmId = fieldName; }
            inline void setFieldNameIsCerenkov(const std::string& fieldName) { fIsCerenkovId = fieldName; }
            inline void setFieldNameModule(const std::string& fieldName) { fModule = fieldName; }


        protected:
            std::string fSystemId;
            std::string fEtaId;
            std::string fPhiId;
            std::string fDepthId;
            std::string fSipmId;
            std::string fIsCerenkovId;
            std::string fModule;

        private:

        };
    }
}

#endif
