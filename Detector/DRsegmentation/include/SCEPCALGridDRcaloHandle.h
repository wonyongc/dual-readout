#ifndef SCEPCALGridDRcaloHandle_h
#define SCEPCALGridDRcaloHandle_h 1

#include "SCEPCALGridDRcalo.h"

#include "DD4hep/Segmentations.h"
#include "DD4hep/detail/SegmentationsInterna.h"

namespace dd4hep {
    class Segmentation;
    template <typename T> class SegmentationWrapper;

    typedef Handle<SegmentationWrapper<DDSegmentation::SCEPCALGridDRcalo>> SCEPCALGridDRcaloHandle;

    class SCEPCALGridDRcalo : public SCEPCALGridDRcaloHandle {
        typedef SCEPCALGridDRcaloHandle::Object Object;

    public:
        SCEPCALGridDRcalo() = default;
        SCEPCALGridDRcalo(const SCEPCALGridDRcalo& e) = default;
        SCEPCALGridDRcalo(const Segmentation& e) : Handle<Object>(e) {}
        SCEPCALGridDRcalo(const Handle<Object>& e) : Handle<Object>(e) {}
        template <typename Q>
        SCEPCALGridDRcalo(const Handle<Q>& e) : Handle<Object>(e) {}
        SCEPCALGridDRcalo& operator=(const SCEPCALGridDRcalo& seg) = default;
        bool operator==(const SCEPCALGridDRcalo& seg) const { return m_element == seg.m_element; }

        inline Position position(const CellID& id) const {
            return Position(access()->implementation->position(id));
            }

        inline dd4hep::CellID cellID(const Position& local, const Position& global, const VolumeID& volID) const {
            return access()->implementation->cellID(local, global, volID);
        }

        inline VolumeID setVolumeID(int System, int Eta, int Phi, int Depth) const { return access()->implementation->setVolumeID(System, Eta,Phi,Depth); }

        inline CellID setCellID(int System, int Eta, int Phi, int Depth) const { return access()->implementation->setCellID(System, Eta, Phi, Depth); }

        inline int System(const CellID& aCellID) const { return access()->implementation->System(aCellID); }
        inline int Eta(const CellID& aCellID) const { return access()->implementation->Eta(aCellID); }
        inline int Phi(const CellID& aCellID) const { return access()->implementation->Phi(aCellID); }
        inline int Depth(const CellID& aCellID) const { return access()->implementation->Depth(aCellID); }
        inline int Sipm(const CellID& aCellID) const { return access()->implementation->Sipm(aCellID); }

        inline bool IsCerenkov(const CellID& aCellID) const { return access()->implementation->IsCerenkov(aCellID); }
        inline bool IsCerenkov(int eta, int phi) const { return access()->implementation->IsCerenkov(eta, phi); }

        inline bool IsTower(const CellID& aCellID) const { return access()->implementation->IsTower(aCellID); }
        inline bool IsSiPM(const CellID& aCellID) const { return access()->implementation->IsSiPM(aCellID); }

        inline int getFirst32bits(const CellID& aCellID) const { return access()->implementation->getFirst32bits(aCellID); }
        inline int getLast32bits(const CellID& aCellID) const { return access()->implementation->getLast32bits(aCellID); }

        inline CellID convertFirst32to64(const int aId32) const { return access()->implementation->convertFirst32to64(aId32); }
        inline CellID convertLast32to64(const int aId32) const { return access()->implementation->convertLast32to64(aId32); }

        inline int System(const int& aId32) const { return access()->implementation->System(aId32); }
        inline int Eta(const int& aId32) const { return access()->implementation->Eta(aId32); }
        inline int Phi(const int& aId32) const { return access()->implementation->Phi(aId32); }
        inline int Depth(const int& aId32) const { return access()->implementation->Depth(aId32); }
        inline int Sipm(const int& aId32) const { return access()->implementation->Sipm(aId32); }
        inline bool IsCerenkov(const int& aId32) const { return access()->implementation->IsCerenkov(aId32); }

        inline bool IsTower(const int& aId32) const { return access()->implementation->IsTower(aId32); }
        inline bool IsSiPM(const int& aId32) const { return access()->implementation->IsSiPM(aId32); }

    };
}

#endif
