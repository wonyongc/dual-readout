#ifndef GridDRcaloHandle_h
#define GridDRcaloHandle_h 1

#include "GridDRcalo.h"

#include "DD4hep/Segmentations.h"
#include "DD4hep/detail/SegmentationsInterna.h"

namespace dd4hep {
  class Segmentation;
  template <typename T> class SegmentationWrapper;

  typedef Handle<SegmentationWrapper<DDSegmentation::GridDRcalo>> GridDRcaloHandle;

  class GridDRcalo : public GridDRcaloHandle {
    typedef GridDRcaloHandle::Object Object;

  public:
    GridDRcalo() = default;
    GridDRcalo(const GridDRcalo& e) = default;
    GridDRcalo(const Segmentation& e) : Handle<Object>(e) {}
    GridDRcalo(const Handle<Object>& e) : Handle<Object>(e) {}
    template <typename Q>
    GridDRcalo(const Handle<Q>& e) : Handle<Object>(e) {}
    GridDRcalo& operator=(const GridDRcalo& seg) = default;
    bool operator==(const GridDRcalo& seg) const { return m_element == seg.m_element; }
    inline Position position(const CellID& id) const { return Position(access()->implementation->position(id)); }
    inline Position localPosition(const CellID& id) const { return Position(access()->implementation->localPosition(id)); }
    inline Position localPosition(int numx, int numy, int x_, int y_) const { return Position(access()->implementation->localPosition(numx,numy,x_,y_)); }

    inline dd4hep::CellID cellID(const Position& local, const Position& global, const VolumeID& volID) const {
      return access()->implementation->cellID(local, global, volID);
    }

    inline VolumeID setVolumeID(int numSystem, int numEta, int numPhi, int numDepth) const { return access()->implementation->setVolumeID(numSystem, numEta,numPhi, numDepth); }
    inline CellID setCellID(int numEta, int numPhi, int x, int y) const { return access()->implementation->setCellID(numEta, numPhi, x, y); }

    inline void setGridSize(double grid) { access()->implementation->setGridSize(grid); }

    inline int numEta(const CellID& aCellID) const { return access()->implementation->numEta(aCellID); }
    inline int numPhi(const CellID& aCellID) const { return access()->implementation->numPhi(aCellID); }

    inline int numX(const CellID& aCellID) const { return access()->implementation->numX(aCellID); }
    inline int numY(const CellID& aCellID) const { return access()->implementation->numY(aCellID); }

    inline int x(const CellID& aCellID) const { return access()->implementation->x(aCellID); } // approx eta direction
    inline int y(const CellID& aCellID) const { return access()->implementation->y(aCellID); } // approx phi direction

    inline bool IsCerenkov(const CellID& aCellID) const { return access()->implementation->IsCerenkov(aCellID); }
    inline bool IsCerenkov(int col, int row) const { return access()->implementation->IsCerenkov(col, row); }

    inline bool IsTower(const CellID& aCellID) const { return access()->implementation->IsTower(aCellID); }
    inline bool IsSiPM(const CellID& aCellID) const { return access()->implementation->IsSiPM(aCellID); }

    inline int getFirst32bits(const CellID& aCellID) const { return access()->implementation->getFirst32bits(aCellID); }
    inline int getLast32bits(const CellID& aCellID) const { return access()->implementation->getLast32bits(aCellID); }

    inline CellID convertFirst32to64(const int aId32) const { return access()->implementation->convertFirst32to64(aId32); }
    inline CellID convertLast32to64(const int aId32) const { return access()->implementation->convertLast32to64(aId32); }

    inline int numEta(const int& aId32) const { return access()->implementation->numEta(aId32); }
    inline int numPhi(const int& aId32) const { return access()->implementation->numPhi(aId32); }

    inline int numX(const int& aId32) const { return access()->implementation->numX(aId32); }
    inline int numY(const int& aId32) const { return access()->implementation->numY(aId32); }

    inline int x(const int& aId32) const { return access()->implementation->x(aId32); }
    inline int y(const int& aId32) const { return access()->implementation->y(aId32); }

    inline bool IsCerenkov(const int& aId32) const { return access()->implementation->IsCerenkov(aId32); }

    inline bool IsTower(const int& aId32) const { return access()->implementation->IsTower(aId32); }
    inline bool IsSiPM(const int& aId32) const { return access()->implementation->IsSiPM(aId32); }
  };
}

#endif
