#include "DD4hep/Factories.h"
#include "DD4hep/detail/SegmentationsInterna.h"


using namespace dd4hep;
using namespace dd4hep::DDSegmentation;

namespace {
  template <typename T> dd4hep::SegmentationObject*
  create_segmentation(const dd4hep::BitFieldCoder* decoder) {
  return new dd4hep::SegmentationWrapper<T>(decoder);
  }
}

#include "GridDRcalo.h"
DECLARE_SEGMENTATION(GridDRcalo, create_segmentation<dd4hep::DDSegmentation::GridDRcalo>)

#include "SCEPCALGridDRcalo.h"
DECLARE_SEGMENTATION(SCEPCALGridDRcalo, create_segmentation<dd4hep::DDSegmentation::SCEPCALGridDRcalo>)
