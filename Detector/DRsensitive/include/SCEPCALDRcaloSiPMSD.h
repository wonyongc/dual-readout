#ifndef SCEPCALDRcaloSiPMSD_h
#define SCEPCALDRcaloSiPMSD_h 1

#include "DRcaloSiPMHit.h"
#include "SCEPCALGridDRcalo.h"
#include "DD4hep/Segmentations.h"

#include "G4VSensitiveDetector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"

namespace drc {
    class SCEPCALDRcaloSiPMSD : public G4VSensitiveDetector {
    public:
        SCEPCALDRcaloSiPMSD(const std::string aName, const std::string aReadoutName, const dd4hep::Segmentation& aSeg);
        ~SCEPCALDRcaloSiPMSD();

        virtual void Initialize(G4HCofThisEvent* HCE) final;
        virtual bool ProcessHits(G4Step* aStep, G4TouchableHistory*) final;

    private:
        DRcaloSiPMHitsCollection* fHitCollection;
        dd4hep::DDSegmentation::SCEPCALGridDRcalo* fSeg;
        G4int fHCID;

        G4int fWavBin;
        G4int fTimeBin;
        G4float fWavlenStart;
        G4float fWavlenEnd;
        G4float fTimeStart;
        G4float fTimeEnd;
        G4float fWavlenStep;
        G4float fTimeStep;

        G4double wavToE(G4double wav) { return h_Planck*c_light/wav; }

        float findWavCenter(G4double en);
        float findTimeCenter(G4double stepTime);
    };
}

#endif
