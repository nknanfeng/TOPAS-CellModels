// TsEllipsoidalShell.hh

#ifndef TSELLIPSOIDALSHELL_HH
#define TSELLIPSOIDALSHELL_HH

#include "TsVGeometryComponent.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"

class TsEllipsoidalShell : public TsVGeometryComponent {
public:
    TsEllipsoidalShell(const G4String& name, TsParameterManager* pM, TsMaterialManager* pMM, TsGeometryManager* pGM, TsVGeometryComponent* parentComponent);
    virtual ~TsEllipsoidalShell() {}

    virtual G4VPhysicalVolume* Construct();

protected:
    G4double fOuterX;
    G4double fOuterY;
    G4double fOuterZ;
    G4double fInnerX;
    G4double fInnerY;
    G4double fInnerZ;
};

#endif
