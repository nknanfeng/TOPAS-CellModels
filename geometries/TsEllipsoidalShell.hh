// TsEllipsoidalShell.hh

#ifndef TSELLIPSOIDALSHELL_HH
#define TSELLIPSOIDALSHELL_HH

#include "TsVGeometryComponent.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"

class TsEllipsoidalShell : public TsVGeometryComponent {
public:
    TsEllipsoidalShell(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
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
