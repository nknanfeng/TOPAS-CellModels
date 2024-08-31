// TsEllipsoidalShell.cc

#include "TsEllipsoidalShell.hh"
#include "G4SystemOfUnits.hh"
#include "G4Ellipsoid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"

TsEllipsoidalShell::TsEllipsoidalShell(const G4String& name, TsParameterManager* pM, TsMaterialManager* pMM, TsGeometryManager* pGM, TsVGeometryComponent* parentComponent)
    : TsVGeometryComponent(name, pM, pMM, pGM, parentComponent)
{
    // 获取外层椭球的三个半长轴
    fOuterX = fPm->GetDoubleParameter(GetFullParmName("OuterX"), "Length");
    fOuterY = fPm->GetDoubleParameter(GetFullParmName("OuterY"), "Length");
    fOuterZ = fPm->GetDoubleParameter(GetFullParmName("OuterZ"), "Length");

    // 获取内层椭球的三个半长轴
    fInnerX = fPm->GetDoubleParameter(GetFullParmName("InnerX"), "Length");
    fInnerY = fPm->GetDoubleParameter(GetFullParmName("InnerY"), "Length");
    fInnerZ = fPm->GetDoubleParameter(GetFullParmName("InnerZ"), "Length");
}

G4VPhysicalVolume* TsEllipsoidalShell::Construct()
{
    BeginConstruction();

    // 创建外层椭球
    G4Ellipsoid* outerEllipsoid = new G4Ellipsoid(fName, fOuterX, fOuterY, fOuterZ);

    // 创建内层椭球
    G4Ellipsoid* innerEllipsoid = new G4Ellipsoid(fName, fInnerX, fInnerY, fInnerZ);

    // 创建椭球壳（外层椭球减去内层椭球）
    G4SubtractionSolid* ellipsoidalShell = new G4SubtractionSolid("EllipsoidalShell", outerEllipsoid, innerEllipsoid);

    // 创建逻辑体
    G4LogicalVolume* logicEllipsoidalShell = CreateLogicalVolume(ellipsoidalShell);

    // 放置物理体
    G4VPhysicalVolume* physEllipsoidalShell = CreatePhysicalVolume(logicEllipsoidalShell);

    InstantiateChildren(physEllipsoidalShell);

    return physEllipsoidalShell;
}
