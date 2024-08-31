// Component for TsEllipsoidalShell 
//
// ********************************************************************
// * This file is based on the TsSphericalCell example                *
// * from the TOPAS-nBio extensions to the TOPAS Simulation Toolkit.  *
// * The TOPAS-nBio extensions are freely available under the license *
// * agreement set forth at: https://topas-nbio.readthedocs.io/       *
// *                                                                  *              
// *  Extended by nknanfeng (2024)                                 *
// *  Please report bugs to hahn@physik.fu-berlin.de                  *
// *  or on https://github.com/BAMresearch/TOPAS-CellModels           *    
// ********************************************************************
//
// A simple spherical cell with nanoparticles can be generated in a fast manner.
// The user has the option of including organelles: nucleus, mitochondria, cell membrane and/or nanoparticles.
// The user can add nanoparticles to the cytosol, to the surface of the nucleus and/or the mitochondria
// Up to 100000 objects can be created in a reasonable time. The time needed for generation of the geometries increases exponentially with the number of objects included in the cell.
// If you use this extension please cite the following literature:
// Hahn, M.B., Zutta Villate, J.M. "Combined cell and nanoparticle models for TOPAS to study radiation dose enhancement in cell organelles." Sci Rep 11, 6721 (2021).
// The extension is described in detail in https://doi.org/10.1038/s41598-021-85964-2

#include "TsEllipsoidalShell.hh"
#include "G4SystemOfUnits.hh"
#include "G4Ellipsoid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "TsParameterManager.hh"

TsEllipsoidalShell::TsEllipsoidalShell(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM, TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
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
    fEnvelopeLog = CreateLogicalVolume(ellipsoidalShell);

    // 放置物理体
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

    InstantiateChildren(fEnvelopePhys);
    
	return fEnvelopePhys;
}
