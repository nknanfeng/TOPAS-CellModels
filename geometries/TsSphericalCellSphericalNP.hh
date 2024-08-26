//
// ********************************************************************
// * This file is based on the TsSphericalCell example                *
// * from the TOPAS-nBio extensions to the TOPAS Simulation Toolkit.  *
// * The TOPAS-nBio extensions are freely available under the license *
// * agreement set forth at: https://topas-nbio.readthedocs.io/       *
// *                                                                  *              
// *  Extended by Marc B. Hahn (2021)                            *
// *  Please report bugs to hahn@physik.fu-berlin.de                 *
// *  or on https://github.com/BAMresearch/TOPAS-CellModels                            *    
// ********************************************************************
//
// A simple spherical cell with nanoparticles can be generated in a fast manner.
// The user has the option of including organelles: nucleus, mitochondria, cell membrane and/or nanoparticles.
// The user can add nanoparticles to the cytosol, to the surface of the nucleus and/or the mitochondria
// Up to 100000 objects can be created in a reasonable time. The time needed for generation of the geometries increases exponentially with the number of objects included in the cell.
// If you use this extension please cite the following literature:
// Hahn, M.B., Zutta Villate, J.M. "Combined cell and nanoparticle models for TOPAS to study radiation dose enhancement in cell organelles." Sci Rep 11, 6721 (2021).
// The extension is described in detail in https://doi.org/10.1038/s41598-021-85964-2

#ifndef TsSphericalCellSphericalNP_hh
#define TsSphericalCellSphericalNP_hh

#include "TsVGeometryComponent.hh"
#include <vector>


class TsSphericalCellSphericalNP : public TsVGeometryComponent
{    
public:
	TsSphericalCellSphericalNP(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsSphericalCellSphericalNP();
	
    // 构建几何模型的关键函数
	G4VPhysicalVolume* Construct();
    
    // 处理参数文件的关键函数
    void ResolveParameters();
    
    
private:
    
    // 数据成员
    G4double CellRadius; // 细胞半径
    G4double NucleusRadius; // 细胞核半径
    G4double MitoRadius; // 线粒体半径
    G4int MitoNumber; // 线粒体数量
    G4double MembraneThickness; // 细胞膜的厚度

    
    G4RotationMatrix* rotationMatrix;  // 旋转矩阵
    G4VPhysicalVolume* pNucleus; // 细胞核的Physics Volume
    
    // 存放了所有在细胞中生成的几何对象的坐标，包括细胞核、线粒体、纳米颗粒等
    // G4double用于表示三维空间中的一个点
    std::vector<std::vector<G4double> >  CellCoordinates; 

    std::vector<G4double>  tmpCoordinates;
    
    // 重叠检查的函数
    G4bool CheckOverlapOfSphereWithGeometryComponents(std::vector<std::vector<G4double> >& Coordinates, G4double r, G4double x, G4double y, G4double z);
    
    // 在细胞中添加球体（线粒体）
    G4ThreeVector* AddSphereToCell(G4double radius);
    // 在细胞中添加纳米颗粒
    G4ThreeVector* AddNanoparticleAtSphereSurface(G4double radius, G4int objectIndex);
    // 将坐标添加到坐标列表中
    void AddCoordinates(std::vector<std::vector<G4double> >& Coordinates, G4double r, G4double x, G4double y, G4double z);   

    
};

#endif
