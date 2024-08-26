// Component for TsSphericalCellSphericalNP
//
// ********************************************************************
// * This file is based on the TsSphericalCell example                *
// * from the TOPAS-nBio extensions to the TOPAS Simulation Toolkit.  *
// * The TOPAS-nBio extensions are freely available under the license *
// * agreement set forth at: https://topas-nbio.readthedocs.io/       *
// *                                                                  *              
// *  Extended by Marc B. Hahn (2021)                                 *
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

#include "TsSphericalCellSphericalNP.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

// 构造函数
TsSphericalCellSphericalNP::TsSphericalCellSphericalNP(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM, TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    // 读取参数文件中的参数
    ResolveParameters();
    
    // 初始化存放这些几何
    tmpCoordinates.resize(4);
    CellCoordinates.resize(0);
   
}


TsSphericalCellSphericalNP::~TsSphericalCellSphericalNP()
{;}

void TsSphericalCellSphericalNP::ResolveParameters() {
    
    // 这一步骤中，只读取了细胞半径，fPm是TsParameterManager的一个对象
    // todo：改成椭球形细胞的话，这里是需要进行修改的
    CellRadius = fPm->GetDoubleParameter(GetFullParmName("CellRadius"), "Length");
    
}


G4VPhysicalVolume* TsSphericalCellSphericalNP::Construct()
{
    // 通常用于启动一个几何体或组件的构建过程
	BeginConstruction();
    
    //***********************************************************************
    //              Envelope Geometry : Spherical cell
    //***********************************************************************
    
    // 构造一个球形的细胞来作为envelop, solid volume,
    // todo：这里需要改成椭球形细胞
    G4Sphere* gCell = new G4Sphere (fName, 0.0, CellRadius, 0., CLHEP::twopi, 0., CLHEP::pi);
    rotationMatrix = new G4RotationMatrix();

    // 创建log volume 和 physics volume, 这两个概念源于geant4中
    fEnvelopeLog = CreateLogicalVolume(gCell);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

    
    //***********************************************************************
    // Optional : include a membrane, nucleus, mitochondria and/or nanoparticles in the cell
    //***********************************************************************
    
    //***************************
    // Subcomponent: Membrane 
    //***************************
    
    // 细胞膜厚度
    G4String nameMembrane = GetFullParmName("Membrane/Thickness");
    if (fPm->ParameterExists(nameMembrane)) {
        
        
        //Membrane thickness for scoring
        MembraneThickness  = fPm->GetDoubleParameter( nameMembrane, "Length" );
        G4ThreeVector* CellPosition = new G4ThreeVector(0,0,0);

        // 球形的细胞膜及相应的solid volume, logic volume 和 physics volume
        // todo:这里需要修改成椭球形
        G4Sphere* gMembrane = new G4Sphere ("Membrane", CellRadius-MembraneThickness, CellRadius, 0., CLHEP::twopi, 0., CLHEP::pi);
        G4LogicalVolume* lMembrane = CreateLogicalVolume("Membrane", gMembrane);
        G4VPhysicalVolume* pMembrane = CreatePhysicalVolume("Membrane", lMembrane, rotationMatrix, CellPosition, fEnvelopePhys);
        
    }
    
    
    //***************************
    // Subcomponent: Nucleus
    //***************************
    
    NucleusRadius = 0.0*um;
    G4String name = GetFullParmName("Nucleus/NucleusRadius");
    if (fPm->ParameterExists(name)) {
        
        NucleusRadius = fPm->GetDoubleParameter(name, "Length");
             
        G4double transNucX = 0 * um;
        G4double transNucY = 0 * um;
        G4double transNucZ = 0 * um;
        
        // 细胞核，x,y,z方向上的偏移量
        G4String name1 = GetFullParmName("Nucleus/translateX");
        G4String name2 = GetFullParmName("Nucleus/translateY");
        G4String name3 = GetFullParmName("Nucleus/translateZ");
        
        if (fPm -> ParameterExists(name1)){
            transNucX = fPm->GetDoubleParameter(name1, "Length");
        }    
        if (fPm -> ParameterExists(name2)){
            transNucY = fPm->GetDoubleParameter(name2, "Length");
        }
        if (fPm -> ParameterExists(name3)){
            transNucZ = fPm->GetDoubleParameter(name3, "Length");
        }
        
        // 简单验证，细胞核中心不能在细胞外
        if ((sqrt(transNucX*transNucX)+(transNucY*transNucY)+(transNucZ*transNucZ)) > (CellRadius-NucleusRadius)) {
                G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
                G4cerr << "Parameter " << name1 << " sets nucleus outside of cell." << G4endl;
                exit(1);
            }
            
        // 细胞核的solid volume, logical volume及 physics volume
        G4ThreeVector* NucPos = new G4ThreeVector(transNucX,transNucY,transNucZ);
        G4Sphere* gNucleus = new G4Sphere ("Nucleus", 0.0, NucleusRadius, 0., CLHEP::twopi, 0., CLHEP::pi);
        G4LogicalVolume* lNucleus = CreateLogicalVolume("Nucleus", gNucleus);
        pNucleus = CreatePhysicalVolume("Nucleus", lNucleus, rotationMatrix, NucPos, fEnvelopePhys);
        
        AddCoordinates(CellCoordinates,NucleusRadius,transNucX,transNucY,transNucX);
        
    //*******************************
    // Subcomponent: Nanoparticles at the Nucleus
    //*******************************
    
    // 细胞核中的纳米粒子
    G4String nameNPNuc = GetFullParmName("Nanoparticle/NumberOfNanoparticlesAtNucleus");
    if (fPm->ParameterExists(nameNPNuc)) {
        
        //number of nanooparticles at nucleus
        const G4int NbOfNP  = fPm->GetIntegerParameter( GetFullParmName("Nanoparticle/NumberOfNanoparticlesAtNucleus") );
        
        //radius of the nanoparticles (default values if none are specified)
        // 纳米粒子的默认半径
        G4double rNP = 10*nanometer;
            
        G4String nameNPNucR=GetFullParmName("Nanoparticle/r");
        if (fPm->ParameterExists(nameNPNucR)){
            rNP = fPm->GetDoubleParameter(GetFullParmName("Nanoparticle/r"), "Length" );
        }
        
        G4Orb* gNP = new G4Orb("Nanoparticle", rNP);
        G4LogicalVolume* lNP = CreateLogicalVolume("Nanoparticle", gNP);
        
        //Randomly distribute mitochondria throughout the cell volume
        for (int m = 0; m < NbOfNP; m++){
            
            G4cout << "** Add NP at Nucleus  " << m  <<  " **" << G4endl;

            // 在AddNanoparticleAtSphereSurface(rNP, 0)函数中实现在指定半径的球体表面上随机放置一个不重合的纳米颗粒
            G4VPhysicalVolume* pNP = CreatePhysicalVolume("Nanoparticle", m, true, lNP, rotationMatrix, AddNanoparticleAtSphereSurface(rNP, 0), fEnvelopePhys);
                      
        }
     }
   }
    

    
    //*******************************
    // Subcomponent: Mitochondria
    //*******************************
    MitoNumber = 0;
    // 线粒体
    G4String name6 = GetFullParmName("Mitochondria/NumberOfMitochondria");
    if (fPm->ParameterExists(name6)) {
        
        //number of mitochondria
        MitoNumber  = fPm->GetIntegerParameter( GetFullParmName("Mitochondria/NumberOfMitochondria") );
        
        //radius of the mitochondria (default values if none are specified)
        // 线粒体默认半径
        G4double MitoRadius = 0.5*micrometer;
        
        G4String name7=GetFullParmName("Mitochondria/r");
        if (fPm->ParameterExists(name7)){MitoRadius = fPm->GetDoubleParameter(GetFullParmName("Mitochondria/r"), "Length" );}
    
        G4Orb* gMito = new G4Orb("Mitochondria", MitoRadius);
        G4LogicalVolume* lMito = CreateLogicalVolume("Mitochondria", gMito);
        
        //Randomly distribute mitochondria throughout the cell volume
        for (int k = 0; k < MitoNumber; k++){
            // AddSphereToCell(MitoRadius)，在细胞内部添加球体（线粒体）确保不会与已有的几何组件重合
            G4VPhysicalVolume* pMito = CreatePhysicalVolume("Mitochondria", k, true, lMito, rotationMatrix, AddSphereToCell(MitoRadius), fEnvelopePhys);
           
           
        }
        
        //*******************************
        // Subcomponent: Nanoparticles at the Mitochondria
        //*******************************
    
        // 线粒体表面的纳米颗粒
        G4String name8 = GetFullParmName("Nanoparticle/NumberOfNanoparticlesAtMitochondria");
        if (fPm->ParameterExists(name8)) {
        
            //number of mitochondria
            const G4int NbOfNP  = fPm->GetIntegerParameter( GetFullParmName("Nanoparticle/NumberOfNanoparticlesAtMitochondria") );
        
            //radius of the nanoparticles (default values if none are specified)
            G4double rNP = 10*nanometer;
            
            name=GetFullParmName("Nanoparticle/r");
            if (fPm->ParameterExists(name)){
                rNP = fPm->GetDoubleParameter(GetFullParmName("Nanoparticle/r"), "Length" );
            }
        
            G4Orb* gNP = new G4Orb("Nanoparticle", rNP);
            G4LogicalVolume* lNP = CreateLogicalVolume("Nanoparticle", gNP);
        
            //Randomly distribute mitochondria throughout the cell volume
            for (int m = 0; m < NbOfNP; m++){
            
                G4int NPatMitochondria = std::round((G4UniformRand()*MitoNumber) + 1.0);
                G4cout << "** Add NP "<<  m << " at Mitochondria  " << NPatMitochondria  <<  " **" << G4endl;
                
                // 在线粒体表面添加纳米颗粒
                G4VPhysicalVolume* pNP = CreatePhysicalVolume("Nanoparticle", m, true, lNP, rotationMatrix, AddNanoparticleAtSphereSurface(rNP, NPatMitochondria), fEnvelopePhys);
                      
            }
        }
    }
    
    
    //*******************************
    // Subcomponent: Nanoparticles
    //*******************************
    
    // 在细胞中的纳米颗粒
    G4String name9 = GetFullParmName("Nanoparticle/NumberOfNanoparticles");
    if (fPm->ParameterExists(name9)) {
        
        //number of nanoparticles
        const G4int NbOfNP  = fPm->GetIntegerParameter( GetFullParmName("Nanoparticle/NumberOfNanoparticles") );
        
        //radius of the nanoparticles (default values if none are specified)
        G4double rNP = 10*nanometer;
            
        G4String name10 =GetFullParmName("Nanoparticle/r");
        if (fPm->ParameterExists(name10)){
            rNP = fPm->GetDoubleParameter(GetFullParmName("Nanoparticle/r"), "Length" );
        }
        
        G4Orb* gNP = new G4Orb("Nanoparticle", rNP);
        G4LogicalVolume* lNP = CreateLogicalVolume("Nanoparticle", gNP);
        
        //Randomly distribute mitochondria throughout the cell volume
        for (int m = 0; m < NbOfNP; m++){
            
            G4cout << "** Add NP  " << m  <<  " **" << G4endl;

            G4VPhysicalVolume* pNP = CreatePhysicalVolume("Nanoparticle", m, true, lNP, rotationMatrix, AddSphereToCell(rNP), fEnvelopePhys);
        }
    }

    G4cout << "*** Total objects in cell  : " << CellCoordinates.size() <<" ***" <<G4endl;
    InstantiateChildren(fEnvelopePhys);
	return fEnvelopePhys;
}


G4ThreeVector* TsSphericalCellSphericalNP::AddSphereToCell(G4double radius){
    
    long unsigned placementAttempts = 0;
    long unsigned placementAttemptsWarning = 10000;
    G4double distanceToMembrane = CellRadius-radius-MembraneThickness;

    while (true){
                    
        G4double u = G4UniformRand()*2*pi;
        G4double v = std::acos(2*G4UniformRand()-1);
        G4double distance = G4UniformRand()*(distanceToMembrane);
                                    
        G4double x =  distance * std::cos(u) * std::sin(v);
        G4double y =  distance * std::sin(u) * std::sin(v);
        G4double z =  distance * std::cos(v);
        
        if (CheckOverlapOfSphereWithGeometryComponents(CellCoordinates, radius,x,y,z)){
            
            placementAttempts ++;
            if (placementAttempts > placementAttemptsWarning){
                G4cerr << "Couldn't find a proper placement position for the current object within "<< placementAttempts <<" attempts. Continuing..."<<  G4endl;
                placementAttemptsWarning = 2*placementAttemptsWarning;
            }
        }
                
        else{
            
            G4ThreeVector* position = new G4ThreeVector(x,y,z);
            AddCoordinates(CellCoordinates,radius,x,y,z);
            return position;
        }
   }
}


G4bool TsSphericalCellSphericalNP::CheckOverlapOfSphereWithGeometryComponents(std::vector<std::vector<G4double> >& Coordinates,G4double r, G4double x, G4double y, G4double z){
    
    for(int i=0; i<Coordinates.size(); i++){
        
        if ( (r+Coordinates[i][0] ) > sqrt( ((x-Coordinates[i][1])*(x-Coordinates[i][1])) + ((y-Coordinates[i][2])*(y-Coordinates[i][2])) + ((z-Coordinates[i][3])*(z-Coordinates[i][3])) ) ) { 
        return true; 
        }
    }
    return false;    
}


// 记录所有这些几何体的半径和位置 
// todo:看起来仅仅是针对球体的，细胞和细胞膜并未被记录
void TsSphericalCellSphericalNP::AddCoordinates(std::vector<std::vector<G4double> >& Coordinates, G4double r, G4double x, G4double y, G4double z){
    tmpCoordinates[0]=r;
    tmpCoordinates[1]=x;
    tmpCoordinates[2]=y;
    tmpCoordinates[3]=z;
    Coordinates.push_back(tmpCoordinates);
}


G4ThreeVector* TsSphericalCellSphericalNP::AddNanoparticleAtSphereSurface(G4double radius, G4int objectIndex){
    long unsigned placementAttemps = 0;
    long unsigned placementAttempsWarning = 10000;
    G4double distance = CellCoordinates[objectIndex][0] + radius *1.01;

    while (true){
                    
        G4double u = G4UniformRand()*2*pi;
        G4double v = std::acos(2*G4UniformRand()-1);
                                    
        G4double x =  CellCoordinates[objectIndex][1] +   distance * std::cos(u) * std::sin(v);
        G4double y =  CellCoordinates[objectIndex][2] +   distance * std::sin(u) * std::sin(v);
        G4double z =  CellCoordinates[objectIndex][3] +   distance * std::cos(v);
        
        if (CheckOverlapOfSphereWithGeometryComponents(CellCoordinates, radius,x,y,z)){
            
            placementAttemps ++;
            if (placementAttemps > placementAttempsWarning){
                G4cerr << "Couldn't find a proper placement position for the current object within "<< placementAttemps <<" attemps. Continuing..."<<  G4endl;
                placementAttempsWarning = 2*placementAttempsWarning;
            }
            if (placementAttemps > 64000){
                objectIndex = std::round((G4UniformRand()*MitoNumber) + 1.0);
                G4cerr << "Couldn't find a proper placement position at the current mitochondria within "<< placementAttemps <<" attemps. Try mitochondria "<< objectIndex << "instead. "<<  G4endl;
                placementAttemps = 0;
            }
        }
                
        else{
            
            G4ThreeVector* position = new G4ThreeVector(x,y,z);
            AddCoordinates(CellCoordinates,radius,x,y,z);
            return position;
        }
    }
}
