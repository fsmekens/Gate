/*----------------------
  Copyright (C): OpenGATE Collaboration

  This software is distributed under the terms of the GNU Lesser General  Public Licence (LGPL)
  See GATE/LICENSE.txt for further details
  ----------------------*/

// This actor is only compiled if ITK is available
#include "GateConfiguration.h"
#ifdef  GATE_USE_ITK

/*
  \class GateThermalActor
  \author vesna.cuplov@gmail.com
  \author fsmekens@gmail.com
  \brief Class GateThermalActor : This actor produces voxelised images of the heat diffusion in tissue.

                                                                    absorption map        heat diffusion map
         laser                _________                               _________              _________
         optical photons     |         |                             |         |            |         |
         ~~~~~~~>            |         |    GATE                     |         |            |   xx    |
         ~~~~~~~>            | phantom |    Simulation Results ==>   |   xx    |     +      |  xxxx   |
         ~~~~~~~>            |         |    (voxelised images)       |   xx    |            |  xxxx   |
         ~~~~~~~>            |         |                             |         |            |   xx    |
                             |_________|                             |_________|            |_________|


  Parameters of the simulation given by the User in the macro:
	- setDiffusivity: tissue thermal diffusivity in mm2/s
	- setTime: diffusion time in s
	- setBloodPerfusionRate: blood perfusion rate in s-1 for the advection term
	- setBloodDensity: blood density (kg/m3)
	- setBloodHeatCapacity: blood heat capacity kJ/(kg C)
	- setTissueDensity: tissue density (kg/m3)
	- setTissueHeatCapacity: tissue heat capacity kJ/(kg C)
	- OPTIONAL: setSimulationScale to get more fluence.
*/

#include <G4VoxelLimits.hh>
#include <G4NistManager.hh>

#include "GateFSThermalActor.hh"
#include "GateMiscFunctions.hh"
#include "G4VProcess.hh"
#include "GateMHDImage.hh"
#include "GateImageT.hh"
#include "GateMiscFunctions.hh"
#include "GateMachine.hh"
#include "GateApplicationMgr.hh"
#include <sys/time.h>
#include <iostream>
#include <string>

//-----------------------------------------------------------------------------

GateFSThermalActor::GateFSThermalActor(G4String name, G4int depth):
  GateVImageActor(name,depth) {
  GateDebugMessageInc("Actor",4,"GateFSThermalActor() -- begin"<<G4endl);

  mCurrentEvent=-1;
  mUserRelaxationTime = -1.0;
  mIsDiffusionActivated = false;
  mIsPerfusionByMaterial = false;
  mCropSize = 5;
  
  mIsPerfusionActivated = false;
  mIsPerfusionByMaterial = false;
  mIsPerfusionByConstant = false;
  mIsPerfusionByImage = false;
  mPerfusionRatio = 0.99;
  mUserPerfusionImageName = "";
  mUserBloodPerfusionRate = 0.0;
  mUserBloodDensity = -1.0;
  mUserBloodHeatCapacity = -1.0;
  mUserTissueHeatCapacity = -1.0;
  mMeasurementFilename = "";
  mIsMeasurementActivated = false;

  mITKfloatReaderFilter = FloatReaderType::New();
  mITKdoubleReaderFilter = DoubleReaderType::New();
  mITKfloatDuplicatorFilter = FloatDuplicatorType::New();
  mITKdoubleDuplicatorFilter = DoubleDuplicatorType::New();
  mITKmultiplyImageFilter = MultiplyFilterType::New();
  mITKaddImageFilter = AddImageFilterType::New();
  mITKsubtractImageFilter = SubtractImageFilterType::New();
  mITKmaxImageFilter = MaxImageFilterType::New();
  mITKgaussianFilterX = GaussianFilterType::New();
  mITKgaussianFilterY = GaussianFilterType::New();
  mITKgaussianFilterZ = GaussianFilterType::New();
  mITKgaussianFilterX->SetDirection(0);
  mITKgaussianFilterY->SetDirection(1);
  mITKgaussianFilterZ->SetDirection(2);
  mITKgaussianFilterX->SetOrder(GaussianFilterType::ZeroOrder);
  mITKgaussianFilterY->SetOrder(GaussianFilterType::ZeroOrder);
  mITKgaussianFilterZ->SetOrder(GaussianFilterType::ZeroOrder);
  mITKgaussianFilterX->SetNormalizeAcrossScale(false);
  mITKgaussianFilterY->SetNormalizeAcrossScale(false);
  mITKgaussianFilterZ->SetNormalizeAcrossScale(false);
  
  pMessenger = new GateFSThermalActorMessenger(this);
  GateDebugMessageDec("Actor",4,"GateFSThermalActor() -- end"<<G4endl);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateFSThermalActor::SetBloodPerfusionByMaterial(G4bool b)
{
  mIsPerfusionActivated = b;
  mIsPerfusionByMaterial = b;
  mIsPerfusionByConstant = false;
  mIsPerfusionByImage = false;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateFSThermalActor::SetBloodPerfusionByConstant(G4double value)
{
  mIsPerfusionActivated = true;
  mIsPerfusionByMaterial = false;
  mIsPerfusionByConstant = true;
  mIsPerfusionByImage = false;
  mUserBloodPerfusionRate = value;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateFSThermalActor::SetBloodPerfusionByImage(G4String string)
{
  mIsPerfusionActivated = true;
  mIsPerfusionByMaterial = false;
  mIsPerfusionByConstant = false;
  mIsPerfusionByImage = true;
  mUserPerfusionImageName = string;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateFSThermalActor::SetMeasurementFilename(G4String string)
{
  mIsMeasurementActivated = true;
  mMeasurementFilename = string;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
/// Destructor
GateFSThermalActor::~GateFSThermalActor()  {
  delete pMessenger;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
/// Constructor
void GateFSThermalActor::Construct() {

  GateDebugMessageInc("Actor", 4, "GateFSThermalActor -- Construct - begin" << G4endl);

  GateVImageActor::Construct();

  // Record the stepHitType
  mUserStepHitType = mStepHitType;

  // Enable callbacks
  EnableBeginOfRunAction(true);
  EnableEndOfRunAction(true); // for save
  EnableBeginOfEventAction(true);
  EnableEndOfEventAction(true);
  EnablePreUserTrackingAction(true);
  EnableUserSteppingAction(true);

  // Output Filenames
  mAbsorptionFilename = G4String(removeExtension(mSaveFilename))+"-AbsorptionMap."+G4String(getExtension(mSaveFilename));
  mHeatDiffusionFilename = G4String(removeExtension(mSaveFilename))+"-HeatDiffusionMap."+G4String(getExtension(mSaveFilename));

  // Set origin, transform, flag
  SetOriginTransformAndFlagToImage(mAbsorptionImage);

  // Resize and allocate images
  mAbsorptionImage.SetResolutionAndHalfSize(mResolution, mHalfSize, mPosition);
  mAbsorptionImage.Allocate();
  mAbsorptionImage.SetFilename(mAbsorptionFilename);

  // initialize ITK energy map from actor energy map
  GateImageDouble *test = dynamic_cast<GateImageDouble *>(&(mAbsorptionImage.GetValueImage()));
  mITKdoubleDuplicatorFilter->SetInputImage(ConvertEnergyImageToITKImage(test));
  mITKdoubleDuplicatorFilter->Update();
  mITKenergyMap = mITKdoubleDuplicatorFilter->GetOutput();
  mITKenergyMap->DisconnectPipeline();
  
  if(mIsMeasurementActivated) { ReadMeasurementFile(mITKenergyMap); }
  
  // construct diffusion masks
  // FS WARNING - needs a voxelised volume as attached volume
  GateVImageVolume *gateVoxelisedMap = dynamic_cast<GateVImageVolume *>(GetVolume());
  if(!gateVoxelisedMap)
  {
    GateError("Error: in its actual version, only voxelised volume can be used as attached volume.");
  }
  else
  {
    ConstructRegionMasks(gateVoxelisedMap);
  }

  mTimeStart = GateApplicationMgr::GetInstance()->GetTimeStart();
  mCurrentTime = mTimeStart;

  // Print information
  GateMessage("Actor", 1,
              "\tThermalActor    = '" << GetObjectName() << "'" << G4endl <<
              "\tAbsorptionFilename      = " << mAbsorptionFilename << G4endl);

  ResetData();
  GateMessageDec("Actor", 4, "GateFSThermalActor -- Construct - end" << G4endl);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
/// Save data
void GateFSThermalActor::SaveData()
{
  mAbsorptionImage.SaveData(mCurrentEvent+1);
}
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
void GateFSThermalActor::ResetData()
{
  mAbsorptionImage.Reset();
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateFSThermalActor::BeginOfRunAction(const G4Run * r)
{  
  GateVActor::BeginOfRunAction(r);

  DD("BeginOfRunAction::Begin");
  
  mCurrentTime = GateApplicationMgr::GetInstance()->GetCurrentTime();
  mTimeStop = GateApplicationMgr::GetInstance()->GetTimeStop();  

  GateDebugMessage("Actor", 3, "GateFSThermalActor -- Begin of Run" << G4endl);

  DD("BeginOfRunAction::End");
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateFSThermalActor::BeginOfEventAction(const G4Event * e) {
  GateVActor::BeginOfEventAction(e);

  mCurrentEvent++;
  GateDebugMessage("Actor", 3, "GateFSThermalActor -- Begin of Event: "<<mCurrentEvent << G4endl);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateFSThermalActor::EndOfEventAction(const G4Event *)
{
  double currentTime = GateApplicationMgr::GetInstance()->GetCurrentTime();  
  double tmpTime = currentTime-mCurrentTime;
  
  if(mIsDiffusionActivated) { ApplyStepDiffusion(tmpTime, false); }
  if(mIsPerfusionActivated) { ApplyStepPerfusion(tmpTime, false); }
  if(mIsMeasurementActivated) { ApplyStepMeasurement(tmpTime, false); }
  
  mCurrentTime = currentTime;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateFSThermalActor::EndOfRunAction(const G4Run* r)
{
  GateVActor::EndOfRunAction(r);

  DD("EndOfRunAction::Begin");
  
  double currentTime = GateApplicationMgr::GetInstance()->GetCurrentTime();  
  double tmpTime = currentTime-mCurrentTime;

  if(mIsDiffusionActivated) { ApplyStepDiffusion(tmpTime, true); }
  if(mIsPerfusionActivated) { ApplyStepPerfusion(tmpTime, true); }
  if(mIsMeasurementActivated) { ApplyStepMeasurement(tmpTime, true); }
  
  mCurrentTime = currentTime;
  
//   G4cout << "start: " << mTimeStart / s << " | current: " << mCurrentTime / s << " | stop:" << mTimeStop / s << G4endl;
  
  CropFilterType::Pointer crop = CropFilterType::New();
  DoubleImageType::SizeType cropSize;
  cropSize[0] = mCropSize;
  cropSize[1] = mCropSize;
  cropSize[2] = mCropSize;
  crop->SetBoundaryCropSize(cropSize);
  crop->SetInput(mITKenergyMap);
  crop->Update();
//   SaveITKimage(crop->GetOutput(), mAbsorptionFilename);
  SaveITKimage(mITKenergyMap, mAbsorptionFilename);
  
  if(mUserRelaxationTime > 0.0) { ApplyUserRelaxation(); }

  DD("EndOfRunAction::End");
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateFSThermalActor::UserPreTrackActionInVoxel(const int /*index*/, const G4Track* track) {

  if(track->GetDefinition()->GetParticleName() == "opticalphoton") { mStepHitType = PostStepHitType; }
  else { mStepHitType = mUserStepHitType; }
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateFSThermalActor::UserSteppingActionInVoxel(const int index, const G4Step* step) {

  GateDebugMessageInc("Actor", 4, "GateFSThermalActor -- UserSteppingActionInVoxel - begin" << G4endl);
  
  // FS: edep according to energyDeposit (not kinetic energy) -> to check
  const double edep = step->GetPostStepPoint()->GetKineticEnergy()/eV;  // in eV
//   const double edep = step->GetTotalEnergyDeposit() / eV;
  
  const G4String process = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  
  // if no energy is deposited or energy is deposited outside image => do nothing
  if (step->GetPostStepPoint()->GetKineticEnergy() == 0) {
    GateDebugMessage("Actor", 5, "edep == 0 : do nothing" << G4endl);
    GateDebugMessageDec("Actor", 4, "GateFSThermalActor -- UserSteppingActionInVoxel -- end" << G4endl);
    return;
  }

  if (index <0) {
    GateDebugMessage("Actor", 5, "index<0 : do nothing" << G4endl);
    GateDebugMessageDec("Actor", 4, "GateFSThermalActor -- UserSteppingActionInVoxel -- end" << G4endl);
    return;
  }

  GateDebugMessage("Actor", 2, "GateFSThermalActor -- UserSteppingActionInVoxel:\tedep = " << G4BestUnit(edep, "Energy") << G4endl);

  if ( process == "NanoAbsorption" || process == "OpticalAbsorption" )
  {
    // add energy in the gate image
    mAbsorptionImage.AddValue(index, edep);

    // add energy in the ITK image (for diffusion and perfusion)
    G4ThreeVector gatePixelCoordinate = mImage.GetCoordinatesFromIndex(index);
    mITKdoubleIndex[0] = gatePixelCoordinate.getX();
    mITKdoubleIndex[1] = gatePixelCoordinate.getY();
    mITKdoubleIndex[2] = gatePixelCoordinate.getZ();
    double pixelValue = mITKenergyMap->GetPixel(mITKdoubleIndex);
    mITKenergyMap->SetPixel(mITKdoubleIndex, pixelValue + edep);
  }
  
  GateDebugMessageDec("Actor", 4, "GateFSThermalActor -- UserSteppingActionInVoxel -- end" << G4endl);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateFSThermalActor::ApplyStepDiffusion(double timeStep, bool forced)
{
  std::map<G4Material *, DiffusionStruct>::iterator itMap;
  std::map<G4Material *, DiffusionStruct>::iterator itBegin = mMaterialToDiffusionStruct.begin();
  std::map<G4Material *, DiffusionStruct>::iterator itEnd = mMaterialToDiffusionStruct.end();
  
  mITKmultiplyImageFilter->SetInput1(mITKenergyMap);
  
  for(itMap = itBegin; itMap != itEnd; ++itMap)
  {
    bool checkDiffusion = itMap->second.CheckDiffusionTime(timeStep, forced);    
    if(checkDiffusion)
    {
      G4cout << "mat "<< itMap->first->GetName() << " is diffusing | sigma = " << itMap->second.sigma / mm << " mm | currentTimeStep = " << itMap->second.currentTimeStep / s << " s | timeStep = " << timeStep / s << " s | currentTime = " << mCurrentTime / s << " s" << G4endl;
      itMap->second.diffusionNumber++;

      // 1. multiply mask with energy map
      mITKmultiplyImageFilter->SetInput2(itMap->second.mask);
      mITKmultiplyImageFilter->Update();
      
      mITKsubtractImageFilter->SetInput1(mITKenergyMap);
      mITKsubtractImageFilter->SetInput2(mITKmultiplyImageFilter->GetOutput());
      mITKsubtractImageFilter->Update();
      
      // 2. apply recursive gaussian filter with corresponding diffusivity
      double sigma = itMap->second.sigma;
      mITKgaussianFilterX->SetInput(mITKmultiplyImageFilter->GetOutput());
      mITKgaussianFilterX->SetSigma(sigma);
      mITKgaussianFilterX->Update();
      mITKgaussianFilterY->SetInput(mITKgaussianFilterX->GetOutput());
      mITKgaussianFilterY->SetSigma(sigma);
      mITKgaussianFilterY->Update();
      mITKgaussianFilterZ->SetInput(mITKgaussianFilterY->GetOutput());
      mITKgaussianFilterZ->SetSigma(sigma);
      mITKgaussianFilterZ->Update();
      
//       DoubleImageType::Pointer tmpImg = mITKmultiplyImageFilter->GetOutput();
      DoubleImageType::Pointer tmpImg = mITKgaussianFilterZ->GetOutput();
      
//       if(mIsPerfusionActivated) {
//         
//         MultiplyFilterType::Pointer multConst = MultiplyFilterType::New();
//         MultiplyFilterType::Pointer multImage = MultiplyFilterType::New();
//         ExpFilterType::Pointer expImage = ExpFilterType::New();
//         
//         multConst->SetInput(mITKperfusionRateMap);
//         multConst->SetConstant(-itMap->second.currentTimeStep);
//         expImage->SetInput(multConst->GetOutput());
//         multImage->SetInput1(tmpImg);
//         multImage->SetInput2(expImage->GetOutput());
//         multImage->Update();
// 
//         mITKdoubleIndex[0] = 117;
//         mITKdoubleIndex[1] = 117;
//         mITKdoubleIndex[2] = 25;
//         DD(mITKperfusionRateMap->GetPixel(mITKdoubleIndex));
//         DD(itMap->second.currentTimeStep);
//         
//         tmpImg = multImage->GetOutput();
//         tmpImg->DisconnectPipeline();
//         
//         
// //         DoubleIteratorType it(mITKgaussianFilterZ->GetOutput(), mITKgaussianFilterZ->GetOutput()->GetRequestedRegion());
// // //         for(it.GoToBegin(); !it.IsAtEnd(); ++it) { it.Set(it.Get() * std::exp(- mITKperfusionRateMap->GetPixel(it.GetIndex()) * timeStep)); }
// //         for(it.GoToBegin(); !it.IsAtEnd(); ++it) { it.Set(it.Get() * std::exp(- 0.142 /s * timeStep)); }
//         
//       }
      
      mITKaddImageFilter->SetInput1(mITKsubtractImageFilter->GetOutput());
      mITKaddImageFilter->SetInput2(tmpImg);
      mITKaddImageFilter->Update();
      
      mITKenergyMap = mITKaddImageFilter->GetOutput();
      mITKenergyMap->DisconnectPipeline();
    }
  }
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

void GateFSThermalActor::ApplyUserRelaxation()
{
  // -------------------------------------------------------------------------------
  // mimic time
  
  double totalTime = mCurrentTime + mUserRelaxationTime;

  DD("BeginUserDiffusion");
  if(mIsDiffusionActivated)
  {
    DD("diffusion YES ; perfusion MAYBE");
    while(mCurrentTime < totalTime)
    {
      ApplyStepDiffusion(mMinTimeStep, false);
      if(mIsPerfusionActivated) { ApplyStepPerfusion(mMinTimeStep, false); }
      if(mIsMeasurementActivated) { ApplyStepMeasurement(mMinTimeStep, false); }
      
//       std::cout << std::setprecision(10) << "start: " << mTimeStart / s << " | current: " << mCurrentTime / s << " | timeStep: " << mMinTimeStep / s << " | stop:" << mTimeStop / s <<" | userStop:" << totalTime / s << std::endl;
      mCurrentTime += mMinTimeStep;
    }

    std::map<G4Material *, DiffusionStruct>::iterator itMap;
    std::map<G4Material *, DiffusionStruct>::iterator itBegin = mMaterialToDiffusionStruct.begin();
    std::map<G4Material *, DiffusionStruct>::iterator itEnd = mMaterialToDiffusionStruct.end();
    for(itMap = itBegin; itMap != itEnd; ++itMap)
    {
      G4cout << "mat " << itMap->first->GetName() <<" has diffused " << itMap->second.diffusionNumber << " times | totalTime = " << itMap->second.totalTime / s << " s" << G4endl;
    }    
  }
  else if(mIsPerfusionActivated) {
    DD("diffusion NO ; perfusion YES");
    ApplyStepPerfusion(mUserRelaxationTime, true);
  }


  CropFilterType::Pointer crop = CropFilterType::New();
  DoubleImageType::SizeType cropSize;
  cropSize[0] = mCropSize;
  cropSize[1] = mCropSize;
  cropSize[2] = mCropSize;
  crop->SetBoundaryCropSize(cropSize);
  crop->SetInput(mITKenergyMap);
  crop->Update();
//   SaveITKimage(crop->GetOutput(), mHeatDiffusionFilename);
  SaveITKimage(mITKenergyMap, mHeatDiffusionFilename);
  
  for(unsigned int i=0; i<mMeasurementPoints.size(); i++)
  {
    std::vector<double> timeList = mMeasurementPoints[i].timeList;
    std::vector<double> measList = mMeasurementPoints[i].measList;

    G4String filename = G4String(removeExtension(mSaveFilename))+"-ROI-"+ std::to_string(mMeasurementPoints[i].label) +".txt";
    std::ofstream os(filename);
    
    os << "time(s) \tenergy(eV) \tNvoxel="<< mMeasurementPoints[i].indexList.size() << std::endl;
    for(unsigned int j=0; j<timeList.size(); j++)
    {
      os.precision(10);
      os << timeList[j] / s << "\t" << measList[j] << "\t" << std::endl;
    }
    os.close();
  }
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateFSThermalActor::ApplyStepPerfusion(double timeStep, bool forced)
{
  mPerfusionTimer += timeStep;
  if(mPerfusionTimer >= mMinPerfusionTimeStep or forced)
  {
    G4cout << "blood perfusion | currentTimeStep = " << mPerfusionTimer / s << " s | timeStep = " << timeStep / s << " s | currentTime = " << mCurrentTime / s << " s" << G4endl;
    
    MultiplyFilterType::Pointer multConst = MultiplyFilterType::New();
    MultiplyFilterType::Pointer multImage = MultiplyFilterType::New();
    ExpFilterType::Pointer expImage = ExpFilterType::New();
    
    multConst->SetInput(mITKperfusionRateMap);
    multConst->SetConstant(-mPerfusionTimer);
    expImage->SetInput(multConst->GetOutput());
    multImage->SetInput1(mITKenergyMap);
    multImage->SetInput2(expImage->GetOutput());
    multImage->Update();

    mITKdoubleIndex[0] = 117;
    mITKdoubleIndex[1] = 117;
    mITKdoubleIndex[2] = 25;
//     DD(mITKperfusionRateMap->GetPixel(mITKdoubleIndex));
//     DD(mPerfusionTimer / s);
    
    mITKenergyMap = multImage->GetOutput();
    mITKenergyMap->DisconnectPipeline();
    
    mPerfusionTimer = 0.0;
  }
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateFSThermalActor::ApplyStepMeasurement(double timeStep, bool forced)
{
  for(unsigned int i=0; i<mMeasurementPoints.size(); i++)
  {
    bool test = mMeasurementPoints[i].CheckMeasurementTime(timeStep, forced);
    if(test)
    {
      DD(mMeasurementPoints[i].totalTime / s);
      double value = 0.0;
      std::vector<DoubleImageType::IndexType> indices = mMeasurementPoints[i].indexList;
      for(unsigned j=0; j<indices.size(); j++) { value += mITKenergyMap->GetPixel(indices[j]); }
      mMeasurementPoints[i].SetValue(mMeasurementPoints[i].totalTime, value);
    }
  }
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
double GateFSThermalActor::GetPropertyFromMaterial(const G4Material *mat, G4String prop, G4double unit)
{
  G4MaterialPropertiesTable *materialPropertyTable = mat->GetMaterialPropertiesTable();
  if(materialPropertyTable)
  {
    return materialPropertyTable->GetConstProperty(prop) * unit;
  }
  else
  {
    return 0.0;
  }
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateFSThermalActor::ConstructRegionMasks(GateVImageVolume *gateImage)
{
  // This function creates a 'diffusion struct' (see .hh) for each material in the
  // considered voxelised volume. Each struct is composed of:
  // 1. diffusivity/period (constant values)
  // 2. image mask (1 if corresponding material, 0 otherwise)
  // 3. timer/sigma (variable values)
  
  // Create an ITK copy of the GATE label image
  GateImage *gateImg = gateImage->GetImage();
  FloatImageType::Pointer originLabelImage = ConvertGateImageToITKImage(gateImg);
  mITKfloatDuplicatorFilter->SetInputImage(originLabelImage);
  mITKfloatDuplicatorFilter->Update();
  FloatImageType::Pointer labelImage = mITKfloatDuplicatorFilter->GetOutput();

  
  // find the maximum resolution in image that will be use to limit the diffusion to the voxel size
  double resolutionMax = 0.0;
  for (unsigned int i = 0; i < 3; i++)
  {
    double spacing = gateImg->GetVoxelSize()[i];
    if (spacing > resolutionMax) { resolutionMax = spacing; }
  }
  resolutionMax = resolutionMax * mm;

//   typedef itk::CropImageFilter<FloatImageType, FloatImageType> tmpCropFilterType;
//   tmpCropFilterType::Pointer crop = tmpCropFilterType::New();
//   FloatImageType::SizeType cropSize;
//   cropSize[0] = 0;
//   cropSize[1] = 0;
//   cropSize[2] = 150;
//   crop->SetBoundaryCropSize(cropSize);
//   crop->SetInput(labelImage);
//   crop->Update();
//   SaveITKimage(crop->GetOutput(), "myCropMouse.mhd");

  // Create a temporary materialToLabel map in order to regroup every voxel with the same material under the same label
  // and, by consequence, minimize the number of image masks when applying the diffusion process.
  std::map<float, G4Material *> labelToNewMaterial;
  std::map<G4Material *, float> materialToNewLabel;
  std::map<G4Material *, float>::iterator itNewLabel;
  gateImage->BuildLabelToG4MaterialVector(mMaterialList);
  int newLabel = 0;
  
  FloatIteratorType itImage(labelImage, labelImage->GetRequestedRegion());
  for(itImage.GoToBegin(); !itImage.IsAtEnd(); ++itImage)
  {
    // 1. get voxel label and material
    int label = itImage.Get();
    G4Material *mat = mMaterialList[label];

    // 2. check if material already exists, create map entry if not
    itNewLabel = materialToNewLabel.find(mat);
    if(itNewLabel == materialToNewLabel.end())
    {
      materialToNewLabel.insert(std::make_pair(mat, (float)newLabel));      
      labelToNewMaterial.insert(std::make_pair((float)newLabel, mat));      
      newLabel++;
    }
    
    // 3. replace old label with the new one
    itNewLabel = materialToNewLabel.find(mat);
    itImage.Set(itNewLabel->second);
    labelImage->Update();
  }

  // Initialisation of the ITK binary filter
  BinaryThresholdFilterType::Pointer binaryThresholdFilter = BinaryThresholdFilterType::New();
  binaryThresholdFilter->SetInput(labelImage);
  binaryThresholdFilter->SetOutsideValue(float(0.0));
  binaryThresholdFilter->SetInsideValue(float(1.0));

  FloatPadFilterType::Pointer pad = FloatPadFilterType::New();
  FloatImageType::SizeType border;
  border.Fill(mCropSize);
  pad->SetPadBound(border);
  
  for(itNewLabel = materialToNewLabel.begin(); itNewLabel != materialToNewLabel.end(); ++itNewLabel)
  {
    newLabel = itNewLabel->second;
    binaryThresholdFilter->SetLowerThreshold(newLabel-0.1);
    binaryThresholdFilter->SetUpperThreshold(newLabel+0.1);
    binaryThresholdFilter->Update();
    
    if(newLabel == 1.0)
    {
      DD("Test succeded !");
      pad->SetConstant(0.0);
    }
    else
    {
      DD("Test failed ...");
      pad->SetConstant(0.0);
    }
    
    pad->SetInput(binaryThresholdFilter->GetOutput());
    pad->Update();
    
    double diffusivity = GetPropertyFromMaterial(itNewLabel->first, "DIFFUSIVITY", mm2/s);
//       FloatImageType::Pointer mask = pad->GetOutput();
    FloatImageType::Pointer mask = binaryThresholdFilter->GetOutput();
    mask->DisconnectPipeline();
    DiffusionStruct newDiffStruct(diffusivity, 1.0 * resolutionMax, mask);
    mMaterialToDiffusionStruct.insert(std::make_pair(itNewLabel->first, newDiffStruct));
  }

  double diffusivityMin = 1.0e10;
  double diffusivityMax = 0.0;
  
  std::map<G4Material *, DiffusionStruct>::iterator itMap;  
  for(itMap = mMaterialToDiffusionStruct.begin(); itMap != mMaterialToDiffusionStruct.end(); ++itMap)
  {
    if(itMap->second.diffusivity > diffusivityMax) { diffusivityMax = itMap->second.diffusivity; }
    if(itMap->second.diffusivity < diffusivityMin and itMap->second.diffusivity > 0) { diffusivityMin = itMap->second.diffusivity; }
  }

  double durationMin = 1.0 * resolutionMax * resolutionMax / (2 * diffusivityMax);
  double durationMax = 1.0 * resolutionMax * resolutionMax / (2 * diffusivityMin);
  G4cout << "diffMax = " << diffusivityMax / (mm2/s) << " durationMin = " << durationMin / s  << " | diffMin = " << diffusivityMin / (mm2/s) << " durationMax = " << durationMax / s << G4endl;
  
  mMinTimeStep = 1.0 * resolutionMax * resolutionMax / (2.0 * diffusivityMax);
  
  // debug
//   std::map<G4Material *, DiffusionStruct>::iterator itTmp;
//   for(itTmp=mMaterialToDiffusionStruct.begin(); itTmp!=mMaterialToDiffusionStruct.end(); ++itTmp)
//   {
//     G4cout << "mat "<< itTmp->first->GetName() << " | diff = " << itTmp->second.diffusivity / (mm2/s)
//            << " mm2.s-1 | period = " << itTmp->second.period / s
//            << " s | timer = " << itTmp->second.timer / s << " s " << G4endl;
//            
//     std::ostringstream temp;
//     temp << itTmp->first->GetName();
//     SaveITKimage(itTmp->second.mask, "output/itkMask_" + temp.str() + ".mhd");
//   }

  double mMinPerfusionCoef = 1.0e10;
  double mMaxPerfusionCoef = 0.0;
  
  // Construct perfusionMap
  if(mIsPerfusionActivated)
  {
    if(mUserBloodDensity<=0.0) { GateError("Error: Please, set the 'bloodDensity' in order to use the blood perfusion process."); }
    if(mUserBloodHeatCapacity<=0.0) { GateError("Error: Please, set the 'bloodHeatCapacity' in order to use the blood perfusion process."); }
//     if(mUserTissueHeatCapacity<=0.0) { GateError("Error: Please, set the 'tissueHeatCapacity' in order to use the blood perfusion process."); }
    DoubleImageType *tmpPerfusionRateMap = DoubleImageType::New();
    
    if(mIsPerfusionByImage) {
      DD("perfByImage");
      mITKdoubleReaderFilter->SetFileName(mUserPerfusionImageName);
      mITKdoubleReaderFilter->Update();
      DD("test1");
      tmpPerfusionRateMap = mITKdoubleReaderFilter->GetOutput();
      DD("test2");
    }

    mITKdoubleDuplicatorFilter->SetInputImage(mITKenergyMap);
    mITKdoubleDuplicatorFilter->Update();
    mITKperfusionRateMap = mITKdoubleDuplicatorFilter->GetOutput();
    
    for(itImage.GoToBegin(); !itImage.IsAtEnd(); ++itImage)
    {
      DoubleImageType::IndexType index = itImage.GetIndex();
      
      // 1. get voxel label and material
      std::map<float, G4Material *>::iterator itLabel = labelToNewMaterial.find(itImage.Get());
      G4Material *mat = itLabel->second;
      
      double perfusionRate = mUserBloodPerfusionRate;
      if(mIsPerfusionByMaterial) { perfusionRate = GetPropertyFromMaterial(mat, "PERFUSIONRATE", 1./s); }
      else if(mIsPerfusionByImage) { perfusionRate = tmpPerfusionRateMap->GetPixel(index) / s; }
      
      double tissueHeatCapacity;
      if(mUserTissueHeatCapacity > 0.0) { tissueHeatCapacity = mUserTissueHeatCapacity; }
      else { tissueHeatCapacity = GetPropertyFromMaterial(mat, "HEATCAPACITY", joule/(kg*kelvin)); }
      
      double tissueDensity = mat->GetDensity();
      
      double perfusionCoef = perfusionRate * (mUserBloodDensity * mUserBloodHeatCapacity) / (tissueDensity * tissueHeatCapacity);
      
//       if(index[2]==15) G4cout <<index[0]<<" "<<index[1]<<" "<<index[2]<<" | perf: "<< perfusionRate * s <<" s-1 | tHC: "<<tissueHeatCapacity / (joule/(kg*kelvin)) <<" J.kg-1.K-1 | td: "<< tissueDensity / (g/cm3) <<" g/cm-3 | bHC: "<< mUserBloodHeatCapacity / (joule/(kg*kelvin)) <<" J.kg-1.K-1 | bd: "<< mUserBloodDensity /(g/cm3)<<" g.cm-3 |"<< G4endl;
      
      if(perfusionCoef < mMinPerfusionCoef) { mMinPerfusionCoef = perfusionCoef; }
      if(perfusionCoef > mMaxPerfusionCoef) { mMaxPerfusionCoef = perfusionCoef; }
      
      mITKperfusionRateMap->SetPixel(index, perfusionCoef);
    }
//     SaveITKimage(mITKperfusionRateMap, "output/perfusionRateMap.mhd");

    mMinPerfusionTimeStep = -log(mPerfusionRatio) / mMaxPerfusionCoef;

    DD(mMinPerfusionCoef * s);
    DD(mMaxPerfusionCoef * s);
    DD(mMinPerfusionTimeStep / s);
  }
  
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
GateFSThermalActor::FloatImageType::Pointer GateFSThermalActor::ConvertGateImageToITKImage(GateImage *gateImg)
{
  mGateToITKImageFilter = ImportFilterType::New();
  ImportFilterType::SizeType  size;
  double origin[3];
  double spacing[3];
  for (unsigned int i = 0; i < 3; i++)
  {
    size[i] = gateImg->GetResolution()[i];
    spacing[i] = gateImg->GetVoxelSize()[i];
    origin[i] = -gateImg->GetHalfSize()[i] + 0.5 * spacing[i];
  }

  ImportFilterType::IndexType start;
  start.Fill(0);
  ImportFilterType::RegionType region;
  region.SetIndex( start );
  region.SetSize(  size  );
 
  mGateToITKImageFilter->SetRegion(region);
  mGateToITKImageFilter->SetOrigin(origin);
  mGateToITKImageFilter->SetSpacing(spacing);
 
  const unsigned int numberOfPixels =  size[0] * size[1] * size[2];
  const bool importImageFilterWillOwnTheBuffer = false;
  mGateToITKImageFilter->SetImportPointer(&*(gateImg->begin()), numberOfPixels, importImageFilterWillOwnTheBuffer);
  mGateToITKImageFilter->Update();

//   SaveITKimage(mGateToITKImageFilter->GetOutput(), "output/GateToITKimageCT.mhd");

  FloatImageType::Pointer output = mGateToITKImageFilter->GetOutput();
  return output;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
GateFSThermalActor::DoubleImageType::Pointer GateFSThermalActor::ConvertEnergyImageToITKImage(GateImageDouble *gateImg)
{
  typedef itk::ImportImageFilter<double, 3> DoubleImportFilterType;
  DoubleImportFilterType::Pointer importFilter = DoubleImportFilterType::New();

  DoubleImportFilterType::SizeType  size;
  double origin[3];
  double spacing[3];
  for (unsigned int i = 0; i < 3; i++)
  {
    size[i] = gateImg->GetResolution()[i];
    spacing[i] = gateImg->GetVoxelSize()[i];
    origin[i] = -gateImg->GetHalfSize()[i] + 0.5 * spacing[i];
  }

  DoubleImportFilterType::IndexType start;
  start.Fill(0);
  DoubleImportFilterType::RegionType region;
  region.SetIndex( start );
  region.SetSize(  size  );

  importFilter->SetRegion(region);
  importFilter->SetOrigin(origin);
  importFilter->SetSpacing(spacing);

  const unsigned int numberOfPixels =  size[0] * size[1] * size[2];
  const bool importImageFilterWillOwnTheBuffer = false;
  importFilter->SetImportPointer(&*(gateImg->begin()), numberOfPixels, importImageFilterWillOwnTheBuffer);
  importFilter->Update();
  
//   SaveITKimage(importFilter->GetOutput(), "output/GateToITKimageEnergy.mhd");

//   DoubleImageType::Pointer output = importFilter->GetOutput();
//   return output;
  
  DoublePadFilterType::Pointer pad = DoublePadFilterType::New();
  DoubleImageType::SizeType border;
  border.Fill(mCropSize);
  pad->SetPadBound(border);
  pad->SetConstant(0.0);
  pad->SetInput(importFilter->GetOutput());
  pad->Update();

//   DoubleImageType::Pointer output = pad->GetOutput();
  DoubleImageType::Pointer output = importFilter->GetOutput();
  
  SaveITKimage(output, "output/GateToITKimageEnergy.mhd");
  
  return output;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateFSThermalActor::ReadMeasurementFile(DoubleImageType::Pointer img)
{
  // Open file
  std::ifstream is;
  OpenFileInput(mMeasurementFilename, is);
  skipComment(is);

  // Use R et/or T ? Time s  | Translation mm | Rotation deg
  double timeUnit=0;
  if (!ReadColNameAndUnit(is, "Time", timeUnit)) {
    GateError("The file '" << mMeasurementFilename << "' need to begin with 'Time'\n");
  }

  // Loop line
  skipComment(is);
  while (is) {
    int label = lrint(ReadDouble(is));
    double timeStep = ReadDouble(is) * timeUnit;
    signed int lx = lrint(ReadDouble(is));
    unsigned int ux = lrint(ReadDouble(is));
    signed int ly = lrint(ReadDouble(is));
    unsigned int uy = lrint(ReadDouble(is));
    signed int lz = lrint(ReadDouble(is));
    unsigned int uz = lrint(ReadDouble(is));
    
    G4cout<<"lx,ux: [" <<lx<<","<<ux<<"] ly,uy: [" <<ly<<","<<uy<<"] lz,uz: [" <<lz<<","<<uz<<"]" << G4endl;
    
    DoubleImageType::SizeType size = img->GetRequestedRegion().GetSize();    
    if(timeStep>0 and lx>-1 and ly>-1 and lz>-1 and ux<size[0] and uy<size[1] and uz<size[2])
    {
      MeasurementStruct newMeasPoint(label, timeStep);
      for(unsigned int i=lx; i<ux+1; i++) {
        for(unsigned int j=ly; j<uy+1; j++) {
          for(unsigned int k=lz; k<uz+1; k++) {
            DoubleImageType::IndexType newIndex;
            newIndex[0] = i;
            newIndex[1] = j;
            newIndex[2] = k;
            newMeasPoint.indexList.push_back(newIndex);
          }
        }        
      }

      mMeasurementPoints.push_back(newMeasPoint);
      
//       G4cout<<"period: "<< newMeasPoint.period <<" | lx: "<<lx<<" ux: "<<ux<<G4endl;
//       DD(newMeasPoint.indexList.size());
//       for(unsigned int i=0; i<newMeasPoint.indexList.size(); i++)
//       {
//         G4cout<<"x: "<< newMeasPoint.indexList[i][0]<<"x: "<< newMeasPoint.indexList[i][1] <<"x: "<< newMeasPoint.indexList[i][2] <<G4endl;
//       }      
    }
  
    skipComment(is);
  }
  
  // End
  is.close();
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateFSThermalActor::SaveITKimage(FloatImageType::Pointer img, G4String name)
{
  itk::ImageFileWriter<FloatImageType>::Pointer writer = itk::ImageFileWriter<FloatImageType>::New();
  writer->SetFileName(name);
  writer->SetInput(img);
  writer->Update();
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateFSThermalActor::SaveITKimage(DoubleImageType::Pointer img, G4String name)
{
  itk::ImageFileWriter<DoubleImageType>::Pointer writer = itk::ImageFileWriter<DoubleImageType>::New();
  writer->SetFileName(name);
  writer->SetInput(img);
  writer->Update();
}
//-----------------------------------------------------------------------------

#endif // end define USE_ITK
