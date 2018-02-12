/*----------------------
  Copyright (C): OpenGATE Collaboration

  This software is distributed under the terms
  of the GNU Lesser General  Public Licence (LGPL)
  See GATE/LICENSE.txt for further details
  ----------------------*/

/*
  \class GateThermalActor
  \author vesna.cuplov@gmail.com
  \author fsmekens@gmail.com
  \brief Class GateThermalActor : This actor produces voxelised images of the heat diffusion in tissue.

*/

#ifndef GATEFSTHERMALACTOR_HH
#define GATEFSTHERMALACTOR_HH

#include <G4NistManager.hh>
#include "GateVImageActor.hh"
#include "GateActorManager.hh"
#include "G4UnitsTable.hh"
#include "GateFSThermalActorMessenger.hh"
#include "GateImageWithStatistic.hh"
#include "GateVImageVolume.hh"

#include "G4Event.hh"
#include <time.h>

// itk
#include "GateConfiguration.h"
#ifdef GATE_USE_ITK
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageIterator.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkExpImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImportImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkCropImageFilter.h"

#endif

class GateFSThermalActor : public GateVImageActor
{
public:
  //-----------------------------------------------------------------------------
  // Actor name
  virtual ~GateFSThermalActor();

  FCT_FOR_AUTO_CREATOR_ACTOR(GateFSThermalActor)

  // typedef for itk
  typedef itk::Image<float, 3> FloatImageType;
  typedef itk::Image<double, 3> DoubleImageType;
  typedef itk::ImageRegionIterator<FloatImageType> FloatIteratorType;
  typedef itk::ImageRegionIterator<DoubleImageType> DoubleIteratorType;
  typedef itk::ImageDuplicator<FloatImageType> FloatDuplicatorType;
  typedef itk::ImageDuplicator<DoubleImageType> DoubleDuplicatorType;
  typedef itk::ImageFileReader<FloatImageType> FloatReaderType;
  typedef itk::ImageFileReader<DoubleImageType> DoubleReaderType;
  typedef itk::ConstantPadImageFilter<FloatImageType, FloatImageType> FloatPadFilterType;
  typedef itk::ConstantPadImageFilter<DoubleImageType, DoubleImageType> DoublePadFilterType;
  typedef itk::CropImageFilter<DoubleImageType, DoubleImageType> CropFilterType;
  typedef itk::MultiplyImageFilter<DoubleImageType, FloatImageType, DoubleImageType> MultiplyFilterType;
  typedef itk::ExpImageFilter<DoubleImageType, FloatImageType> ExpFilterType;
  typedef itk::AddImageFilter<DoubleImageType, DoubleImageType, DoubleImageType> AddImageFilterType;
  typedef itk::SubtractImageFilter<DoubleImageType, DoubleImageType> SubtractImageFilterType;
  typedef itk::BinaryThresholdImageFilter<FloatImageType, FloatImageType> BinaryThresholdFilterType;
  typedef itk::MaximumImageFilter<FloatImageType, FloatImageType, FloatImageType> MaxImageFilterType;
  typedef itk::RecursiveGaussianImageFilter<DoubleImageType, DoubleImageType> GaussianFilterType;
  typedef itk::DiscreteGaussianImageFilter<DoubleImageType, DoubleImageType> DGaussianFilterType;
  typedef itk::ImportImageFilter<float, 3> ImportFilterType;
  
  //-----------------------------------------------------------------------------
  struct DiffusionStruct {
    double diffusivity;
    double period;
    double sigma;
    FloatImageType::Pointer mask;
    double timer;
    double currentTimeStep;
    double totalTime;
    int diffusionNumber;
    
    DiffusionStruct(double c, double s, FloatImageType::Pointer m) {
      diffusivity = c;
      sigma = s;
      mask = m;
      period = sigma * sigma / (2.0 * diffusivity);
      timer = 0.0;
      currentTimeStep = 0.0;
      totalTime = 0.0;
      diffusionNumber = 0;
    }
    
    bool CheckDiffusionTime(double stepTime, bool forced) {
      timer += stepTime;
      totalTime += stepTime;
//       std::cout << std::setprecision(10) << "timer: " << timer / s << " | period: " << period / s << std::endl;
      sigma = sqrt(2.0 * timer * diffusivity);
      if((timer >= period or forced) and sigma > 0.0) {
        currentTimeStep = timer;
        timer = 0.0;
        return true;
      }
      else { return false; }
    }
  };
  //-----------------------------------------------------------------------------
  struct MeasurementStruct {
    int label;
    double period;
    double timer;
    double totalTime;
    std::vector<DoubleImageType::IndexType> indexList;
    std::vector<double> timeList;
    std::vector<double> measList;
    
    MeasurementStruct(int l, double t) {
      label = l;
      period = t;
      timer = 0.0;
      totalTime = 0.0;
      indexList.clear();
      timeList.clear();
      measList.clear();
    }
    
    bool CheckMeasurementTime(double stepTime, bool forced) {
      timer += stepTime;
      totalTime += stepTime;
      if(timer >= period or forced)
      {
        timer = 0.0;
        return true;
      }
      else { return false; }
    }
    
    void SetValue(double t, double m) {
      timeList.push_back(t);
      measList.push_back(m);
    }
    
  };
  //-----------------------------------------------------------------------------
  // Constructs the sensor
  virtual void Construct();

  virtual void BeginOfRunAction(const G4Run *r);
  virtual void EndOfRunAction(const G4Run *); // default action (save)
  virtual void BeginOfEventAction(const G4Event *event);
  virtual void EndOfEventAction(const G4Event *event);
  virtual void UserSteppingActionInVoxel(const int index, const G4Step *step);
  virtual void UserPreTrackActionInVoxel(const int /*index*/, const G4Track *track);
  virtual void UserPostTrackActionInVoxel(const int /*index*/, const G4Track * /*t*/) {}

  //  Saves the data collected to the file
  virtual void SaveData();
  virtual void ResetData();

  // Scorer related
  virtual void Initialize(G4HCofThisEvent*){}
  virtual void EndOfEvent(G4HCofThisEvent*){}

  void ConstructRegionMasks(GateVImageVolume *);
  void ConstructPerfusionMap();
  double GetPropertyFromMaterial(const G4Material *, G4String, G4double);
  void ApplyStepPerfusion(double, bool);
  void ApplyStepDiffusion(double, bool);
  void ApplyStepMeasurement(double, bool);
  void ApplyUserRelaxation();
  
  FloatImageType::Pointer ConvertGateImageToITKImage(GateImage *);
  DoubleImageType::Pointer ConvertEnergyImageToITKImage(GateImageDouble *);
    
  void SaveITKimage(FloatImageType::Pointer, G4String);
  void SaveITKimage(DoubleImageType::Pointer, G4String);
  
  void SetBloodPerfusionByMaterial(G4bool);
  void SetBloodPerfusionByConstant(G4double);
  void SetBloodPerfusionByImage(G4String);
  void SetMeasurementFilename(G4String);
  void ReadMeasurementFile(DoubleImageType::Pointer);

  void setRelaxationTime(G4double t) { mUserRelaxationTime = t; }
  void setDiffusivity(G4double d) { mUserMaterialDiffusivity = d; }
  void setBloodDensity(G4double d) { mUserBloodDensity = d; }
  void setBloodHeatCapacity(G4double c) { mUserBloodHeatCapacity = c; }
  void setTissueDensity(G4double d) { mUserTissueDensity = d; }
  void setTissueHeatCapacity(G4double c) { mUserTissueHeatCapacity = c; }
  void setScale(G4double s) { mUserSimulationScale = s; }
  void enableStepDiffusion(G4bool b) { mIsDiffusionActivated = b; }
  
protected:

  G4double mTimeNow;

  GateFSThermalActor(G4String name, G4int depth=0);
  GateFSThermalActorMessenger * pMessenger;

  int mCurrentEvent;
  int counter;

  StepHitType mUserStepHitType;

  GateImageWithStatistic mAbsorptionImage;

  G4String mAbsorptionFilename;
  G4String mHeatDiffusionFilename;

  double mUserRelaxationTime;
  double mUserMaterialDiffusivity;
  double mUserSimulationScale;
  
  double mTimeStart;
  double mTimeStop;
  double mCurrentTime;
  double mPreviousTime;
  
  bool mIsDiffusionActivated;
  double mMinTimeStep;
  int mCropSize;
  std::vector<G4Material *> mMaterialList;
  std::map<int, DiffusionStruct> mLabelToDiffusionStruct;
  std::map<G4Material *, DiffusionStruct> mMaterialToDiffusionStruct;
  
  bool mIsMeasurementActivated;
  G4String mMeasurementFilename;
  std::vector<MeasurementStruct> mMeasurementPoints;
  
  bool mIsPerfusionActivated;
  bool mIsPerfusionByMaterial;
  bool mIsPerfusionByConstant;
  bool mIsPerfusionByImage;

  double deltaT;
  G4String mUserPerfusionImageName;
  double mUserBloodPerfusionRate;
  double mUserBloodDensity;
  double mUserBloodHeatCapacity;
  double mUserTissueDensity;
  double mUserTissueHeatCapacity;

  double mMinPerfusionCoef;
  double mMaxPerfusionCoef;
  double mPerfusionRatio;
  double mMinPerfusionTimeStep;
  double mPerfusionTimer;
  
  DoubleImageType::IndexType mITKdoubleIndex;
  DoubleImageType::Pointer mITKenergyMap;
  DoubleImageType::Pointer mITKperfusionRateMap;

  // ITK filters
  FloatReaderType::Pointer mITKfloatReaderFilter;
  DoubleReaderType::Pointer mITKdoubleReaderFilter;
  FloatDuplicatorType::Pointer mITKfloatDuplicatorFilter;
  DoubleDuplicatorType::Pointer mITKdoubleDuplicatorFilter;
  ImportFilterType::Pointer mGateToITKImageFilter;
  MultiplyFilterType::Pointer mITKmultiplyImageFilter;
  AddImageFilterType::Pointer mITKaddImageFilter;
  SubtractImageFilterType::Pointer mITKsubtractImageFilter;
  MaxImageFilterType::Pointer mITKmaxImageFilter;
  GaussianFilterType::Pointer mITKgaussianFilterX;
  GaussianFilterType::Pointer mITKgaussianFilterY;
  GaussianFilterType::Pointer mITKgaussianFilterZ;
};

MAKE_AUTO_CREATOR_ACTOR(FSThermalActor,GateFSThermalActor)

#endif /* end #define GATESIMULATIONSTATISTICACTOR_HH */
