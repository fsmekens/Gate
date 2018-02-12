/*----------------------
  Copyright (C): OpenGATE Collaboration

  This software is distributed under the terms
  of the GNU Lesser General  Public Licence (LGPL)
  See GATE/LICENSE.txt for further details
  ----------------------*/

/*
  \class GateThermalActorMessenger
  \brief This class is the GateThermalActor messenger. 
  \author vesna.cuplov@gmail.com
  \author fsmekens@gmail.com
*/

#ifndef GATEFSTHERMALACTORMESSENGER_HH
#define GATEFSTHERMALACTORMESSENGER_HH

#include "G4UIcmdWithABool.hh"
#include "GateImageActorMessenger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"

class GateFSThermalActor;
class GateFSThermalActorMessenger : public GateImageActorMessenger
{
public:
  GateFSThermalActorMessenger(GateFSThermalActor* sensor);
  virtual ~GateFSThermalActorMessenger();

  void BuildCommands(G4String base);
  void SetNewValue(G4UIcommand*, G4String);

protected:
  GateFSThermalActor * pThermalActor;
  G4UIcmdWithADoubleAndUnit* pRelaxationTimeCmd;
  G4UIcmdWithADouble* pDiffusivityCmd;
  
  G4UIcmdWithABool* pSetPerfusionRateByMaterialCmd;
  G4UIcmdWithADouble* pSetPerfusionRateByConstantCmd;
  G4UIcmdWithAString* pSetPerfusionRateByImageCmd;
  G4UIcmdWithADoubleAndUnit* pBloodDensityCmd;
  G4UIcmdWithADouble* pBloodHeatCapacityCmd;
  G4UIcmdWithADouble* pTissueHeatCapacityCmd;
  G4UIcmdWithADouble* pScaleCmd;
  G4UIcmdWithABool* pEnableWeightCmd;
  G4UIcmdWithABool* pEnableStepDiffusionCmd;
  G4UIcmdWithAString* pSetMeasurementFilenameCmd;
};

#endif /* end #define GATETHERMALACTORMESSENGER_HH*/

