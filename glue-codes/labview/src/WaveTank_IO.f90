!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2025  National Renewable Energy Laboratory
!
!    This file is a module specific to an experimental wave tank at NREL.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!**********************************************************************************************************************************
!
!  This code is designed to connect with LabView for a specific wave tank test case and likely will not work for other purposes.
!
!  For this test, a physical platform is deployed in a wave tank with cable acutators that are controlled through LabView.  This
!  module is called to provide some loads that are not present in the physical tank setup.  These include the following:
!     - rotor loading from a fixed RPM MHK rotor from AeroDyn.  This is calculated from either steady current provided by SeaState
!     - Mooring loads from MoorDyn
!
!
!**********************************************************************************************************************************
MODULE WaveTank_IO

   use NWTC_Library
   use NWTC_IO
   use WaveTank_Types

   implicit none

contains

subroutine ParseInputFile(FileInfo_In, SimSettings, ErrStat, ErrMsg)
   type(FileInfoType),        intent(in   )  :: FileInfo_In       !< The derived type for holding the file information.
   type(SimSettingsType),     intent(  out)  :: SimSettings
   integer(IntKi),            intent(  out)  :: ErrStat
   character(*),              intent(  out)  :: ErrMsg

   ! Local variables
   integer                                   :: CurLine
   character(1024), target                   :: TmpPath
   character(1024)                           :: FileName
   integer(IntKi)                            :: ErrStat2             ! local status of error message
   character(ErrMsgLen)                      :: ErrMsg2              ! local error message if errStat /= ErrID_None
   character(*), parameter                   :: RoutineName = 'WaveTankTesting.ParseInputFile'

   ErrStat = ErrID_None
   ErrMsg  = " "

   CurLine = 1
   ! Separator/header line skipped automatically
   ! ----- Simulation control -------------
   call ParseVar( FileInfo_In, CurLine, 'DT',               SimSettings%Sim%DT,                ErrStat2, ErrMsg2); if(Failed()) return;  ! timestep (unused)
   call ParseVar( FileInfo_In, CurLine, 'TMax',             SimSettings%Sim%TMax,              ErrStat2, ErrMsg2); if(Failed()) return;  ! Max sim time (unused)
   call ParseVar( FileInfo_In, CurLine, 'MHK',              SimSettings%Sim%MHK,               ErrStat2, ErrMsg2); if(Failed()) return;  ! MHK turbine type (switch) {0=Not an MHK turbine; 1=Fixed MHK turbine; 2=Floating MHK turbine}
   call ParseVar( FileInfo_In, CurLine, 'InterpOrd',        SimSettings%Sim%InterpOrd,         ErrStat2, ErrMsg2); if(Failed()) return;  ! Interpolation order (unused)
   call ParseVar( FileInfo_In, CurLine, 'DebugLevel',       SimSettings%Sim%DebugLevel,        ErrStat2, ErrMsg2); if(Failed()) return;  ! 0: none, 1: I/O summary, 2: +positions/orientations passed, 3:, 4: +all meshes
   call ParseVar( FileInfo_In, CurLine, 'OutRootName',      SimSettings%Sim%OutRootName,       ErrStat2, ErrMsg2); if(Failed()) return;  ! Root name for any summary or other files
   ! -------- Environment ----------------
   call ParseVar( FileInfo_In, CurLine, 'Gravity',          SimSettings%Env%Gravity,           ErrStat2, ErrMsg2); if(Failed()) return;  ! Gravitational acceleration (m/s^2)
   call ParseVar( FileInfo_In, CurLine, 'WtrDens',          SimSettings%Env%WtrDens,           ErrStat2, ErrMsg2); if(Failed()) return;  ! Water density (kg/m^3)
   call ParseVar( FileInfo_In, CurLine, 'WtrVisc',          SimSettings%Env%WtrVisc,           ErrStat2, ErrMsg2); if(Failed()) return;  ! Kinematic viscosity of working fluid (m^2/s)
   call ParseVar( FileInfo_In, CurLine, 'SpdSound',         SimSettings%Env%SpdSound,          ErrStat2, ErrMsg2); if(Failed()) return;  ! Speed of sound in working fluid (m/s)
   call ParseVar( FileInfo_In, CurLine, 'Patm',             SimSettings%Env%Patm,              ErrStat2, ErrMsg2); if(Failed()) return;  ! Atmospheric pressure (Pa) [used only for an MHK turbine cavitation check]
   call ParseVar( FileInfo_In, CurLine, 'Pvap',             SimSettings%Env%Pvap,              ErrStat2, ErrMsg2); if(Failed()) return;  !  Vapour pressure of working fluid (Pa) [used only for an MHK turbine cavitation check]
   call ParseVar( FileInfo_In, CurLine, 'WtrDpth',          SimSettings%Env%WtrDpth,           ErrStat2, ErrMsg2); if(Failed()) return;  ! Water depth (m)
   call ParseVar( FileInfo_In, CurLine, 'MSL2SWL',          SimSettings%Env%MSL2SWL,           ErrStat2, ErrMsg2); if(Failed()) return;  ! Offset between still-water level and mean sea level (m) [positive upward]
   ! -------- SeaState -------------------
   ! -------- MoorDyn --------------------
   ! -------- AeroDyn + InflowWind -------
   ! -------- Turbine Configuration ------
   call ParseVar( FileInfo_In, CurLine, 'NumBl',            SimSettings%TCfg%NumBl,            ErrStat2, ErrMsg2); if(Failed()) return;  ! Number of blades (-)
   call ParseVar( FileInfo_In, CurLine, 'TipRad',           SimSettings%TCfg%TipRad,           ErrStat2, ErrMsg2); if(Failed()) return;  ! The distance from the rotor apex to the blade tip (meters)
   call ParseVar( FileInfo_In, CurLine, 'HubRad',           SimSettings%TCfg%HubRad,           ErrStat2, ErrMsg2); if(Failed()) return;  ! The distance from the rotor apex to the blade root (meters)
   call ParseVar( FileInfo_In, CurLine, 'PreCone',          SimSettings%TCfg%PreCone,          ErrStat2, ErrMsg2); if(Failed()) return;  ! Blade cone angle (degrees)
   call ParseVar( FileInfo_In, CurLine, 'OverHang',         SimSettings%TCfg%OverHang,         ErrStat2, ErrMsg2); if(Failed()) return;  ! Distance from yaw axis to rotor apex [3 blades] or teeter pin [2 blades] (meters)
   call ParseVar( FileInfo_In, CurLine, 'ShftGagL',         SimSettings%TCfg%ShftGagL,         ErrStat2, ErrMsg2); if(Failed()) return;  ! Distance from rotor apex [3 blades] or teeter pin [2 blades] to shaft strain gages [positive for upwind rotors] (meters)
   call ParseVar( FileInfo_In, CurLine, 'ShftTilt',         SimSettings%TCfg%ShftTilt,         ErrStat2, ErrMsg2); if(Failed()) return;  ! Rotor shaft tilt angle (degrees)
   call ParseVar( FileInfo_In, CurLine, 'Twr2Shft',         SimSettings%TCfg%Twr2Shft,         ErrStat2, ErrMsg2); if(Failed()) return;  ! Vertical distance from the tower-top to the rotor shaft (meters)
   call ParseVar( FileInfo_In, CurLine, 'TowerHt',          SimSettings%TCfg%TowerHt,          ErrStat2, ErrMsg2); if(Failed()) return;  ! Height of tower relative MSL
   call ParseVar( FileInfo_In, CurLine, 'TowerBsHt',        SimSettings%TCfg%TowerBsHt,        ErrStat2, ErrMsg2); if(Failed()) return;  ! Height of tower base relative to ground level [onshore], MSL [floating MHK] (meters)
   call ParseAry( FileInfo_In, CurLine, 'PtfmRef',          SimSettings%TCfg%PtfmRef,    3,    ErrStat2, ErrMsg2); if(Failed()) return;  ! Location of platform reference point, relative to MSL.  Motions and loads all connect to this point
   ! -------- Turbine Operating Point ----
   call ParseVar( FileInfo_In, CurLine, 'RotSpeed',         SimSettings%TOp%RotSpeed,          ErrStat2, ErrMsg2); if(Failed()) return;  ! Rotational speed of rotor in rotor coordinates (rpm)
   call ParseVar( FileInfo_In, CurLine, 'NacYaw',           SimSettings%TOp%NacYaw,            ErrStat2, ErrMsg2); if(Failed()) return;  ! Initial or fixed nacelle-yaw angle (degrees)
   call ParseVar( FileInfo_In, CurLine, 'BldPitch',         SimSettings%TOp%BldPitch,          ErrStat2, ErrMsg2); if(Failed()) return;  ! Blade 1 pitch (deg)
   ! -------- Output ---------------------
   call ParseVar( FileInfo_In, CurLine, 'SendScreenToFile', SimSettings%Outs%SendScreenToFile, ErrStat2, ErrMsg2); if(Failed()) return;  ! send to file <OutRootName>.screen.log if true
   call ParseVar( FileInfo_In, CurLine, 'OutFile',          SimSettings%Outs%OutFile,          ErrStat2, ErrMsg2); if(Failed()) return;  ! 0: no output file of channels, 1: output file in text format (at default DT) 
   call ParseVar( FileInfo_In, CurLine, 'OutFmt',           SimSettings%Outs%OutFmt,           ErrStat2, ErrMsg2); if(Failed()) return;  ! Format used for text tabular output, excluding the time channel. (quoted string)
   ! -------- VTK output -----------------
   call ParseVar( FileInfo_In, CurLine, 'WrVTK_Dir',        SimSettings%Viz%WrVTK_Dir,         ErrStat2, ErrMsg2); if(Failed()) return;  ! output directory for visualization
   call ParseVar( FileInfo_In, CurLine, 'WrVTK',            SimSettings%Viz%WrVTK,             ErrStat2, ErrMsg2); if(Failed()) return;  ! VTK visualization data output: (switch) {0=none; 1=initialization data only; 2=animation; 3=mode shapes}
   call ParseVar( FileInfo_In, CurLine, 'WrVTK_type',       SimSettings%Viz%WrVTK_type,        ErrStat2, ErrMsg2); if(Failed()) return;  ! Type of VTK visualization data: (switch) {1=surfaces; 2=basic meshes (lines/points); 3=all meshes (debug)} [unused if WrVTK=0]
   call ParseVar( FileInfo_In, CurLine, 'WrVTK_DT',         SimSettings%Viz%WrVTK_DT,          ErrStat2, ErrMsg2); if(Failed()) return;  ! DT for writing VTK files
   call ParseAry( FileInfo_In, CurLine, 'VTKNacDim',        SimSettings%Viz%VTKNacDim,   6,    ErrStat2, ErrMsg2); if(Failed()) return;  !  Nacelle dimension passed in for VTK surface rendering [0,y0,z0,Lx,Ly,Lz] (m)

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine



subroutine ValidateInputFile(SimSettings, ErrStat, ErrMsg)
   type(SimSettingsType),     intent(inout)  :: SimSettings
   integer(IntKi),            intent(  out)  :: ErrStat
   character(*),              intent(  out)  :: ErrMsg
   integer(IntKi)                            :: ErrStat2
   character(ErrMsgLen)                      :: ErrMsg2
   character(*), parameter                   :: RoutineName = 'WaveTankTesting.ValidateInputFile'

   ErrStat = ErrID_None
   ErrMsg  = ""

   !------------------------
   ! Sim Control
   !------------------------
   if (SimSettings%Sim%MHK /= 2_c_int) call SetErrStat(ErrID_Fatal, "WaveTank module only works for floating MHK turbines at present (MHK=2).",ErrStat,ErrMsg,RoutineName)

   !------------------------
   ! Environment
   !------------------------

   !------------------------
   ! Turbine Config
   !------------------------

   !------------------------
   ! Turbine Operating point
   !------------------------

   !------------------------
   ! Output
   !------------------------

   !------------------------
   ! VTK
   !------------------------

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine

END MODULE WaveTank_IO
