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
MODULE WaveTankTesting

   use ISO_C_BINDING
   use NWTC_Library
   use SeaState_C_Binding, ONLY: SeaSt_C_PreInit, SeaSt_C_Init, SeaSt_C_CalcOutput, SeaSt_C_End, MaxOutPts, SeaSt_C_GetWaveFieldPointer
   use SeaSt_WaveField_Types, ONLY: SeaSt_WaveFieldType
   use AeroDyn_Inflow_C_BINDING, ONLY: ADI_C_PreInit, ADI_C_SetupRotor, ADI_C_Init, ADI_C_End, MaxADIOutputs, ADI_C_SetRotorMotion, ADI_C_UpdateStates, ADI_C_CalcOutput, ADI_C_GetRotorLoads
   use MoorDyn_C, ONLY: MD_C_Init, MD_C_End, MD_C_SetWaveFieldData, MD_C_UpdateStates, MD_C_CalcOutput
   use NWTC_C_Binding, ONLY: IntfStrLen, SetErrStat_C, SetErrStat_F2C, ErrMsgLen_C, StringConvert_F2C, FileNameFromCString, AbortErrLev_C
   use WaveTank_Types
   use WaveTank_IO

   implicit none
   save

   public :: WaveTank_Init
   public :: WaveTank_CalcStep
   public :: WaveTank_End

   ! output to screen or to file
   integer(IntKi) :: ScreenLogOutput_Un = -1
   character(1024):: ScreenLogOutput_File

   ! Simulation data storage
   type(SimSettingsType), target           :: SimSettings

   ! file output

!FIXME: replace all this with meshes
!   REAL(C_FLOAT), DIMENSION(3,3) :: FloaterPositions = 0.0_C_FLOAT
!   REAL(C_FLOAT), DIMENSION(2,6) :: FloaterVelocities = 0.0_C_FLOAT
!   REAL(C_FLOAT), DIMENSION(1,6) :: FloaterAccelerations = 0.0_C_FLOAT
!   REAL(C_FLOAT), DIMENSION(3,3) :: NacellePositions = 0.0_C_FLOAT
!   REAL(C_FLOAT), DIMENSION(2,6) :: NacelleVelocities = 0.0_C_FLOAT
!   REAL(C_FLOAT), DIMENSION(1,6) :: NacelleAccelerations = 0.0_C_FLOAT
!   REAL(C_FLOAT), ALLOCATABLE :: BladeRootPositions(:,:)
!   REAL(C_FLOAT), ALLOCATABLE :: BladeRootVelocities(:,:)
!   REAL(C_FLOAT), ALLOCATABLE :: BladeRootAccelerations(:,:)
!   REAL(C_FLOAT), ALLOCATABLE :: BladeMeshPositions(:,:)
!   REAL(C_FLOAT), ALLOCATABLE :: BladeMeshVelocities(:,:)
!   REAL(C_FLOAT), ALLOCATABLE :: BladeMeshAccelerations(:,:)
!FIXME: this is temporary until meshes are properly setup.  Move all this to some temporary storage that gets allocated on init
real(c_double) :: tmpIdent9(9) = (/ 1.0_c_float, 0.0_c_float, 0.0_c_float,  0.0_c_float, 1.0_c_float, 0.0_c_float,  0.0_c_float, 0.0_c_float, 1.0_c_float /)
real(c_double) :: tmpBldRootOri(18) = (/ 1.0_c_float, 0.0_c_float, 0.0_c_float,  0.0_c_float, 1.0_c_float, 0.0_c_float,  0.0_c_float, 0.0_c_float, 1.0_c_float, &
                                         1.0_c_float, 0.0_c_float, 0.0_c_float,  0.0_c_float, 1.0_c_float, 0.0_c_float,  0.0_c_float, 0.0_c_float, 1.0_c_float /)
real(c_double) :: tmpInitMeshOri(18)
real(c_float)  :: tmpBldRootPos(6), tmpHubPos(3), tmpNacPos(3), tmpInitMeshPos(6)
integer(c_int) :: tmpNumMeshPts = 2
integer(c_int) :: tmpMeshPtToBladeNum(2) = (/ 1_c_int, 2_c_int /)

!FIXME: this is temporary -- move to a derived type for clarity
   INTEGER(C_INT) :: SS_NumChannels_C
   INTEGER(C_INT) :: MD_NumChannels_C
   INTEGER(C_INT) :: ADI_NumChannels_C
   INTEGER(C_INT) :: NumBlades_C
   INTEGER(C_INT) :: NumMeshPts_C

CONTAINS

SUBROUTINE WaveTank_Init(  &
   WT_InputFile_C,         &
   RootName_C,             &
   VTKdir_C,               &
   ErrStat_C,              &
   ErrMsg_C                &
) BIND (C, NAME='WaveTank_Init')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: WaveTank_Init
!GCC$ ATTRIBUTES DLLEXPORT :: WaveTank_Init
#endif

   character(c_char),          intent(in   ) :: WT_InputFile_C(IntfStrLen)
   character(kind=c_char),     intent(  out) :: RootName_C(IntfStrLen)
   character(kind=c_char),     intent(  out) :: VTKdir_C(IntfStrLen)
   integer(c_int),             intent(  out) :: ErrStat_C
   character(kind=c_char),     intent(  out) :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   integer(c_int)                          :: ErrStat_C2
   character(kind=c_char, len=ErrMsgLen_C) :: ErrMsg_C2
   integer(IntKi)                          :: ErrStat_F2
   character(ErrMsgLen)                    :: ErrMsg_F2
   character(1024)                         :: InputFile
   integer(IntKi)                          :: i
   type(FileInfoType)                      :: FileInfo_In   !< The derived type for holding the full input file for parsing -- we may pass this in the future
   real(c_float)                           :: InitPtfmPosOri(6)
   ! local C variables for transferring names
   character(kind=c_char) :: WrVTK_Dir_C(IntfStrLen)
   character(kind=c_char) :: OutRootName_C(IntfStrLen)

   ! The length of these arrays much match what is set in the corresponding C binding modules
   character(kind=c_char) :: SS_OutputChannelNames_C(ChanLen*MaxOutPts+1)
   character(kind=c_char) :: SS_OutputChannelUnits_C(ChanLen*MaxOutPts+1)
   character(kind=c_char) :: MD_OutputChannelNames_C(100000)
   character(kind=c_char) :: MD_OutputChannelUnits_C(100000)
   character(kind=c_char) :: ADI_OutputChannelNames_C(ChanLen*MaxADIOutputs+1)
   character(kind=c_char) :: ADI_OutputChannelUnits_C(ChanLen*MaxADIOutputs+1)

   ! Filename conversions -- read in as fortran strings, but sent to other modules as c_char arrays
   character(kind=c_char)         :: SS_InputFile_C(IntfStrLen)
   character(kind=c_char), target :: MD_InputFile_C(IntfStrLen)
   character(kind=c_char), target :: AD_InputFile_C(IntfStrLen)
   character(kind=c_char), target :: IfW_InputFile_C(IntfStrLen)


   ! Initialize error handling
   ErrStat_C = ErrID_None
   ErrMsg_C  = " "//C_NULL_CHAR

   InputFile = transfer(WT_InputFile_C, InputFile)
   i = index(InputFile, char(0))
   InputFile = InputFile(1:i)
   call ProcessComFile(InputFile,  FileInfo_In, ErrStat_F2, ErrMsg_F2); if (Failed()) return
   call ParseInputFile(FileInfo_In, SimSettings, ErrStat_F2, ErrMsg_F2); if (Failed()) return

   call ValidateInputFile(SimSettings, ErrStat_F2, ErrMsg_F2); if (Failed()) return

   ! transfer filenames for passing to modules
   SS_InputFile_C  = c_null_char
   MD_InputFile_C  = c_null_char
   AD_InputFile_C  = c_null_char
   IfW_InputFile_C = c_null_char
   SS_InputFile_C  = transfer(trim(SimSettings%IptFile%SS_InputFile ), SS_InputFile_C )
   MD_InputFile_C  = transfer(trim(SimSettings%IptFile%MD_InputFile ), MD_InputFile_C )
   AD_InputFile_C  = transfer(trim(SimSettings%IptFile%AD_InputFile ), AD_InputFile_C )
   IfW_InputFile_C = transfer(trim(SimSettings%IptFile%IfW_InputFile), IfW_InputFile_C)

   ! return rootname
   RootName_C = c_null_char
   RootName_C = transfer(trim(SimSettings%Sim%OutRootName),RootName_C)

   ! If SendScreenToFile - send to file <OutRootName>.screen.log if true
   if (SimSettings%Outs%SendScreenToFile) then
      call GetNewUnit(ScreenLogOutput_Un, ErrStat_F2, ErrMsg_F2); if (Failed()) return
      ScreenLogOutput_File = trim(SimSettings%Sim%OutRootName)//'.screen.log'
      call OpenFOutFile(ScreenLogOutput_Un, ScreenLogOutput_File, ErrStat_F2, ErrMsg_F2); if (Failed()) return
      call SetConsoleUnit(ScreenLogOutput_Un)   ! this will redirect all screen output to a file instead
   endif

   ! debugging
   if (SimSettings%Sim%DebugLevel > 0_c_int) call ShowPassedData()

   ! VTK directory
   WrVTK_Dir_C = c_null_char
   WrVTK_Dir_C = transfer( trim(SimSettings%Viz%WrVTK_Dir), WrVTK_Dir_C )
   ! return VTKdir
   VTKdir_C = c_null_char
   if (SimSettings%Viz%WrVTK > 0_c_int) VTKdir_C = WrVTK_Dir_C


!FIXME: build struct model


   !------------------------------
   ! Setup and initialize SeaState
   !------------------------------
   call SeaSt_C_PreInit(            &
      SimSettings%Env%Gravity,      &
      SimSettings%Env%WtrDens,      &
      SimSettings%Env%WtrDpth,      &
      SimSettings%Env%MSL2SWL,      &
      SimSettings%Sim%DebugLevel,   &
      WrVTK_Dir_C,                  &
      SimSettings%Viz%WrVTK,        &
      SimSettings%Viz%WrVTK_DT,     &
      ErrStat_C2, ErrMsg_C2         &
   )
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'SeaSt_C_PreInit')
   if (ErrStat_C >= AbortErrLev_C) return

   OutRootName_C = transfer(trim(SimSettings%Sim%OutRootName)//'.SeaSt'//c_null_char, OutRootName_C)
   call SeaSt_C_Init(            &    
      SS_InputFile_C,            &
      OutRootName_C,             &
      1000_c_int,                & !FIXME: do we need the number of timesteps???
      SimSettings%Sim%DT,        &
      SS_NumChannels_C,          &
      SS_OutputChannelNames_C,   &
      SS_OutputChannelUnits_C,   &
      ErrStat_C2, ErrMsg_C2      &
   )
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'SeaSt_C_Init')
   if (ErrStat_C >= AbortErrLev_C) return

   !------------------------------
   ! Set the SeaState Wave Field pointer onto MoorDyn
   !------------------------------
   call WaveTank_SetWaveFieldPointer(ErrStat_C2, ErrMsg_C2)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_SetWaveFieldPointer')
   if (ErrStat_C >= AbortErrLev_C) return

   !------------------------------
   ! Setup and initialize MoorDyn
   !------------------------------
!FIXME: this interface will change!!! -- Split with PreInit
!FIXME: add WrVTK_Dir_C,  SimSettings%Viz%WrVTK, SimSettings%Viz%WrVTK_DT
   !SimSettings%TCfg%PtfmRef
!FIXME: 6 DOF with 3xPos, 3xEulerAngle
   InitPtfmPosOri = 0.0_c_float
   OutRootName_C = transfer(trim(SimSettings%Sim%OutRootName)//'.MD'//c_null_char, OutRootName_C)
   call MD_C_Init(                           &
      0_c_int,                               &   !< InputFilePassed: 0 for file, 1 for string
      c_loc(MD_InputFile_C(1)),              &
      int(IntfStrLen,c_int),                 &   !< InputFileStringLength_C
      SimSettings%Sim%DT,                    &
      SimSettings%Env%Gravity,               &
      SimSettings%Env%WtrDens,               &
      SimSettings%Env%WtrDpth,               &
      InitPtfmPosOri,                        &
      SimSettings%Sim%InterpOrd,             &
      MD_NumChannels_C,                      &
      MD_OutputChannelNames_C,               &
      MD_OutputChannelUnits_C,               &
      ErrStat_C2, ErrMsg_C2                  &
   )
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'MD_C_Init')
   if (ErrStat_C >= AbortErrLev_C) return

   !------------------------------
   ! Setup and initialize MoorDyn
   !------------------------------
   call ADI_C_PreInit(                       &
      1_c_int,                               &     ! only one turbine
      1_c_int,                               &     ! transpose DCM inside ADI (true)
      1_c_int,                               &     ! PointLoadOutput - use line to point load mapping -- necessary for mapping to blade root without an actual blade structure
      SimSettings%Env%Gravity,               &
      SimSettings%Env%WtrDens,               &
      SimSettings%Env%WtrVisc,               &
      SimSettings%Env%SpdSound,              &
      SimSettings%Env%Patm,                  &
      SimSettings%Env%Pvap,                  &
      SimSettings%Env%WtrDpth,               &
      SimSettings%Env%MSL2SWL,               &
      SimSettings%Sim%MHK,                   &
      SimSettings%Sim%DebugLevel,            &
      ErrStat_C2, ErrMsg_C2                  &
   )
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_PreInit')
   if (ErrStat_C >= AbortErrLev_C) return


   !------------------------------
   ! Setup and initialize MoorDyn
   !------------------------------
!FIXME: temporary location info until meshes set up
tmpNacPos = SimSettings%TCfg%TowerHt
tmpHubPos = tmpNacPos + SimSettings%TCfg%OverHang + SimSettings%TCfg%Twr2Shft     ! missing angles
tmpBldRootPos(1:3) = tmpHubPos   ! blade 1
tmpBldRootPos(4:6) = tmpHubPos   ! blade 2
tmpInitMeshPos = tmpBldRootPos
tmpInitMeshOri = tmpBldRootOri      !FIXME: blade pitch
   call ADI_C_SetupRotor(                    &
      1_c_int,                               &  ! iWT  -- turbine number (only one allowed)
      1_c_int,                               &  ! TurbineIsHAWT = true
      InitPtfmPosOri,                        &  ! TurbOrigin - FIXME: is this at the platform reference point? or tower base height?
      tmpHubPos          ,     &  ! HubPos
      tmpIdent9          ,     &  ! HubOri
      tmpNacPos          ,     &  ! NacPos
      tmpIdent9          ,     &  ! NacOri
      SimSettings%TCfg%NumBl,                &  ! NumBlades
      tmpBldRootPos      ,     &  ! BldRootPos
      tmpBldRootOri      ,     &  ! BldRootOri
      tmpNumMeshPts      ,     &  ! NumMeshPts
      tmpInitMeshPos     ,     &  ! InitMeshPos
      tmpInitMeshOri     ,     &  ! InitMeshOri
      tmpMeshPtToBladeNum,     &  ! MeshPtToBladeNum
      ErrStat_C2, ErrMsg_C2                   &
   )
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_SetupRotor')
   if (ErrStat_C >= AbortErrLev_C) return

   OutRootName_C = transfer(trim(SimSettings%Sim%OutRootName)//'.ADI'//c_null_char, OutRootName_C)
   call ADI_C_Init(                          &
      0,                                     &  ! ADinputFilePassed; 0 for file, 1 for string
      c_loc(AD_InputFile_C(1)),              &  ! ADinputFileString_C; Input file as a single string with lines delineated by C_NULL_CHAR
      IntfStrLen,                            &  ! ADinputFileStringLength_C; length of the input file string
      0,                                     &  ! IfWinputFilePassed; 0 for file, 1 for string
      c_loc(IfW_InputFile_C(1)),             &  ! IfWinputFileString_C; Input file as a single string with lines delineated by C_NULL_CHAR
      IntfStrLen,                            &  ! IfWinputFileStringLength_C; length of the input file string
      OutRootName_C,                         &  ! Root name to use for echo files and other
      WrVTK_Dir_C,                           &  ! vtk directory to use
      SimSettings%Sim%InterpOrd,             &  ! interpolation order for extrap/interp
      SimSettings%Sim%DT,                    &  ! DT for simulation (used in checks only)
      SimSettings%Sim%TMax,                  &  ! Max time for simulation (not used here)
      0_c_int,                               &  ! storeHHVel - Store hub height time series from IfW -- set to false since not used here
      SimSettings%Viz%WrVTK,                 &  ! VTK visualization data output: (switch) {0=none; 1=initialization data only; 2=animation; 3=mode shapes}
      SimSettings%Viz%WrVTK_Type,            &  ! Type of VTK visualization data: (switch) {1=surfaces; 2=basic meshes (lines/points); 3=all meshes (debug)} [unused if WrVTK=0]
      SimSettings%Viz%WrVTK_DT,              &  ! timestep of VTK writing
      SimSettings%Viz%VTKNacDim,             &  ! Nacelle dimension passed in for VTK surface rendering [0,y0,z0,Lx,Ly,Lz] (m)
      SimSettings%TCfg%HubRad,               &  ! Hub radius for VTK surface rendering
      1_c_int,                               &  ! wrOuts_C -- Write ADI output file -- hard code to true for now
      SimSettings%Sim%DT,                    &  ! Timestep to write output file from ADI
      ADI_NumChannels_C,                     &
      ADI_OutputChannelNames_C,              &
      ADI_OutputChannelUnits_C,              &
      ErrStat_C2, ErrMsg_C2                  &
   )
   CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_Init')
   IF (ErrStat_C >= AbortErrLev_C) RETURN

contains
   logical function Failed()
      call SetErrStat_F2C(ErrStat_F2, ErrMsg_F2, ErrStat_C, ErrMsg_C)
      Failed = ErrStat_C >= AbortErrLev_C
      if (Failed)    call Cleanup()
   end function Failed
   subroutine Cleanup()
      call NWTC_Library_Destroyfileinfotype(FileInfo_In, ErrStat_F2, ErrMsg_F2)  ! ignore error from this
      if (ScreenLogOutput_Un > 0)   close(ScreenLogOutput_Un)
!      if (ErrStat_C >= AbortErrLev_C) call DeallocEverything()
   end subroutine Cleanup
   subroutine ShowPassedData()
      character(IntfStrLen) :: tmpPath
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  WaveTank_Init")
      call WrScr("   --------------------------------------------------------")
      call WrScr("   WT_InputFile_C         -> "//trim(InputFile))
      call WrScr("-----------------------------------------------------------")
   end subroutine ShowPassedData
END SUBROUTINE WaveTank_Init

!subroutine DeallocEverything()
!end subroutine DeallocEverything



subroutine WaveTank_CalcStep( &
   time,                      &
   pos,                       &
   vel,                       &
   acc,                       &
   loads,                     &
   ErrStat_C,                 &
   ErrMsg_C                   &
) BIND (C, NAME='WaveTank_CalcStep')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: WaveTank_CalcStep
!GCC$ ATTRIBUTES DLLEXPORT :: WaveTank_CalcStep
#endif
   real(c_double),         intent(in   ) :: time
   real(c_float),          intent(in   ) :: pos(6)       ! [x,y,z,roll,pitch,yaw]
   real(c_float),          intent(in   ) :: vel(6)       ! [x_dot,y_dot,z_dot,roll_dot,pitch_dot,yaw_dot]
   real(c_float),          intent(in   ) :: acc(6)       ! [x_ddot,y_ddot,z_ddot,roll_ddot,pitch_ddot,yaw_ddot]
   real(c_float),          intent(  out) :: loads(6)     ! [Fx,Fy,Fz,Mx,My,Mz]
   integer(c_int),         intent(  out) :: ErrStat_C
   character(c_char),      intent(  out) :: ErrMsg_C(ErrMsgLen_C)
   integer(c_int)                        :: ErrStat_C2(ErrMsgLen_C)
   character(c_char)                     :: ErrMsg_C2

   ! Initialize error handling
   ErrStat_C = ErrID_None
   ErrMsg_C  = " "//C_NULL_CHAR

   ! zero loads in case of error
   loads = 0.0_c_float

!FIXME: setup new method for handling this
!FIXME: turb off HHVel
!   real(c_float) :: md_outputs(MD_NumChannels_C)
!   real(c_float) :: adi_outputs(ADI_NumChannels_C)
!   real(c_float) :: ADI_HHVel_C(3)                    !< Wind speed array [Vx,Vy,Vz]                      -- (m/s) (global)

    ! ! ADI
    ! ! SetRotorMotion
!    REAL(C_FLOAT)      :: ADI_HubPos_C( 3 )                 !< Hub position
!    REAL(C_DOUBLE)     :: ADI_HubOri_C( 9 )                 !< Hub orientation
!    REAL(C_FLOAT)      :: ADI_HubVel_C( 6 )                 !< Hub velocity
!    REAL(C_FLOAT)      :: ADI_HubAcc_C( 6 )                 !< Hub acceleration
!    REAL(C_FLOAT)      :: ADI_NacPos_C( 3 )                 !< Nacelle position
!    REAL(C_DOUBLE)     :: ADI_NacOri_C( 9 )                 !< Nacelle orientation
!    REAL(C_FLOAT)      :: ADI_NacVel_C( 6 )                 !< Nacelle velocity
!    REAL(C_FLOAT)      :: ADI_NacAcc_C( 6 )                 !< Nacelle acceleration
!    REAL(C_FLOAT)      :: ADI_BldRootPos_C( 3*NumBlades_C ) !< Blade root positions
!    REAL(C_DOUBLE)     :: ADI_BldRootOri_C( 9*NumBlades_C ) !< Blade root orientations
!    REAL(C_FLOAT)      :: ADI_BldRootVel_C( 6*NumBlades_C ) !< Blade root velocities
!    REAL(C_FLOAT)      :: ADI_BldRootAcc_C( 6*NumBlades_C ) !< Blade root accelerations
!    ! Blade mesh nodes
!    REAL(C_FLOAT)      :: ADI_MeshPos_C( 3*NumMeshPts_C )   !< A 3xNumMeshPts_C array [x,y,z]
!    REAL(C_DOUBLE)     :: ADI_MeshOri_C( 9*NumMeshPts_C )   !< A 9xNumMeshPts_C array [r11,r12,r13,r21,r22,r23,r31,r32,r33]
!    REAL(C_FLOAT)      :: ADI_MeshVel_C( 6*NumMeshPts_C )   !< A 6xNumMeshPts_C array [x,y,z]
!    REAL(C_FLOAT)      :: ADI_MeshAcc_C( 6*NumMeshPts_C )   !< A 6xNumMeshPts_C array [x,y,z]


!   ! Shift the positions and velocities over one index
!   FloaterPositions(1,:) = FloaterPositions(2,:)
!   FloaterPositions(2,:) = FloaterPositions(3,:)
!   NacellePositions(1,:) = NacellePositions(2,:)
!   NacellePositions(2,:) = NacellePositions(3,:)
!   BladeRootPositions(1,:) = BladeRootPositions(2,:)
!   BladeRootPositions(2,:) = BladeRootPositions(3,:)
!   BladeMeshPositions(1,:) = BladeMeshPositions(2,:)
!   BladeMeshPositions(2,:) = BladeMeshPositions(3,:)
!   FloaterVelocities(1,:) = FloaterVelocities(2,:)
!   NacelleVelocities(1,:) = NacelleVelocities(2,:)
!   BladeRootVelocities(1,:) = BladeRootVelocities(2,:)
!   BladeMeshVelocities(1,:) = BladeMeshVelocities(2,:)

!   ! Load the new positions
!   FloaterPositions(3,:) = (/ positions_x, positions_y, positions_z /)
!   DeltaS = FloaterPositions(3,:) - FloaterPositions(2,:)
!   ! TODO: rigid body rotation is missing on moment arm
!   NacellePositions(3,:) = NacellePositions(2,:) + DeltaS

    ! TODO:
    ! - Create mesh for platform point (WaveTank_Init)
    ! - Create mesh for hub position (WaveTank_Init)
    !   - Need rotational velocity for hub
    ! - Mesh map hub to blade roots
    ! - Mesh map blade roots to blade nodes (blade pitch?)
    ! - Create entire structural model (copy from AeroDyn driver)

    ! Stride
    ! Lower bound: (I-1)*3+1 is the first of the three position components for the current blade
    !               +1 is because Fortran starts indexing at 1
    ! Upper bound: (I-1)*3+1+2 is the last of the three position components for the current blade
    !               +2 is because Fortran includes the last index in the range
!   DO I=1,NumBlades_C
!       I0 = (I-1)*3+1
!       I1 = (I-1)*3+1+2
!       BladeRootPositions(3,I0:I1) = BladeRootPositions(2,I0:I1) + DeltaS
!   END DO
!   DO I=1,NumMeshPts_C
!       I0 = (I-1)*3+1
!       I1 = (I-1)*3+1+2
!       BladeMeshPositions(3,I0:I1) = BladeMeshPositions(2,I0:I1) + DeltaS
!   END DO

!   ! Calculate velocities and acceleration
!   FloaterVelocities(1,:) = (/ (FloaterPositions(2,:) - FloaterPositions(1,:)) / REAL(DT, C_FLOAT), 0.0_C_FLOAT, 0.0_C_FLOAT, 0.0_C_FLOAT /)
!   FloaterVelocities(2,:) = (/ (FloaterPositions(3,:) - FloaterPositions(2,:)) / REAL(DT, C_FLOAT), 0.0_C_FLOAT, 0.0_C_FLOAT, 0.0_C_FLOAT /)
!   FloaterAccelerations(1,:) = (/ (FloaterVelocities(2,:) - FloaterVelocities(1,:)) / REAL(DT, C_FLOAT) /)

!   NacelleVelocities(1,:) = (/ (NacellePositions(2,:) - NacellePositions(1,:)) / REAL(DT, C_FLOAT), 0.0_C_FLOAT, 0.0_C_FLOAT, 0.0_C_FLOAT /)
!   NacelleVelocities(2,:) = (/ (NacellePositions(3,:) - NacellePositions(2,:)) / REAL(DT, C_FLOAT), 0.0_C_FLOAT, 0.0_C_FLOAT, 0.0_C_FLOAT /)
!   NacelleAccelerations(1,:) = (/ (NacelleVelocities(2,:) - NacelleVelocities(1,:)) / REAL(DT, C_FLOAT) /)

!   DO I=1,NumBlades_C
!       I0 = (I-1)*6+1
!       I1 = (I-1)*6+1+5
!       BladeRootVelocities(1,I0:I1) = (/ (BladeRootPositions(2,(I-1)*3+1:(I-1)*3+1+2) - BladeRootPositions(1,(I-1)*3+1:(I-1)*3+1+2)) / REAL(DT, C_FLOAT), 0.0_C_FLOAT, 0.0_C_FLOAT, 0.0_C_FLOAT /)
!       BladeRootVelocities(2,I0:I1) = (/ (BladeRootPositions(3,(I-1)*3+1:(I-1)*3+1+2) - BladeRootPositions(2,(I-1)*3+1:(I-1)*3+1+2)) / REAL(DT, C_FLOAT), 0.0_C_FLOAT, 0.0_C_FLOAT, 0.0_C_FLOAT /)
!       BladeRootAccelerations(1,I0:I1) = (BladeRootVelocities(2,I0:I1) - BladeRootVelocities(1,I0:I1)) / REAL(DT, C_FLOAT)
!   END DO

!   DO I=1,NumMeshPts_C
!       I0 = (I-1)*6+1
!       I1 = (I-1)*6+1+5
!       BladeMeshVelocities(1,I0:I1) = (/ (BladeMeshPositions(2,(I-1)*3+1:(I-1)*3+1+2) - BladeMeshPositions(1,(I-1)*3+1:(I-1)*3+1+2)) / REAL(DT, C_FLOAT), 0.0_C_FLOAT, 0.0_C_FLOAT, 0.0_C_FLOAT /)
!       BladeMeshVelocities(2,I0:I1) = (/ (BladeMeshPositions(3,(I-1)*3+1:(I-1)*3+1+2) - BladeMeshPositions(2,(I-1)*3+1:(I-1)*3+1+2)) / REAL(DT, C_FLOAT), 0.0_C_FLOAT, 0.0_C_FLOAT, 0.0_C_FLOAT /)
!       BladeMeshAccelerations(1,I0:I1) = (BladeMeshVelocities(2,I0:I1) - BladeMeshVelocities(1,I0:I1)) / REAL(DT, C_FLOAT)
!   END DO

!   ! Get loads from MoorDyn
!   ! NOTE: MD_C_UpdateStates and MD_C_CalcOutput do not use the positions, velocities, and accelerations.
!   !       They're passed here just for consistency, but we should not let that interface drive
!   !       the design of this module.
!   ! TODO: get angles from mesh
!   CALL MD_C_UpdateStates(                 &
!       time,                               &
!       REAL(time + DT, C_DOUBLE),          &
!       (/ FloaterPositions(3,:), 0.0_C_FLOAT, 0.0_C_FLOAT, 0.0_C_FLOAT /), &
!       (/ FloaterVelocities(2,:) /),       &
!       (/ FloaterAccelerations(1,:) /),    &
!       ErrStat_C2, ErrMsg_C2               &
!   )
!   CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'MD_C_UpdateStates')
!   IF (ErrStat_C >= AbortErrLev_C) RETURN

!   ! TODO: get angles from mesh
!   CALL MD_C_CalcOutput(                   &
!       time,                               &
!       (/ FloaterPositions(3,:), 0.0_C_FLOAT, 0.0_C_FLOAT, 0.0_C_FLOAT /), &
!       (/ FloaterVelocities(2,:) /),       &
!       (/ FloaterAccelerations(1,:) /),    &
!       MD_Forces_C,                        &
!       md_outputs,                         &
!       ErrStat_C2, ErrMsg_C2               &
!   )
!   CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'MD_C_CalcOutput')
!   IF (ErrStat_C >= AbortErrLev_C) RETURN

!   ! Get loads from ADI

!   ! All components are rigidly connected so they share velocities, accelerations, and orientation
!   ! Positions are set by the input file geometry and the calling code
!   ADI_HubPos_C = FloaterPositions(3,:)
!   ADI_HubOri_C = floater_rotation_matrix
!   ADI_HubVel_C = FloaterVelocities(2,:)
!   ADI_HubAcc_C = FloaterAccelerations(1,:)

!   ADI_NacPos_C = NacellePositions(3,:)
!   ADI_NacOri_C = floater_rotation_matrix
!   ADI_NacVel_C = NacelleVelocities(2,:)
!   ADI_NacAcc_C = NacelleAccelerations(1,:)

!   ADI_BldRootPos_C = BladeRootPositions(3,:)
!   ADI_BldRootOri_C = blade_rotation_matrix
!   ADI_BldRootVel_C = BladeRootVelocities(2,:)
!   ADI_BldRootAcc_C = BladeRootAccelerations(1,:)

!   ADI_MeshPos_C = BladeMeshPositions(3,:)
!   ADI_MeshOri_C = blade_rotation_matrix
!   ADI_MeshVel_C = BladeMeshVelocities(2,:)
!   ADI_MeshAcc_C = BladeMeshAccelerations(1,:)

!   CALL ADI_C_SetRotorMotion(              &
!       iWT_c,                              & !< Wind turbine / rotor number
!       ADI_HubPos_C,                       & !< Hub position
!       ADI_HubOri_C,                       & !< Hub orientation
!       ADI_HubVel_C,                       & !< Hub velocity
!       ADI_HubAcc_C,                       & !< Hub acceleration
!       ADI_NacPos_C,                       & !< Nacelle position
!       ADI_NacOri_C,                       & !< Nacelle orientation
!       ADI_NacVel_C,                       & !< Nacelle velocity
!       ADI_NacAcc_C,                       & !< Nacelle acceleration
!       ADI_BldRootPos_C,                   & !< Blade root positions, 3xNumBlades_C
!       ADI_BldRootOri_C,                   & !< Blade root orientations, 9xNumBlades_C
!       ADI_BldRootVel_C,                   & !< Blade root velocities, 6xNumBlades_C
!       ADI_BldRootAcc_C,                   & !< Blade root accelerations, 6xNumBlades_C
!       NumMeshPts_C,                       & !< Number of mesh points we are transferring motions to and output loads to
!       ADI_MeshPos_C,                      & !< A 3xNumMeshPts_C array [x,y,z]
!       ADI_MeshOri_C,                      & !< A 9xNumMeshPts_C array [r11,r12,r13,r21,r22,r23,r31,r32,r33]
!       ADI_MeshVel_C,                      & !< A 6xNumMeshPts_C array [x,y,z]
!       ADI_MeshAcc_C,                      & !< A 6xNumMeshPts_C array [x,y,z]
!       ErrStat_C2, ErrMsg_C2               &
!   )
!   CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_SetRotorMotion')
!   IF (ErrStat_C >= AbortErrLev_C) RETURN

!   CALL ADI_C_UpdateStates(                &
!       time,                               &
!       REAL(time + DT, C_DOUBLE),          &
!       ErrStat_C2, ErrMsg_C2               &
!   )
!   CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_UpdateStates')
!   IF (ErrStat_C >= AbortErrLev_C) RETURN

!   CALL ADI_C_CalcOutput(                  &
!       time,                               &
!       adi_outputs,                        &
!       ErrStat_C2, ErrMsg_C2               &
!   )
!   CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_CalcOutput')
!   IF (ErrStat_C >= AbortErrLev_C) RETURN

!   CALL ADI_C_GetRotorLoads(               &
!       iWT_c,                              & !< Wind turbine / rotor number
!       NumMeshPts_C,                       & !< Number of mesh points we are transfering motions to and output loads to
!       ADI_MeshFrc_C,                      & !< A 6xNumMeshPts_C array [Fx,Fy,Fz,Mx,My,Mz]       -- forces and moments (global)
!       ADI_HHVel_C,                        & !< Wind speed array [Vx,Vy,Vz]                      -- (m/s) (global)
!       ErrStat_C2, ErrMsg_C2               &
!   )
!   CALL SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_GetRotorLoads')
!   IF (ErrStat_C >= AbortErrLev_C) RETURN

end subroutine


SUBROUTINE WaveTank_End(ErrStat_C, ErrMsg_C) bind (C, NAME="WaveTank_End")
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: WaveTank_End
!GCC$ ATTRIBUTES DLLEXPORT :: WaveTank_End
#endif

   integer(c_int),         intent(  out) :: ErrStat_C
   character(kind=c_char), intent(  out) :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   integer(c_int)                          :: ErrStat_C2
   character(kind=c_char, len=ErrMsgLen_C) :: ErrMsg_C2
   integer(IntKi)                          :: ErrStat_F2
   character(ErrMsgLen)                    :: ErrMsg_F2

   ErrStat_C = ErrID_None
   ErrMsg_C  = " "//C_NULL_CHAR

   ! in case we were writing to a file instead of the screen
   if (ScreenLogOutput_Un > 0)   close(ScreenLogOutput_Un)

   call MD_C_END(ErrStat_C2, ErrMsg_C2)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'MD_C_END')

   call SeaSt_C_END(ErrStat_C2, ErrMsg_C2)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'SeaSt_C_END')

   call ADI_C_END(ErrStat_C2, ErrMsg_C2)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_END')

!FIXME: close output file here

END SUBROUTINE

subroutine WaveTank_SetWaveFieldPointer(ErrStat_C, ErrMsg_C) bind (C, NAME="WaveTank_SetWaveFieldPointer")
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: WaveTank_SetWaveFieldPointer
!GCC$ ATTRIBUTES DLLEXPORT :: WaveTank_SetWaveFieldPointer
#endif
   integer(c_int),            intent(  out) :: ErrStat_C
   character(kind=c_char),    intent(  out) :: ErrMsg_C(ErrMsgLen_C)
   integer(c_int)                           :: ErrStat_C2
   character(kind=c_char,  len=ErrMsgLen_C) :: ErrMsg_C2

   ! Set the SeaState FlowField pointer onto MoorDyn
   type(c_ptr)                              :: WaveFieldPointer_C
   type(SeaSt_WaveFieldType),       pointer :: WaveFieldPointer_F => NULL()      ! used only in sanity check

   ErrStat_C = ErrID_None
   ErrMsg_C  = " "//C_NULL_CHAR

   call SeaSt_C_GetWaveFieldPointer(WaveFieldPointer_C, ErrStat_C2, ErrMsg_C2)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_SetWaveFieldPointer')
   if (ErrStat_C >= AbortErrLev_C) return

   call C_F_POINTER(WaveFieldPointer_C, WaveFieldPointer_F)
   ! Verify that the data in the WaveField pointer has been set
   if (WaveFieldPointer_F%WtrDpth == 0) then
       ErrStat_C2 = ErrID_Fatal
       ErrMsg_C2 = "SeaState WaveFieldPointer is WtrDpth is 0.0, so it it probably not initialized."
       call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_SetWaveFieldPointer')
       return
   endif

!FIXME: MD needs error handling on this
   call MD_C_SetWaveFieldData(WaveFieldPointer_C) ! ErrStat_C2, ErrMsg_C2
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_SetWaveFieldPointer')
   if (ErrStat_C >= AbortErrLev_C) return

   ! Probably doesn't matter, but clear this pointer just in case
   WaveFieldPointer_F => NULL()

end subroutine

END MODULE WaveTankTesting
