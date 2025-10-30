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
   use SeaState_C_Binding, ONLY: SeaSt_C_PreInit, SeaSt_C_Init, SeaSt_C_CalcOutput, SeaSt_C_End, MaxOutPts, SeaSt_C_GetWaveFieldPointer, SeaSt_C_GetSurfElev
   use SeaSt_WaveField_Types, ONLY: SeaSt_WaveFieldType
   use AeroDyn_Inflow_C_BINDING, ONLY: ADI_C_PreInit, ADI_C_SetupRotor, ADI_C_Init, ADI_C_End, MaxADIOutputs, ADI_C_SetRotorMotion, ADI_C_UpdateStates, ADI_C_CalcOutput, ADI_C_GetRotorLoads
   use MoorDyn_C, ONLY: MD_C_Init, MD_C_End, MD_C_SetWaveFieldData, MD_C_UpdateStates, MD_C_CalcOutput
   use NWTC_C_Binding, ONLY: IntfStrLen, SetErrStat_C, SetErrStat_F2C, ErrMsgLen_C, StringConvert_F2C, FileNameFromCString, AbortErrLev_C
   use WaveTank_Types
   use WaveTank_IO
   use WaveTank_Struct

   implicit none
   save

   public :: WaveTank_Init
   public :: WaveTank_CalcStep
   public :: WaveTank_End

   ! output to screen or to file (LabView doesn't capture console output nicely)
   integer(IntKi)    :: ScreenLogOutput_Un = -1
   character(1024)   :: ScreenLogOutput_File

   ! Simulation data storage
   type(SimSettingsType), target :: SimSettings

   ! IO data storage for CalcStep
   type(CalcStepIOdataType)      :: CalcStepIO

   ! Output file writing: headers, units, data, filename, fileunit etc.
   type(WrOutputDataType)        :: WrOutputData

   ! Structural model data storage
   type(MeshesMotionType)  :: MeshMotions    ! motion meshes (inputs)
   type(MeshesLoadsType )  :: MeshLoads      ! load meshes   (output)
   type(MeshesMapsType  )  :: MeshMaps       ! mappings
   type(StructTmpType   )  :: StructTmp      ! temporary data - avoids reallocation

   ! time stuff
   integer(IntKi)          :: VTKn_Global    ! global timestep for VTK
   integer(IntKi)          :: VTKn_last      ! last global timestep for VTK

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

   ! temporary storage of output channels
   real(c_float), allocatable :: OutData_SS_c(:)
   real(c_float), allocatable :: OutData_MD_c(:)
   real(c_float), allocatable :: OutData_ADI_c(:)

!TODO:
!     - add echo file
!     - add summary file
!     - add scaling
!        - Input for scaling already in place
!        - add info into summary file on scaling
!        - add unscaled interface IO outputs to file as well as the regular IO currently in there
!        - add pre and post scaling routines for time, pos, vel, acc, force/moment


contains

subroutine WaveTank_Init(  &
   WT_InputFile_C,         &
   RootName_C,             &
   VTKdir_C,               &
   ErrStat_C,              &
   ErrMsg_C                &
) bind (C, name='WaveTank_Init')
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

   ! The length of these arrays much match what is set in the corresponding C binding modules, or be larger
   character(kind=c_char) :: SS_WriteOutputHdr_C(ChanLen*MaxOutPts+1)
   character(kind=c_char) :: SS_WriteOutputUnt_C(ChanLen*MaxOutPts+1)
   character(kind=c_char) :: MD_WriteOutputHdr_C(ChanLen*1000)               ! probably oversized
   character(kind=c_char) :: MD_WriteOutputUnt_C(ChanLen*1000)               ! probably oversized
   character(kind=c_char) :: ADI_WriteOutputHdr_C(ChanLen*MaxADIOutputs+1)
   character(kind=c_char) :: ADI_WriteOutputUnt_C(ChanLen*MaxADIOutputs+1)

   ! Filename conversions -- read in as fortran strings, but sent to other modules as c_char arrays
   character(kind=c_char)         :: SS_InputFile_C(IntfStrLen)
   character(kind=c_char), target :: MD_InputFile_C(IntfStrLen)
   character(kind=c_char), target :: AD_InputFile_C(IntfStrLen)
   character(kind=c_char), target :: IfW_InputFile_C(IntfStrLen)

   ! temporary storage of number of output channels
   integer(c_int)    :: SS_NumChannels_C
   integer(c_int)    :: MD_NumChannels_C
   integer(c_int)    :: ADI_NumChannels_C

   ! Initialize error handling
   ErrStat_C = ErrID_None
   ErrMsg_C  = " "//C_NULL_CHAR

   InputFile = transfer(WT_InputFile_C, InputFile)
   i = index(InputFile, char(0))
   InputFile = InputFile(1:i)
   call ProcessComFile(InputFile,  FileInfo_In, ErrStat_F2, ErrMsg_F2); if (Failed()) return
   call ParseInputFile(FileInfo_In, SimSettings, ErrStat_F2, ErrMsg_F2); if (Failed()) return

   call ValidateInputFile(SimSettings, ErrStat_F2, ErrMsg_F2); if (Failed()) return

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
   if (SimSettings%Sim%DebugLevel > 2_c_int) call Print_FileInfo_Struct(CU,FileInfo_In)

   ! VTK directory
   WrVTK_Dir_C = c_null_char
   WrVTK_Dir_C = transfer( trim(SimSettings%Viz%WrVTK_Dir), WrVTK_Dir_C )
   ! return VTKdir
   VTKdir_C = c_null_char
   if (SimSettings%Viz%WrVTK > 0_c_int) VTKdir_C = WrVTK_Dir_C


   !------------------------------
   ! Build struct model
   !------------------------------
   call StructCreate(SimSettings, MeshMotions, MeshLoads, MeshMaps, StructTmp, ErrStat_F2, ErrMsg_F2)
   if (Failed()) return

   ! output VTK for struct model (if requested)
   if (SimSettings%Viz%WrVTK > 0_c_int) then
      ! create directory if doesn't exist
      call MKDIR( trim(SimSettings%Viz%WrVTK_Dir) )
      ! write mesh refs
      call WrVTK_Struct_Ref(SimSettings, MeshMotions, MeshLoads, ErrStat_F2, ErrMsg_F2)
      if (Failed()) return
   endif


   !------------------------------
   ! Setup and initialize SeaState
   !------------------------------
   call SeaSt_C_PreInit(            &
      SimSettings%Env%Gravity,      &
      SimSettings%Env%WtrDens,      &
      SimSettings%Env%WtrDpth,      &
      SimSettings%Env%MSL2SWL,      &
      SimSettings%Sim%DebugLevel-1, &     ! adjust down 1 for SS
      WrVTK_Dir_C,                  &
      SimSettings%Viz%WrVTK,        &
      SimSettings%Viz%WrVTK_DT,     &
      ErrStat_C2, ErrMsg_C2         &
   )
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'SeaSt_C_PreInit')
   if (ErrStat_C >= AbortErrLev_C) then
      call CleanUp()
      return
   endif

   SS_InputFile_C = c_null_char
   SS_InputFile_C = transfer(trim(SimSettings%ModSettings%SS_InputFile ), SS_InputFile_C )
   OutRootName_C  = transfer(trim(SimSettings%Sim%OutRootName)//'.SeaSt'//c_null_char, OutRootName_C)
   call SeaSt_C_Init(            &
      SS_InputFile_C,            &
      OutRootName_C,             &
      1000_c_int,                & !FIXME: do we need the number of timesteps???
      SimSettings%Sim%DT,        &
      SimSettings%ModSettings%WaveTimeShift, & 
      SS_NumChannels_C,          &
      SS_WriteOutputHdr_C,       &
      SS_WriteOutputUnt_C,       &
      ErrStat_C2, ErrMsg_C2      &
   )
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'SeaSt_C_Init')
   if (ErrStat_C >= AbortErrLev_C) then
      call CleanUp()
      return
   endif

   ! store channel info
   WrOutputData%NumChans_SS = int(SS_NumChannels_c,IntKi)
   call TransferOutChanNamesUnits(WrOutputData%NumChans_SS, 'WriteOutputHdr_SS', SS_WriteOutputHdr_c, WrOutputData%WriteOutputHdr_SS,ErrStat_F2,ErrMsg_F2); if (Failed()) return
   call TransferOutChanNamesUnits(WrOutputData%NumChans_SS, 'WriteOutputUnt_SS', SS_WriteOutputUnt_c, WrOutputData%WriteOutputUnt_SS,ErrStat_F2,ErrMsg_F2); if (Failed()) return


   !------------------------------
   ! Set the SeaState Wave Field pointer onto MoorDyn
   !------------------------------
   call WaveTank_SetWaveFieldPointer(ErrStat_C2, ErrMsg_C2)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_SetWaveFieldPointer')
   if (ErrStat_C >= AbortErrLev_C) then
      call CleanUp()
      return
   endif


   !------------------------------
   ! Setup and initialize MoorDyn
   !------------------------------
!FIXME: this interface will change!!! -- Split with PreInit
!FIXME: add WrVTK_Dir_C,  SimSettings%Viz%WrVTK, SimSettings%Viz%WrVTK_DT
   !SimSettings%TrbCfg%PtfmRef
!FIXME: 6 DOF with 3xPos, 3xEulerAngle
   InitPtfmPosOri = 0.0_c_float
   MD_InputFile_C = c_null_char
   MD_InputFile_C = transfer(trim(SimSettings%ModSettings%MD_InputFile ), MD_InputFile_C )
   OutRootName_C  = transfer(trim(SimSettings%Sim%OutRootName)//'.MD'//c_null_char, OutRootName_C)
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
      MD_WriteOutputHdr_C,                   &
      MD_WriteOutputUnt_C,                   &
      ErrStat_C2, ErrMsg_C2                  &
   )
!FIXME: add this when updating MD interface
!      SimSettings%Sim%DebugLevel-1, &     ! adjust down 1 for MD
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'MD_C_Init')
   if (ErrStat_C >= AbortErrLev_C) then
      call CleanUp()
      return
   endif

   ! store channel info
   WrOutputData%NumChans_MD = int(MD_NumChannels_c,IntKi)
   call TransferOutChanNamesUnits(WrOutputData%NumChans_MD, 'WriteOutputHdr_MD', MD_WriteOutputHdr_c, WrOutputData%WriteOutputHdr_MD,ErrStat_F2,ErrMsg_F2); if (Failed()) return
   call TransferOutChanNamesUnits(WrOutputData%NumChans_MD, 'WriteOutputUnt_MD', MD_WriteOutputUnt_c, WrOutputData%WriteOutputUnt_MD,ErrStat_F2,ErrMsg_F2); if (Failed()) return

   !------------------------------
   ! Setup and initialize AeroDyn+Inflow
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
      SimSettings%Sim%DebugLevel-1,          &     ! adjust down 1 for ADI
      ErrStat_C2, ErrMsg_C2                  &
   )
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_PreInit')
   if (ErrStat_C >= AbortErrLev_C) then
      call CleanUp()
      return
   endif

!FIXME: temporary location info until meshes set up
tmpNacPos = SimSettings%TrbCfg%TowerHt
tmpHubPos = tmpNacPos + SimSettings%TrbCfg%OverHang + SimSettings%TrbCfg%Twr2Shft     ! missing angles
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
      int(SimSettings%TrbCfg%NumBl,c_int),   &  ! NumBlades
      tmpBldRootPos      ,     &  ! BldRootPos
      tmpBldRootOri      ,     &  ! BldRootOri
      tmpNumMeshPts      ,     &  ! NumMeshPts
      tmpInitMeshPos     ,     &  ! InitMeshPos
      tmpInitMeshOri     ,     &  ! InitMeshOri
      tmpMeshPtToBladeNum,     &  ! MeshPtToBladeNum
      ErrStat_C2, ErrMsg_C2                  &
   )
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_SetupRotor')
   if (ErrStat_C >= AbortErrLev_C) then
      call CleanUp()
      return
   endif

   AD_InputFile_C  = c_null_char
   AD_InputFile_C  = transfer(trim(SimSettings%ModSettings%AD_InputFile ), AD_InputFile_C )
   IfW_InputFile_C = c_null_char
   IfW_InputFile_C = transfer(trim(SimSettings%ModSettings%IfW_InputFile), IfW_InputFile_C)
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
      SimSettings%TrbCfg%HubRad,             &  ! Hub radius for VTK surface rendering
      1_c_int,                               &  ! wrOuts_C -- Write ADI output file -- hard code to true for now
      SimSettings%Sim%DT,                    &  ! Timestep to write output file from ADI
      ADI_NumChannels_C,                     &
      ADI_WriteOutputHdr_C,                  &
      ADI_WriteOutputUnt_C,                  &
      ErrStat_C2, ErrMsg_C2                  &
   )
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_Init')
   if (ErrStat_C >= AbortErrLev_C) then
      call CleanUp()
      return
   endif

   ! store channel info
   WrOutputData%NumChans_ADI = int(ADI_NumChannels_c,IntKi)
   call TransferOutChanNamesUnits(WrOutputData%NumChans_ADI, 'WriteOutputHdr_ADI', ADI_WriteOutputHdr_c, WrOutputData%WriteOutputHdr_ADI, ErrStat_F2, ErrMsg_F2); if (Failed()) return
   call TransferOutChanNamesUnits(WrOutputData%NumChans_ADI, 'WriteOutputUnt_ADI', ADI_WriteOutputUnt_c, WrOutputData%WriteOutputUnt_ADI, ErrStat_F2, ErrMsg_F2); if (Failed()) return


   !------------------------------
   ! Assemble data for output file
   !------------------------------
   if (SimSettings%Outs%OutFile > 0_IntKi) then
      WrOutputData%OutName = trim(SimSettings%Sim%OutRootName)//'.out'
      call InitOutputFile(WrOutputData,ErrStat_F2,ErrMsg_F2); if (Failed()) return

      ! allocate storage for output channels from each of the modules, and c_float versions
      call AllocAry(WrOutputData%OutData_SS,    WrOutputData%Numchans_SS,  'OutData_SS',    ErrStat_F2, ErrMsg_F2); if (Failed()) return
      call AllocAry(WrOutputData%OutData_SS_c,  WrOutputData%Numchans_SS,  'OutData_SS_c',  ErrStat_F2, ErrMsg_F2); if (Failed()) return
      call AllocAry(WrOutputData%OutData_MD,    WrOutputData%Numchans_MD,  'OutData_MD',    ErrStat_F2, ErrMsg_F2); if (Failed()) return
      call AllocAry(WrOutputData%OutData_MD_c,  WrOutputData%Numchans_MD,  'OutData_MD_c',  ErrStat_F2, ErrMsg_F2); if (Failed()) return
      call AllocAry(WrOutputData%OutData_ADI,   WrOutputData%Numchans_ADI, 'OutData_ADI',   ErrStat_F2, ErrMsg_F2); if (Failed()) return
      call AllocAry(WrOutputData%OutData_ADI_c, WrOutputData%Numchans_ADI, 'OutData_ADI_c', ErrStat_F2, ErrMsg_F2); if (Failed()) return
   endif

   !------------------------------
   ! Final cleanup
   !------------------------------
   ! Initialize time counting for VTK
   VTKn_Global = 0_IntKi
   call ShowReturnData()

contains
   logical function Failed()
      call SetErrStat_F2C(ErrStat_F2, ErrMsg_F2, ErrStat_C, ErrMsg_C)
      Failed = ErrStat_C >= AbortErrLev_C
      if (Failed)    call Cleanup()
   end function Failed
   subroutine Cleanup()
      call NWTC_Library_DestroyFileInfoType(FileInfo_In, ErrStat_F2, ErrMsg_F2)  ! ignore error from this
      if (ScreenLogOutput_Un > 0)   close(ScreenLogOutput_Un)
   end subroutine Cleanup
   subroutine ShowPassedData()
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  WaveTank_Init input values")
      call WrScr("   --------------------------------------------------------")
      call WrScr("   WT_InputFile_C         -> "//trim(InputFile))
      call WrScr("   --------------------------------------------------------")
   end subroutine ShowPassedData
   subroutine ShowReturnData()
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  WaveTank_Init returned values")
      call WrScr("   --------------------------------------------------------")
      call WrScr("   RootName_C             <- "//trim(SimSettings%Sim%OutRootName))
      call WrScr("   WrVTK_Dir_C            <- "//trim(SimSettings%Viz%WrVTK_Dir))
      call WrScr("-----------------------------------------------------------")
   end subroutine
end subroutine WaveTank_Init

!subroutine DeallocEverything()
!end subroutine DeallocEverything



subroutine WaveTank_CalcStep( &
   time_c,                    &
   pos_c,                     &
   vel_c,                     &
   acc_c,                     &
   loads_c,                   &
   buoyWaveElev_c,            &
   ErrStat_C,                 &
   ErrMsg_C                   &
) BIND (C, NAME='WaveTank_CalcStep')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: WaveTank_CalcStep
!GCC$ ATTRIBUTES DLLEXPORT :: WaveTank_CalcStep
#endif
   real(c_double),         intent(in   ) :: time_c
   real(c_float),          intent(in   ) :: pos_c(6)        ! [x,y,z,roll,pitch,yaw]
   real(c_float),          intent(in   ) :: vel_c(6)        ! [x_dot,y_dot,z_dot,roll_dot,pitch_dot,yaw_dot]
   real(c_float),          intent(in   ) :: acc_c(6)        ! [x_ddot,y_ddot,z_ddot,roll_ddot,pitch_ddot,yaw_ddot]
   real(c_float),          intent(  out) :: loads_c(6)      ! [Fx,Fy,Fz,Mx,My,Mz]
   real(c_float),          intent(  out) :: buoyWaveElev_c  ! wave elevation at buoy
   integer(c_int),         intent(  out) :: ErrStat_C
   character(c_char),      intent(  out) :: ErrMsg_C(ErrMsgLen_C)
   integer(c_int)                        :: ErrStat_C2
   character(c_char)                     :: ErrMsg_C2(ErrMsgLen_C)
   integer(IntKi)                        :: ErrStat_F2
   character(ErrMsgLen)                  :: ErrMsg_F2
   real(c_float)                         :: tmpPos_C(2)     ! temporary for wave buoy position

   ! Initialize error handling
   ErrStat_C = ErrID_None
   ErrMsg_C  = " "//C_NULL_CHAR


   ! debugging
   if (SimSettings%Sim%DebugLevel > 0_c_int) call ShowPassedData()

   ! zero loads in case of error
   loads_c = 0.0_c_float

   ! Transfer CalcStepIO data (storing for output to file)
   CalcStepIO%Time_c    = time_C
   CalcStepIO%PosAng_c  = pos_c
   CalcStepIO%Vel_c     = vel_c
   CalcStepIO%Acc_c     = acc_c


   !--------------------------------------
   ! Update motion meshes
   !--------------------------------------
!FIXME: add this


   !--------------------------------------
   ! Wave elevation at buoy, update buoy
   !--------------------------------------
   tmpPos_C = real(SimSettings%WaveBuoy%XYLoc, c_float)
   call SeaSt_C_GetSurfElev(Time_C, tmpPos_C, buoyWaveElev_c, ErrStat_C2, ErrMsg_C2)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_CalcStep::SeaSt_C_GetSurfElev')
   if (ErrStat_C >= AbortErrLev_C) return
   MeshMotions%WaveBuoyMotion%TranslationDisp(:,1) = (/ 0.0_ReKi, 0.0_ReKi, real(buoyWaveElev_c, ReKi) /)

 
   !--------------------------------------
   ! Write VTK if requested
   !     Do this here in case failed calcs
   !--------------------------------------
   if (SimSettings%Viz%WrVTK > 0_c_int) then
      ! only write on desired time interval (same logic used in c-binding modules)
      VTKn_Global = nint(Time_C / SimSettings%Viz%WrVTK_DT)
      if (VTKn_Global /= VTKn_last) then   ! already wrote this one
         VTKn_last = VTKn_Global           ! store the current number to make sure we don't write it twice
         call WrVTK_Struct(VTKn_Global, SimSettings, MeshMotions, MeshLoads, ErrStat_F2, ErrMsg_F2)
         if (Failed()) return
      endif
   endif

  
   !--------------------------------------
   ! call SeaState_Calc (writes vis)
   !--------------------------------------
   call SeaSt_C_CalcOutput(Time_C, WrOutputData%OutData_SS_c, ErrStat_C, ErrMsg_C)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_CalcStep::SeaSt_C_CalcOutput')
   if (ErrStat_C >= AbortErrLev_C) return
   ! transfer data for writing out
   WrOutputData%OutData_SS = real(WrOutputData%OutData_SS_c, ReKi)

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


   ! debugging
   if (SimSettings%Sim%DebugLevel > 0_c_int) call ShowReturnData()

   ! Store returned values
   CalcStepIO%FrcMom_C = 0.0_c_float
   loads_c = CalcStepIO%FrcMom_C

   ! Transfer outputs and write to file
   if (SimSettings%Outs%OutFile > 0_IntKi) then
      ! store mapped loads from modules at the platform point
!FIXME
!      CalcStepIO%FrcMom_SS
!      CalcStepIO%FrcMom_MD
!      CalcStepIO%FrcMom_ADI

      ! store outputs from modules
!FIXME
      !WrOutputData%OutData_SS,
      !WrOutputData%OutData_MD,
      !WrOutputData%OutData_ADI
      ! output to file
      call WriteOutputLine(SimSettings%Outs%OutFmt, CalcStepIO, StructTmp, WrOutputData, ErrStat_F2, ErrMsg_F2)
      if (Failed()) return
   endif
contains
   logical function Failed()
      call SetErrStat_F2C(ErrStat_F2, ErrMsg_F2, ErrStat_C, ErrMsg_C)
      Failed = ErrStat_C >= AbortErrLev_C
      !if (Failed)    call Cleanup()
   end function Failed
   subroutine ShowPassedData()
      character(120) :: TmpStr
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  WaveTank_CalcStep input values")
      call WrScr("   --------------------------------------------------------")
      call WrScr("   time_c    -> "//trim(Num2LStr(time_c)))
      write(TmpStr,'("(", *(f10.5, :, ","))') pos_c;  TmpStr=trim(TmpStr)//" )"
      call WrScr("   pos_c     -> "//trim(TmpStr))
      write(TmpStr,'("(", *(f10.5, :, ","))') vel_c;  TmpStr=trim(TmpStr)//" )"
      call WrScr("   vel_c     -> "//trim(TmpStr))
      write(TmpStr,'("(", *(f10.5, :, ","))') acc_c;  TmpStr=trim(TmpStr)//" )"
      call WrScr("   acc_c     -> "//trim(TmpStr))
      call WrScr("   --------------------------------------------------------")
   end subroutine ShowPassedData
   subroutine ShowReturnData()
      character(120) :: TmpStr
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  WaveTank_CalcStep returned values")
      call WrScr("   --------------------------------------------------------")
      write(TmpStr,'("(", *(f10.5, :, ","))') loads_c;  TmpStr=trim(TmpStr)//" )"
      call WrScr("   loads_c        <- "//trim(TmpStr))
      call WrScr("   buoyWaveElev_c <- "//trim(Num2LStr(buoyWaveElev_c)))
      call WrScr("-----------------------------------------------------------")
   end subroutine
end subroutine


subroutine WaveTank_End(ErrStat_C, ErrMsg_C) bind (C, NAME="WaveTank_End")
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

   ! destroy mesh info
   call StructDestroy(MeshMotions, MeshLoads, MeshMaps, StructTmp, ErrStat_F2, ErrMsg_F2)
   call SetErrStat_F2C(ErrStat_F2, ErrMsg_F2, ErrStat_C, ErrMsg_C)

   ! in case we were writing to a file instead of the screen
   if (ScreenLogOutput_Un > 0)   close(ScreenLogOutput_Un)

   call MD_C_END(ErrStat_C2, ErrMsg_C2)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'MD_C_END')

   call SeaSt_C_END(ErrStat_C2, ErrMsg_C2)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'SeaSt_C_END')

   call ADI_C_END(ErrStat_C2, ErrMsg_C2)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_END')

   ! close output file
   if (WrOutputData%OutUn > 0) then
      close(WrOutputData%OutUn, iostat=ErrStat_F2)
      if (ErrStat_F2 /= 0_IntKi) call SetErrStat_C(int(ErrID_Fatal,c_int), 'could no close output file '//trim(WrOutputData%OutName), ErrStat_C, ErrMsg_C, 'ADI_C_END')
      WrOutputData%OutUn = -1    ! mark as closed - prevents faults
   endif

end subroutine


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

   ! There isn't a good way to check for an error here.  Will get caught at init
   call MD_C_SetWaveFieldData(WaveFieldPointer_C)

   ! Probably doesn't matter, but clear the fortran pointer just in case
   WaveFieldPointer_F => NULL()
end subroutine

END MODULE WaveTankTesting
