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
MODULE WaveTank_Types

   use ISO_C_BINDING
   use NWTC_Library
   use NWTC_C_Binding, ONLY: IntfStrLen, SetErrStat_C, SetErrStat_F2C, ErrMsgLen_C, StringConvert_F2C, FileNameFromCString

   implicit none
   type WaveTank_InitInput
      real(c_double)      :: DT

      ! SeaState variables
      character(c_char)   :: SS_OutRootName_C(IntfStrLen)
      real(c_float)       :: SS_Gravity_C
      real(c_float)       :: SS_WtrDens_C
      real(c_float)       :: SS_WtrDpth_C
      real(c_float)       :: SS_MSL2SWL_C
      integer(c_int)      :: SS_NSteps_C
      real(c_float)       :: SS_TimeInterval_C
      integer(c_int)      :: SS_WaveElevSeriesFlag_C
      integer(c_int)      :: SS_WrWvKinMod_C

      ! MD variables
      ! real(c_double)    :: MD_DT_C                                !< Using global DT
      real(c_float)       :: MD_G_C
      real(c_float)       :: MD_RHO_C
      real(c_float)       :: MD_DEPTH_C
      real(c_float)       :: MD_PtfmInit_C(6)
      integer(c_int)      :: MD_InterpOrder_C

      ! ADI variables
      ! Preinit
      integer(c_int)  :: NumTurbines_C
      integer(c_int)  :: TransposeDCM
      integer(c_int)  :: PointLoadOutput
      real(c_float)   :: ADI_gravity_C
      real(c_float)   :: ADI_defFldDens_C
      real(c_float)   :: ADI_defKinVisc_C
      real(c_float)   :: ADI_defSpdSound_C
      real(c_float)   :: ADI_defPatm_C
      real(c_float)   :: ADI_defPvap_C
      real(c_float)   :: ADI_WtrDpth_C
      real(c_float)   :: ADI_MSL2SWL_C
      integer(c_int)  :: MHK
      integer(c_int)  :: DebugLevel
      ! SetupRotor
      integer(c_int)  :: iWT_c                                    !< Wind turbine / rotor number
      integer(c_int)  :: TurbineIsHAWT_c                          !< true for HAWT, false for VAWT
      real(c_float)   :: TurbOrigin_C(3)                          !< turbine origin (tower base). Gets added to all meshes to shift turbine position.
      real(c_float)   :: HubPos_C( 3 )                            !< Hub position
      real(c_double)  :: HubOri_C( 9 )                            !< Hub orientation
      real(c_float)   :: NacPos_C( 3 )                            !< Nacelle position
      real(c_double)  :: NacOri_C( 9 )                            !< Nacelle orientation
      integer(c_int)  :: NumBlades_C
      real(c_float), dimension(:), allocatable :: BldRootPos_C    !< Blade root positions; 3xNumBlades_C
      real(c_double), dimension(:), allocatable :: BldRootOri_C   !< Blade root orientations; 9xNumBlades_C
      ! Initial nodes
      integer(c_int)  :: NumMeshPts_C                             !< Number of mesh points we are transferring motions and outputting loads to
      real(c_float), dimension(:), allocatable :: InitMeshPos_C   !< A 3xNumMeshPts_C array [x,y,z]
      real(c_double), dimension(:), allocatable :: InitMeshOri_C  !< A 9xNumMeshPts_C array [r11,r12,r13,r21,r22,r23,r31,r32,r33]
      integer(c_int), dimension(:), allocatable :: MeshPtToBladeNum_C !< A NumMeshPts_C array of blade numbers associated with each mesh point
      ! Init
      character(kind=c_char) :: ADI_OutRootName_C(IntfStrLen)     !< Root name to use for echo files and other
      character(kind=c_char) :: ADI_OutVTKDir_C(IntfStrLen)       !< Directory to put all vtk output
      ! Interpolation
      integer(c_int) :: ADI_InterpOrder_C                         !< Interpolation order to use (must be 1 or 2)
      ! Time
      ! real(c_double) :: ADI_DT_C                                  !< Timestep used with AD for stepping forward from t to t+dt.  Must be constant.; Using global DT
      real(c_double) :: ADI_TMax_C                                !< Maximum time for simulation
      ! Flags
      integer(c_int) :: ADI_storeHHVel                            !< Store hub height time series from IfW
      ! VTK
      integer(c_int) :: ADI_WrVTK_in
      integer(c_int) :: ADI_WrVTK_inType
      real(c_double) :: ADI_WrVTK_inDT
      real(c_float)  :: ADI_VTKNacDim_in(6)
      real(c_float)  :: ADI_VTKHubrad_in
      integer(c_int) :: ADI_wrOuts_C
      real(c_double) :: ADI_DT_Outs_C
   end type WaveTank_InitInput
END MODULE WaveTank_Types
