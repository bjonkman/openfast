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

subroutine ReadInput(InputFilePath, InitInp, ErrStat, ErrMsg)

   character(*),               intent(in   )   :: InputFilePath
   type(WaveTank_InitInput),   intent(  out)   :: InitInp
   integer(IntKi),             intent(  out)   :: ErrStat
   character(*),               intent(  out)   :: ErrMsg

   ! Local variables
   integer :: UnIn = -1
   character(1024), target                     :: TmpPath
   character(1024)                             :: FileName

   integer(IntKi)                              :: ErrStat2             ! local status of error message
   character(ErrMsgLen)                        :: ErrMsg2              ! local error message if errStat /= ErrID_None

   character(*), parameter                     :: RoutineName = 'WaveTankTesting.ReadInput'

   ErrStat = ErrID_None
   ErrMsg  = " "

   FileName = TRIM(InputFilePath)
   call GetNewUnit( UnIn )
   call OpenFInpFile( UnIn, FileName, ErrStat2, ErrMsg2); if(Failed()) return;

   call ReadCom( UnIn, FileName, 'Init comment', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%DT, 'DT', 'DT', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadCom( UnIn, FileName, 'SeaState Init comment', ErrStat2, ErrMsg2); if(Failed()) return;
   if (ErrStat >= AbortErrLev) return

   call ReadVar( UnIn, FileName, TmpPath, 'SS_OutRootName_C', 'SS_OutRootName_C', ErrStat2, ErrMsg2); if(Failed()) return;
   InitInp%SS_OutRootName_C = transfer(TmpPath, InitInp%SS_OutRootName_C)
   call ReadVar( UnIn, FileName, InitInp%SS_Gravity_C, 'SS_Gravity_C', 'SS_Gravity_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%SS_WtrDens_C, 'SS_WtrDens_C', 'SS_WtrDens_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%SS_WtrDpth_C, 'SS_WtrDpth_C', 'SS_WtrDpth_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%SS_MSL2SWL_C, 'SS_MSL2SWL_C', 'SS_MSL2SWL_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%SS_NSteps_C, 'SS_NSteps_C', 'SS_NSteps_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%SS_TimeInterval_C, 'SS_TimeInterval_C', 'SS_TimeInterval_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%SS_WaveElevSeriesFlag_C, 'SS_WaveElevSeriesFlag_C', 'SS_WaveElevSeriesFlag_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%SS_WrWvKinMod_C, 'SS_WrWvKinMod_C', 'SS_WrWvKinMod_C', ErrStat2, ErrMsg2); if(Failed()) return;
   if (ErrStat >= AbortErrLev) return

   call ReadCom( UnIn, FileName, 'MoorDyn Init comment', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%MD_G_C, 'MD_G_C', 'MD_G_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%MD_RHO_C, 'MD_RHO_C', 'MD_RHO_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%MD_DEPTH_C, 'MD_DEPTH_C', 'MD_DEPTH_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadAry( UnIn, FileName, InitInp%MD_PtfmInit_C, 6, 'MD_PtfmInit_C', 'MD_PtfmInit_C', ErrStat2,  ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%MD_InterpOrder_C, 'MD_InterpOrder_C', 'MD_InterpOrder_C', ErrStat2, ErrMsg2); if(Failed()) return;
   if (ErrStat >= AbortErrLev) return

   call ReadCom( UnIn, FileName, 'ADI Init comment', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%NumTurbines_C, 'NumTurbines_C', 'NumTurbines_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%TransposeDCM, 'TransposeDCM', 'TransposeDCM', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%PointLoadOutput, 'PointLoadOutput', 'PointLoadOutput', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%ADI_gravity_C, 'ADI_gravity_C', 'ADI_gravity_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%ADI_defFldDens_C, 'ADI_defFldDens_C', 'ADI_defFldDens_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%ADI_defKinVisc_C, 'ADI_defKinVisc_C', 'ADI_defKinVisc_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%ADI_defSpdSound_C, 'ADI_defSpdSound_C', 'ADI_defSpdSound_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%ADI_defPatm_C, 'ADI_defPatm_C', 'ADI_defPatm_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%ADI_defPvap_C, 'ADI_defPvap_C', 'ADI_defPvap_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%ADI_WtrDpth_C, 'ADI_WtrDpth_C', 'ADI_WtrDpth_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%ADI_MSL2SWL_C, 'ADI_MSL2SWL_C', 'ADI_MSL2SWL_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%MHK, 'MHK', 'MHK', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%DebugLevel, 'DebugLevel', 'DebugLevel', ErrStat2, ErrMsg2); if(Failed()) return;
   if (ErrStat >= AbortErrLev) return

   call ReadVar( UnIn, FileName, InitInp%iWT_c, 'iWT_c', 'iWT_c', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%TurbineIsHAWT_c, 'TurbineIsHAWT_c', 'TurbineIsHAWT_c', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadAry( UnIn, FileName, InitInp%TurbOrigin_C, 3, 'TurbOrigin_C', 'TurbOrigin_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadAry( UnIn, FileName, InitInp%HubPos_C, 3, 'HubPos_C', 'HubPos_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadAry( UnIn, FileName, InitInp%HubOri_C, 9, 'HubOri_C', 'HubOri_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadAry( UnIn, FileName, InitInp%NacPos_C, 3, 'NacPos_C', 'NacPos_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadAry( UnIn, FileName, InitInp%NacOri_C, 9, 'NacOri_C', 'NacOri_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%NumBlades_C, 'NumBlades_C', 'NumBlades_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call AllocAry(InitInp%BldRootPos_C, 3*InitInp%NumBlades_C, 'BldRootPos_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call AllocAry(InitInp%BldRootOri_C, 9*InitInp%NumBlades_C, 'BldRootPos_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadAry( UnIn, FileName, InitInp%BldRootPos_C, 3*InitInp%NumBlades_C, 'BldRootPos_C', 'BldRootPos_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadAry( UnIn, FileName, InitInp%BldRootOri_C, 9*InitInp%NumBlades_C, 'BldRootOri_C', 'BldRootOri_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%NumMeshPts_C, 'NumMeshPts_C', 'NumMeshPts_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call AllocAry(InitInp%InitMeshPos_C, 3*InitInp%NumMeshPts_C, 'InitMeshPos_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call AllocAry(InitInp%InitMeshOri_C, 9*InitInp%NumMeshPts_C, 'InitMeshOri_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call AllocAry(InitInp%MeshPtToBladeNum_C, InitInp%NumMeshPts_C, 'MeshPtToBladeNum_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadAry( UnIn, FileName, InitInp%InitMeshPos_C, 3*InitInp%NumMeshPts_C, 'InitMeshPos_C', 'InitMeshPos_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadAry( UnIn, FileName, InitInp%InitMeshOri_C, 9*InitInp%NumMeshPts_C, 'InitMeshOri_C', 'InitMeshOri_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadAry( UnIn, FileName, InitInp%MeshPtToBladeNum_C, InitInp%NumMeshPts_C, 'MeshPtToBladeNum_C', 'MeshPtToBladeNum_C', ErrStat2, ErrMsg2); if(Failed()) return;
   if (ErrStat >= AbortErrLev) return

   call ReadVar( UnIn, FileName, TmpPath, 'ADI_OutRootName_C', 'ADI_OutRootName_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call StringConvert_F2C(TmpPath, InitInp%ADI_OutRootName_C)
   call ReadVar( UnIn, FileName, TmpPath, 'ADI_OutVTKDir_C', 'ADI_OutVTKDir_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call StringConvert_F2C(TmpPath, InitInp%ADI_OutVTKDir_C)
   call ReadVar( UnIn, FileName, InitInp%ADI_InterpOrder_C, 'ADI_InterpOrder_C', 'ADI_InterpOrder_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%ADI_TMax_C, 'ADI_TMax_C', 'ADI_TMax_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%ADI_storeHHVel, 'ADI_storeHHVel', 'ADI_storeHHVel', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%ADI_WrVTK_in, 'ADI_WrVTK_in', 'ADI_WrVTK_in', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%ADI_WrVTK_inType, 'ADI_WrVTK_inType', 'ADI_WrVTK_inType', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%ADI_WrVTK_inDT, 'ADI_WrVTK_inDT', 'ADI_WrVTK_inDT', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadAry( UnIn, FileName, InitInp%ADI_VTKNacDim_in, 6, 'ADI_VTKNacDim_in', 'ADI_VTKNacDim_in', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%ADI_VTKHubrad_in, 'ADI_VTKHubrad_in', 'ADI_VTKHubrad_in', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%ADI_wrOuts_C, 'ADI_wrOuts_C', 'ADI_wrOuts_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ReadVar( UnIn, FileName, InitInp%ADI_DT_Outs_C, 'ADI_DT_Outs_C', 'ADI_DT_Outs_C', ErrStat2, ErrMsg2); if(Failed()) return;
   if (ErrStat >= AbortErrLev) return

   if(UnIn>0) close( UnIn )

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine


END MODULE WaveTank_IO
