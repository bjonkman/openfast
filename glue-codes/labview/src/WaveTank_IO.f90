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

subroutine ParseInputFile(FileInfo_In, InitInp, ErrStat, ErrMsg)
   type(FileInfoType),         intent(in   )  :: FileInfo_In       !< The derived type for holding the file information.
   type(WaveTank_InitInput),   intent(  out)  :: InitInp
   integer(IntKi),             intent(  out)  :: ErrStat
   character(*),               intent(  out)  :: ErrMsg

   ! Local variables
   integer                                    :: CurLine
   character(1024), target                    :: TmpPath
   character(1024)                            :: FileName

   integer(IntKi)                             :: ErrStat2             ! local status of error message
   character(ErrMsgLen)                       :: ErrMsg2              ! local error message if errStat /= ErrID_None

   character(*), parameter                    :: RoutineName = 'WaveTankTesting.ParseInputFile'

   ErrStat = ErrID_None
   ErrMsg  = " "

   CurLine = 1
   ! Separator line skipped
   call ParseVar( FileInfo_In, CurLine, 'DT', InitInp%DT, ErrStat2, ErrMsg2); 

   call ParseVar( FileInfo_In, CurLine, 'SS_OutRootName_C', TmpPath, ErrStat2, ErrMsg2); if(Failed()) return;
   InitInp%SS_OutRootName_C = transfer(TmpPath, InitInp%SS_OutRootName_C)
   call ParseVar( FileInfo_In, CurLine, 'SS_Gravity_C', InitInp%SS_Gravity_C, ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'SS_WtrDens_C', InitInp%SS_WtrDens_C, ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'SS_WtrDpth_C', InitInp%SS_WtrDpth_C, ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'SS_MSL2SWL_C', InitInp%SS_MSL2SWL_C, ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'SS_NSteps_C',  InitInp%SS_NSteps_C,  ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'SS_TimeInterval_C',       InitInp%SS_TimeInterval_C,       ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'SS_WaveElevSeriesFlag_C', InitInp%SS_WaveElevSeriesFlag_C, ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'SS_WrWvKinMod_C',         InitInp%SS_WrWvKinMod_C,         ErrStat2, ErrMsg2); if(Failed()) return;

   ! Separator line skipped
   call ParseVar( FileInfo_In, CurLine, 'MD_G_C',            InitInp%MD_G_C,              ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'MD_RHO_C',          InitInp%MD_RHO_C,            ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'MD_DEPTH_C',        InitInp%MD_DEPTH_C,          ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseAry( FileInfo_In, CurLine, 'MD_PtfmInit_C',     InitInp%MD_PtfmInit_C,    6, ErrStat2,  ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'MD_InterpOrder_C',  InitInp%MD_InterpOrder_C,    ErrStat2, ErrMsg2); if(Failed()) return;

   ! Separator line skipped
   call ParseVar( FileInfo_In, CurLine, 'NumTurbines_C',     InitInp%NumTurbines_C,     ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'TransposeDCM',      InitInp%TransposeDCM,      ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'PointLoadOutput',   InitInp%PointLoadOutput,   ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'ADI_gravity_C',     InitInp%ADI_gravity_C,     ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'ADI_defFldDens_C',  InitInp%ADI_defFldDens_C,  ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'ADI_defKinVisc_C',  InitInp%ADI_defKinVisc_C,  ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'ADI_defSpdSound_C', InitInp%ADI_defSpdSound_C, ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'ADI_defPatm_C',     InitInp%ADI_defPatm_C,     ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'ADI_defPvap_C',     InitInp%ADI_defPvap_C,     ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'ADI_WtrDpth_C',     InitInp%ADI_WtrDpth_C,     ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'ADI_MSL2SWL_C',     InitInp%ADI_MSL2SWL_C,     ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'MHK',               InitInp%MHK,               ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'DebugLevel',        InitInp%DebugLevel,        ErrStat2, ErrMsg2); if(Failed()) return;

   call ParseVar( FileInfo_In, CurLine, 'iWT_c',            InitInp%iWT_c,              ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'TurbineIsHAWT_c',  InitInp%TurbineIsHAWT_c,    ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseAry( FileInfo_In, CurLine, 'TurbOrigin_C',     InitInp%TurbOrigin_C,    3, ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseAry( FileInfo_In, CurLine, 'HubPos_C',         InitInp%HubPos_C,        3, ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseAry( FileInfo_In, CurLine, 'HubOri_C',         InitInp%HubOri_C,        9, ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseAry( FileInfo_In, CurLine, 'NacPos_C',         InitInp%NacPos_C,        3, ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseAry( FileInfo_In, CurLine, 'NacOri_C',         InitInp%NacOri_C,        9, ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'NumBlades_C',      InitInp%NumBlades_C,        ErrStat2, ErrMsg2); if(Failed()) return;

   call AllocAry(InitInp%BldRootPos_C, 3*InitInp%NumBlades_C, 'BldRootPos_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call AllocAry(InitInp%BldRootOri_C, 9*InitInp%NumBlades_C, 'BldRootPos_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseAry( FileInfo_In, CurLine, 'BldRootPos_C', InitInp%BldRootPos_C, 3*InitInp%NumBlades_C, ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseAry( FileInfo_In, CurLine, 'BldRootOri_C', InitInp%BldRootOri_C, 9*InitInp%NumBlades_C, ErrStat2, ErrMsg2); if(Failed()) return;

   call ParseVar( FileInfo_In, CurLine, 'NumMeshPts_C',       InitInp%NumMeshPts_C,                                ErrStat2, ErrMsg2); if(Failed()) return;
   call AllocAry(InitInp%InitMeshPos_C, 3*InitInp%NumMeshPts_C, 'InitMeshPos_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call AllocAry(InitInp%InitMeshOri_C, 9*InitInp%NumMeshPts_C, 'InitMeshOri_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call AllocAry(InitInp%MeshPtToBladeNum_C, InitInp%NumMeshPts_C, 'MeshPtToBladeNum_C', ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseAry( FileInfo_In, CurLine, 'InitMeshPos_C',      InitInp%InitMeshPos_C,       3*InitInp%NumMeshPts_C, ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseAry( FileInfo_In, CurLine, 'InitMeshOri_C',      InitInp%InitMeshOri_C,       9*InitInp%NumMeshPts_C, ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseAry( FileInfo_In, CurLine, 'MeshPtToBladeNum_C', InitInp%MeshPtToBladeNum_C,  InitInp%NumMeshPts_C,   ErrStat2, ErrMsg2); if(Failed()) return;

   call ParseVar( FileInfo_In, CurLine, 'ADI_OutRootName_C', TmpPath, ErrStat2, ErrMsg2); if(Failed()) return;
   call StringConvert_F2C(TmpPath, InitInp%ADI_OutRootName_C)
   call ParseVar( FileInfo_In, CurLine, 'ADI_OutVTKDir_C',   TmpPath, ErrStat2, ErrMsg2); if(Failed()) return;
   call StringConvert_F2C(TmpPath, InitInp%ADI_OutVTKDir_C)
   call ParseVar( FileInfo_In, CurLine, 'ADI_InterpOrder_C', InitInp%ADI_InterpOrder_C,   ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'ADI_TMax_C',        InitInp%ADI_TMax_C,          ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'ADI_storeHHVel',    InitInp%ADI_storeHHVel,      ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'ADI_WrVTK_in',      InitInp%ADI_WrVTK_in,        ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'ADI_WrVTK_inType',  InitInp%ADI_WrVTK_inType,    ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'ADI_WrVTK_inDT',    InitInp%ADI_WrVTK_inDT,      ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseAry( FileInfo_In, CurLine, 'ADI_VTKNacDim_in',  InitInp%ADI_VTKNacDim_in, 6, ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'ADI_VTKHubrad_in',  InitInp%ADI_VTKHubrad_in,    ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'ADI_wrOuts_C',      InitInp%ADI_wrOuts_C,        ErrStat2, ErrMsg2); if(Failed()) return;
   call ParseVar( FileInfo_In, CurLine, 'ADI_DT_Outs_C',     InitInp%ADI_DT_Outs_C,       ErrStat2, ErrMsg2); if(Failed()) return;

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine


END MODULE WaveTank_IO
