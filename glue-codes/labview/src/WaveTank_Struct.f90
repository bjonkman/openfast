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
!  This module provides structural model for the wavetank interface
!
!**********************************************************************************************************************************
module WaveTank_Struct
   use ISO_C_BINDING
   use NWTC_Library
   use WaveTank_Types

   implicit none
   private

   save

   public :: StructCreate
   public :: StructDestroy
   public :: WrVTK_Struct_Ref
   public :: WrVTK_Struct

contains


!> create the structural model, allocate temp data storage, setup mesh mappings
subroutine StructCreate(SimSettings, MeshMotions, MeshLoads, MeshMaps, StructTmp, ErrStat, ErrMsg)
   type(SimSettingsType),  target,  intent(in   )  :: SimSettings
   type(MeshesMotionType),          intent(inout)  :: MeshMotions
   type(MeshesLoadsType ),          intent(inout)  :: MeshLoads
   type(MeshesMapsType  ),          intent(inout)  :: MeshMaps
   type(StructTmpType   ),          intent(inout)  :: StructTmp
   integer(IntKi),                  intent(  out)  :: ErrStat
   character(ErrMsgLen),            intent(  out)  :: ErrMsg
   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2
   character(*),           parameter               :: RoutineName = 'WaveTank::StructCreate'
   real(ReKi)                                      :: TmpPos(3)
   real(DbKi)                                      :: Orient(3,3) ! temporary orientation
   type(TurbConfigType),   pointer                 :: TrbCfg      ! to shorten notation
   type(TurbInitCondType), pointer                 :: TrbInit     ! to shorten notation
   ErrStat = ErrID_None
   ErrMsg  = ''

   TrbCfg    => SimSettings%TrbCfg
   TrbInit   => SimSettings%TrbInit

   ! Set some state information
   StructTmp%RotSpeed = TrbInit%RotSpeed
   StructTmp%BldPitch = TrbInit%BldPitch
   StructTmp%NacYaw   = TrbInit%NacYaw
   StructTmp%Azimuth  = TrbInit%Azimuth


   !-------------------------------
   ! Wave measurement buoy
   TmpPos = 0.0_ReKi
   TmpPos(1:2) = SimSettings%WaveBuoy%XYLoc(1:2)
   call Eye(Orient, ErrStat2, ErrMsg2);   if (Failed()) return
   call CreateInputPointMesh(MeshMotions%WaveBuoyMotion, TmpPos, Orient, ErrStat2, ErrMsg2, hasMotion=.true., hasLoads=.false.); if (Failed()) return


   !-------------------------------
   ! create PRP platform mesh point
   TmpPos = TrbCfg%PtfmRefPos
   Orient=WT_EulerToDCM_fromInput(TrbCfg%PtfmRefOrient)
   call CreateInputPointMesh(MeshMotions%PtfmPtMotion, TmpPos, Orient, ErrStat2, ErrMsg2, hasMotion=.true., hasLoads=.false.); if (Failed()) return

   !-------------------------------
   ! create 2 point tower mesh
   call MeshCreate ( BlankMesh = MeshMotions%TowerMotion, IOS=COMPONENT_INPUT, Nnodes=2, ErrStat=ErrStat2, ErrMess=ErrMsg2,  &
               Orientation = .true., TranslationDisp = .true., TranslationVel = .true., TranslationAcc  = .TRUE. )
   if (Failed()) return

   ! Tower bottom
   TmpPos = real(TrbCfg%TowerBsPt,ReKi)   ! c_float to ReKi
   call MeshPositionNode(MeshMotions%TowerMotion, 1, TmpPos, errStat2, errMsg2)  ! orientation is identity by default
   if (Failed()) return

   ! Tower top -- assumes vertical tower
   TmpPos(3) = real(TrbCfg%TowerHt,ReKi)   ! c_float to ReKi
   call MeshPositionNode(MeshMotions%TowerMotion, 2, TmpPos, errStat2, errMsg2)  ! orientation is identity by default
   if (Failed()) return

   ! create line element
   call MeshConstructElement( MeshMotions%TowerMotion, ELEMENT_LINE2, errStat2, errMsg2, p1=1, p2=2 )
   if (Failed()) return
  
   ! commit mesh          
   call MeshCommit(MeshMotions%TowerMotion, errStat2, errMsg2 )

   ! initialize location
   MeshMotions%TowerMotion%Orientation     = MeshMotions%TowerMotion%RefOrientation
   MeshMotions%TowerMotion%TranslationDisp = 0.0_R8Ki
   MeshMotions%TowerMotion%TranslationVel  = 0.0_ReKi
 

   !-------------------------------



 
contains
   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine


!> destroy all structural model related info
subroutine StructDestroy(MeshMotions, MeshLoads, MeshMaps, StructTmp, ErrStat,ErrMsg)  ! We are actually ignoring all errors from here
   type(MeshesMotionType), intent(inout)  :: MeshMotions
   type(MeshesLoadsType ), intent(inout)  :: MeshLoads
   type(MeshesMapsType  ), intent(inout)  :: MeshMaps
   type(StructTmpType   ), intent(inout)  :: StructTmp
   integer(IntKi),         intent(  out)  :: ErrStat
   character(ErrMsgLen),   intent(  out)  :: ErrMsg
   integer(IntKi)                         :: ErrStat2
   character(ErrMsgLen)                   :: ErrMsg2
   character(*),           parameter      :: RoutineName = 'WaveTank::StructDestroy'
   ErrStat = ErrID_None
   ErrMsg  = ''
   call WT_DestroyMeshesMotionType(MeshMotions, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call WT_DestroyMeshesLoadsType(MeshLoads, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call WT_DestroyMeshesMapsType(MeshMaps, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call WT_DestroyStructTmpType(StructTmp, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
end subroutine



!> Convert an Euler angle set of Roll, Pitch, Yaw ordering to a DCM.
!! this routine exists for two reasons
!!    1. ordering may be different
!!    2. incoming Euler angle is c_float instead of R8Ki
!! NOTE: no Euler angles are exported, so we stick with the OF convention
!!    for all internal conversions
function WT_EulerToDCM_fromInput(Ang) result(DCM)
   real(c_float), intent(in   )  :: Ang(3)
   real(R8Ki)                    :: DCM(3,3)
   !>>> Select one of the two following orders
   ! 3-2-1 intrinsic rotation sequence of the 3 Tait-Bryan angles (1-2-3 extrinsic rotation)
   DCM = EulerConstruct(real(Ang, DbKi))
   !! 1-2-3 intrinsic rotation sequence of the 3 Tait-Bryan angles (3-2-1 extrinsic rotation)
   !DCM = EulerConstructZYX(real(Ang, DbKi))
end function



subroutine WrVTK_Struct_Ref(SimSettings, MeshMotions, MeshLoads, ErrStat, ErrMsg)
   type(SimSettingsType),  target,  intent(in   )  :: SimSettings
   type(MeshesMotionType),          intent(in   )  :: MeshMotions
   type(MeshesLoadsType ),          intent(in   )  :: MeshLoads
   integer(IntKi),                  intent(  out)  :: ErrStat
   character(ErrMsgLen),            intent(  out)  :: ErrMsg
   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2
   character(*),           parameter               :: RoutineName = 'WaveTank::WrVTK_Struct_Ref'
   character(1024)                                 :: DirRootName
   real(SiKi)                                      :: RefPt(3)
   ErrStat = ErrID_None
   ErrMsg  = ''
   RefPt   = (/ 0.0_SiKi, 0.0_SiKi, 0.0_SiKi /)
   DirRootName = trim(SimSettings%Viz%WrVTK_dir)//PathSep//trim(SimSettings%Sim%OutRootName)
   ! Wave elevation measurement buoy
   call MeshWrVTKreference(RefPt, MeshMotions%WaveBuoyMotion, trim(DirRootName)//'.WaveBuoyMotion', ErrStat2, ErrMsg2); if (Failed()) return
   ! Platform point
   call MeshWrVTKreference(RefPt, MeshMotions%PtfmPtMotion, trim(DirRootName)//'.Struct'//'.PtfmPtMotion',   ErrStat2, ErrMsg2); if (Failed()) return
   ! Tower
   call MeshWrVTKreference(RefPt, MeshMotions%TowerMotion,  trim(DirRootName)//'.Struct'//'.TowerMotion',    ErrStat2, ErrMsg2); if (Failed()) return
contains
   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine

subroutine WrVTK_Struct(n_Global, SimSettings, MeshMotions, MeshLoads, ErrStat, ErrMsg)
   integer(IntKi),                  intent(in   )  :: n_Global
   type(SimSettingsType),  target,  intent(in   )  :: SimSettings
   type(MeshesMotionType),          intent(in   )  :: MeshMotions
   type(MeshesLoadsType ),          intent(in   )  :: MeshLoads
   integer(IntKi),                  intent(  out)  :: ErrStat
   character(ErrMsgLen),            intent(  out)  :: ErrMsg
   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2
   character(*),           parameter               :: RoutineName = 'WaveTank::WrVTK_Struct'
   character(1024)                                 :: DirRootName
   real(SiKi)                                      :: RefPt(3)
   ErrStat = ErrID_None
   ErrMsg  = ''
   RefPt   = (/ 0.0_SiKi, 0.0_SiKi, 0.0_SiKi /)
   DirRootName = trim(SimSettings%Viz%WrVTK_dir)//PathSep//trim(SimSettings%Sim%OutRootName)
   ! Wave elevation measurement buoy
   call MeshWrVTK(RefPt, MeshMotions%WaveBuoyMotion, trim(DirRootName)//'.WaveBuoyMotion', n_Global, .true., ErrStat2, ErrMsg2, Twidth=SimSettings%Viz%Twidth); if (Failed()) return
   ! Platform point
   call MeshWrVTK(RefPt, MeshMotions%PtfmPtMotion, trim(DirRootName)//'.Struct'//'.PtfmPtMotion',   n_Global, .true., ErrStat2, ErrMsg2, Twidth=SimSettings%Viz%Twidth); if (Failed()) return
   ! Tower
   call MeshWrVTK(RefPt, MeshMotions%TowerMotion,  trim(DirRootName)//'.Struct'//'.TowerMotion',    n_Global, .true., ErrStat2, ErrMsg2, Twidth=SimSettings%Viz%Twidth); if (Failed()) return
contains
   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine

end module
