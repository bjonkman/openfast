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


   ! create PRP platform point
   Orient=WT_EulerToDCM_fromInput(TrbCfg%PtfmRefOrient)
   call CreateInputPointMesh(MeshMotions%PtfmPtMotion, TrbCfg%PtfmRefPos, Orient, ErrStat2, ErrMsg2, hasMotion=.true., hasLoads=.false.); if (Failed()) return

   !



   ! Rotor

 
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
   real(ReKi)                                      :: RefPt(3)
   ErrStat = ErrID_None
   ErrMsg  = ''
   RefPt   = (/ 0.0_ReKi, 0.0_ReKi, 0.0_ReKi /)
   ! Platform point
   call MeshWrVTKreference(RefPt, MeshMotions%PtfmPtMotion, trim(SimSettings%Viz%WrVTK_dir)//PathSep//trim(SimSettings%Sim%OutRootName)//'.Struct'//'.PtfmPtMotion', ErrStat2, ErrMsg2); if (Failed()) return
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
   real(ReKi)                                      :: RefPt(3)
   ErrStat = ErrID_None
   ErrMsg  = ''
   RefPt   = (/ 0.0_ReKi, 0.0_ReKi, 0.0_ReKi /)
   ! Platform point
   call MeshWrVTK(RefPt, MeshMotions%PtfmPtMotion, trim(SimSettings%Viz%WrVTK_dir)//PathSep//trim(SimSettings%Sim%OutRootName)//'.Struct'//'.PtfmPtMotion', n_Global, .true., ErrStat2, ErrMsg2, Twidth=SimSettings%Viz%Twidth); if (Failed()) return
contains
   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine

end module
