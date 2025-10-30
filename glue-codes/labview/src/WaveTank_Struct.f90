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
   public :: StructCreateMeshMaps
   public :: StructDestroy
   public :: StructMotionUpdate
   public :: WrVTK_Struct_Ref
   public :: WrVTK_Struct
   public :: FroudeScaleM2F_Disp
   public :: FroudeScaleM2F_TVel
   public :: FroudeScaleM2F_RVel
   public :: FroudeScaleM2F_TAcc
   public :: FroudeScaleM2F_RAcc
   public :: FroudeScaleM2F_Time
   public :: FroudeScaleF2M_Frc
   public :: FroudeScaleF2M_Mom

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
   real(DbKi)                                      :: TmpAng(3)   ! temporary euler angle
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
   TmpPos = real(TrbCfg%PtfmRefPos, ReKi)
   Orient=WT_EulerToDCM_fromInput(TrbCfg%PtfmRefOrient)
   call CreateInputPointMesh(MeshMotions%PtfmPtMotion, TmpPos, Orient, ErrStat2, ErrMsg2, hasMotion=.true., hasLoads=.false.); if (Failed()) return
   MeshMotions%PtfmPtMotion%RemapFlag = .false.

   !-------------------------------
   ! create 2 point tower mesh
   call MeshCreate ( BlankMesh = MeshMotions%TowerMotion, IOS=COMPONENT_INPUT, Nnodes=2, ErrStat=ErrStat2, ErrMess=ErrMsg2,  &
               Orientation = .true., TranslationDisp = .true., TranslationVel = .true., TranslationAcc  = .TRUE. )
   if (Failed()) return

   ! Tower bottom
   TmpPos(1:2) = real(TrbCfg%TowerBsPt(1:2) + TrbCfg%PtfmRefPos(1:2), ReKi)      ! relative to PtfmRefPos in (x,y)
   TmpPos(3)   = real(TrbCfg%TowerBsPt(3), ReKi)                                 ! relative to MSL in (z)
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
   MeshMotions%TowerMotion%RemapFlag = .false.


   !-------------------------------
   ! create hub mesh
   TmpPos(1:3) = MeshMotions%TowerMotion%Position(1:3,2)                                        ! Tower top
   TmpPos(1)   = TmpPos(1) + cos(TrbInit%NacYaw) * TrbCfg%OverHang                              ! X, nacelle yaw, and overhang
   TmpPos(2)   = TmpPos(2) + sin(TrbInit%NacYaw) * TrbCfg%OverHang                              ! Y, nacelle yaw, and overhang
   TmpPos(3)   = TmpPos(3) + TrbCfg%Twr2Shft     + abs(TrbCfg%OverHang) * tan(TrbCfg%ShftTilt)  ! Z, shaft height above tower top, and shaft tilt

   TmpAng = (/ 0.0_DbKi, real(TrbCfg%ShftTilt,DbKi), real(TrbInit%NacYaw + Pi,DbKi)  /)       ! Hub/rotor azimuth is zero for reference. Rotate about Z by 180 so hub X points in negative X
   Orient = EulerConstruct(TmpAng)
   call CreateInputPointMesh(MeshMotions%HubMotion, TmpPos, Orient, ErrStat2, ErrMsg2, hasMotion=.true., hasLoads=.false.); if (Failed()) return
   MeshMotions%HubMotion%RemapFlag = .false.

contains
   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine

!> create mesh mappings
subroutine StructCreateMeshMaps(MeshMotions, MeshLoads, MeshMaps, ErrStat, ErrMsg)
   type(MeshesMotionType),          intent(inout)  :: MeshMotions
   type(MeshesLoadsType ),          intent(inout)  :: MeshLoads
   type(MeshesMapsType  ),          intent(inout)  :: MeshMaps
   integer(IntKi),                  intent(  out)  :: ErrStat
   character(ErrMsgLen),            intent(  out)  :: ErrMsg
   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2
   character(*),           parameter               :: RoutineName = 'WaveTank::StructCreateMeshMaps'
   ErrStat = ErrID_None
   ErrMsg  = ''
   !-------------------------------
   ! Mesh motion mappings
   call MeshMapCreate(MeshMotions%PtfmPtMotion, MeshMotions%TowerMotion, MeshMaps%Motion_PRP_2_Twr, errStat2, errMsg2); if(Failed())return
   call MeshMapCreate(MeshMotions%PtfmPtMotion, MeshMotions%HubMotion,   MeshMaps%Motion_PRP_2_Hub, errStat2, errMsg2); if(Failed())return
 
contains
   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine



!> updates the structural meshes
subroutine StructMotionUpdate(SimSettings,CalcStepIO, MeshMotions, MeshMaps, StructTmp, ErrStat, ErrMsg)
   type(SimSettingsType),  target,  intent(in   )  :: SimSettings
   type(CalcStepIOdataType),        intent(in   )  :: CalcStepIO
   type(MeshesMotionType), target,  intent(inout)  :: MeshMotions
   type(MeshesMapsType  ),          intent(inout)  :: MeshMaps
   type(StructTmpType   ),          intent(inout)  :: StructTmp
   integer(IntKi),                  intent(  out)  :: ErrStat
   character(ErrMsgLen),            intent(  out)  :: ErrMsg
   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2
   character(*),           parameter               :: RoutineName = 'WaveTank::StructMotionUpdate'
   real(R8Ki)                                      :: TmpTransDisp(3)
   real(DbKi)                                      :: TmpAng(3)   ! temporary euler angle
   real(R8Ki)                                      :: Orient(3,3) ! temporary orientation
   type(TurbConfigType),   pointer                 :: TrbCfg      ! to shorten notation
   type(TurbInitCondType), pointer                 :: TrbInit     ! to shorten notation
   type(MeshType),         pointer                 :: Ptfm        ! to shorten notation
   type(MeshType),         pointer                 :: Twr         ! to shorten notation
   type(MeshType),         pointer                 :: Hub         ! to shorten notation
   real(c_float)                                   :: ScaleFact   ! to shorten notation
   ErrStat = ErrID_None
   ErrMsg  = ''

   TrbCfg    => SimSettings%TrbCfg
   TrbInit   => SimSettings%TrbInit

   ! scaling factor
   ScaleFact = SimSettings%Sim%ScaleFact

   ! update PtfmPtMotion
   Ptfm => MeshMotions%PtfmPtMotion
   Ptfm%TranslationDisp(1:3,1) = FroudeScaleM2F_Disp(ScaleFact, CalcStepIO%PosAng_c(1:3), Ptfm%Position(1:3,1))
   Ptfm%Orientation(1:3,1:3,1) = WT_EulerToDCM_fromInput(CalcStepIO%PosAng_c(4:6))        ! angles don't scale
   Ptfm%TranslationVel(1:3,1)  = FroudeScaleM2F_TVel(ScaleFact, CalcStepIO%Vel_c(1:3))
   Ptfm%RotationVel(1:3,1)     = FroudeScaleM2F_RVel(ScaleFact, CalcStepIO%Vel_c(4:6))
   Ptfm%TranslationAcc(1:3,1)  = FroudeScaleM2F_TAcc(ScaleFact, CalcStepIO%Acc_c(1:3))
   Ptfm%RotationAcc(1:3,1)     = FroudeScaleM2F_RAcc(ScaleFact, CalcStepIO%Acc_c(4:6))

   ! transfer Ptfm to Tower
   Twr => MeshMotions%TowerMotion
   call Transfer_Point_to_Line2( Ptfm, Twr, MeshMaps%Motion_PRP_2_Twr, ErrStat2, ErrMsg2 ); if (Failed()) return;

   ! transfer Ptfm to hub (tower is rigid)
   Hub => MeshMotions%HubMotion
   call Transfer_Point_to_Point( Ptfm, Hub, MeshMaps%Motion_PRP_2_Hub, ErrStat2, ErrMsg2 ); if (Failed()) return;

   ! rotor azimuth
   StructTmp%Azimuth = modulo(real(CalcStepIO%Time_c,ReKi)*StructTmp%RotSpeed + TrbInit%Azimuth, TwoPi )
  
   ! update hub azimuth -- include initial azimuth
   TmpAng = (/ real(StructTmp%Azimuth,DbKi), 0.0_DbKi, 0.0_DbKi /)
   Orient = EulerConstruct(TmpAng)
   Hub%Orientation(1:3,1:3,1) = matmul(Orient,Hub%Orientation(1:3,1:3,1))


   ! hub to blades



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
   call MeshWrVTKreference(RefPt, MeshMotions%PtfmPtMotion, trim(DirRootName)//'.Struct'//'.PtfmPtMotion', ErrStat2, ErrMsg2); if (Failed()) return
   ! Tower
   call MeshWrVTKreference(RefPt, MeshMotions%TowerMotion,  trim(DirRootName)//'.Struct'//'.TowerMotion',  ErrStat2, ErrMsg2); if (Failed()) return
   ! hub point
   call MeshWrVTKreference(RefPt, MeshMotions%HubMotion,    trim(DirRootName)//'.Struct'//'.HubMotion',    ErrStat2, ErrMsg2); if (Failed()) return
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
   call MeshWrVTK(RefPt, MeshMotions%PtfmPtMotion, trim(DirRootName)//'.Struct'//'.PtfmPtMotion', n_Global, .true., ErrStat2, ErrMsg2, Twidth=SimSettings%Viz%Twidth); if (Failed()) return
   ! Tower
   call MeshWrVTK(RefPt, MeshMotions%TowerMotion,  trim(DirRootName)//'.Struct'//'.TowerMotion',  n_Global, .true., ErrStat2, ErrMsg2, Twidth=SimSettings%Viz%Twidth); if (Failed()) return
   ! Hub point
   call MeshWrVTK(RefPt, MeshMotions%HubMotion,    trim(DirRootName)//'.Struct'//'.HubMotion',    n_Global, .true., ErrStat2, ErrMsg2, Twidth=SimSettings%Viz%Twidth); if (Failed()) return
contains
   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine



!-----------------------------------------------
! Froude scaling from here: https://home.hvl.no/ansatte/gste/ftp/MarinLab_files/Litteratur/NTNU_Scaling_Laws.pdf, page 21
! notation below:
!     model scale: _m
!     full scale:  _f
!     ScaleFact:   length_f/length_m = lambda
!     DensFact:    rho_f/rho_m

!> scale model displacements to full scale 
!! length_full = length_model * lambda
function FroudeScaleM2F_Disp(ScaleFact, Pos_m, refPos_f) result(transDisp_f)
   real(c_float), intent(in   ) :: ScaleFact
   real(c_float), intent(in   ) :: Pos_m(3)
   real(ReKi),    intent(in   ) :: refPos_f(3)
   real(R8Ki)                   :: transdisp_f(3)
   transDisp_f = real(ScaleFact*Pos_m,R8Ki) - real(refPos_f,R8Ki)
end function

!> scale model translational velocity to full scale 
!! TVel_full = TVel_model * sqrt( lambda )         TODO: check this!!!! 
function FroudeScaleM2F_TVel(ScaleFact, TVel_m) result(TVel_f)
   real(c_float), intent(in   ) :: ScaleFact
   real(c_float), intent(in   ) :: TVel_m(3)
   real(ReKi)                   :: TVel_f(3)
   TVel_f = sqrt(real(ScaleFact,ReKi)) * real(TVel_m,ReKi)
end function

!> scale model rotational velocity to full scale 
!! RVel_full = RVel_model * sqrt(lambda)           TODO: check this!!!! 
function FroudeScaleM2F_RVel(ScaleFact, RVel_m) result(RVel_f)
   real(c_float), intent(in   ) :: ScaleFact
   real(c_float), intent(in   ) :: RVel_m(3)
   real(ReKi)                   :: RVel_f(3)
   RVel_f =  real(RVel_m,ReKi) / sqrt(real(ScaleFact,ReKi))
end function

!> scale model translational acceleration to full scale 
!! TAcc_full = TAcc_model ---> no scaling applied
function FroudeScaleM2F_TAcc(ScaleFact, TAcc_m) result(TAcc_f)
   real(c_float), intent(in   ) :: ScaleFact
   real(c_float), intent(in   ) :: TAcc_m(3)
   real(ReKi)                   :: TAcc_f(3)
   TAcc_f = real(TAcc_m,ReKi)
end function

!> scale model rotational acceleration to full scale 
!! RAcc_full = RAcc_model / lambda                 TODO: check this!!!! 
function FroudeScaleM2F_RAcc(ScaleFact, RAcc_m) result(RAcc_f)
   real(c_float), intent(in   ) :: ScaleFact
   real(c_float), intent(in   ) :: RAcc_m(3)
   real(ReKi)                   :: RAcc_f(3)
   RAcc_f = real(RAcc_m,ReKi) / real(ScaleFact,ReKi)
end function

!> scale model time to full scale
!! sqrt(lambda) = sqrt(Length_full/ Length_model)
function FroudeScaleM2F_Time(ScaleFact, Time_m) result(Time_f)
   real(c_float), intent(in   ) :: ScaleFact
   real(c_double),intent(in   ) :: Time_m
   real(R8Ki)                   :: Time_f
   Time_f = sqrt(real(ScaleFact,R8Ki)) * real(Time_m,R8Ki)
end function

!> scale full scale force to model
!! lambda^3 * DensFact * Frc_model = Frc_full
function FroudeScaleF2M_Frc(ScaleFact, DensFact, Frc_f) result(Frc_m)
   real(c_float), intent(in   ) :: ScaleFact
   real(c_float), intent(in   ) :: DensFact
   real(ReKi),    intent(in   ) :: Frc_f(3)
   real(c_float)                :: Frc_m(3)
   Frc_m = real(Frc_f, c_float) / (ScaleFact**3 * DensFact)
end function

!> scale full scale moment to model
!! lambda^4 * DensFact * Mom_model = Mom_full
function FroudeScaleF2M_Mom(ScaleFact, DensFact, Mom_f) result(Mom_m)
   real(c_float), intent(in   ) :: ScaleFact
   real(c_float), intent(in   ) :: DensFact
   real(ReKi),    intent(in   ) :: Mom_f(3)
   real(c_float)                :: Mom_m(3)
   Mom_m = real(Mom_f, c_float) / (ScaleFact**4 * DensFact)
end function


end module
