!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2016-2018  Envision Energy USA, LTD
!
!    This file is part of AeroDyn.
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
!
!**********************************************************************************************************************************
module AeroDyn_Driver_Subs
   
   use AeroDyn_Driver_Types   
   use AeroDyn
   use VersionInfo

   implicit none   
   
   TYPE(ProgDesc), PARAMETER   :: version   = ProgDesc( 'AeroDyn_driver', '', '' )  ! The version number of this program.
                                                    
   contains

!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_Init(DvrData,errStat,errMsg )

   type(Dvr_SimData),            intent(  out) :: DvrData       ! driver data
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None

      ! local variables
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   character(*), parameter                     :: RoutineName = 'Dvr_Init'

   CHARACTER(1000)                             :: inputFile     ! String to hold the file name.
   CHARACTER(200)                              :: git_commit    ! String containing the current git commit hash
   CHARACTER(20)                               :: FlagArg       ! flag argument from command line


   ErrStat = ErrID_None
   ErrMsg  = ""


   DvrData%OutFileData%unOutFile   = -1
   
   CALL NWTC_Init( ProgNameIN=version%Name )

   InputFile = ""  ! initialize to empty string to make sure it's input from the command line
   CALL CheckArgs( InputFile, Flag=FlagArg )
   IF ( LEN( TRIM(FlagArg) ) > 0 ) CALL NormStop()

      ! Display the copyright notice
   CALL DispCopyrightLicense( version%Name )
      ! Obtain OpenFAST git commit hash
   git_commit = QueryGitVersion()
      ! Tell our users what they're running
   CALL WrScr( ' Running '//TRIM( version%Name )//' a part of OpenFAST - '//TRIM(git_Commit)//NewLine//' linked with '//TRIM( NWTC_Ver%Name )//NewLine )
         
      ! Read the AeroDyn driver input file
   call Dvr_ReadInputFile(inputFile, DvrData, errStat2, errMsg2 )
      call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName) 
      if (errStat >= AbortErrLev) return
      
      ! validate the inputs
   call ValidateInputs(DvrData, errStat2, errMsg2)      
      call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName) 
      
end subroutine Dvr_Init 
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Init_AeroDyn(iCase, DvrData, AD, dt, errStat, errMsg)

   integer(IntKi),               intent(in   ) :: iCase         ! driver case
   type(Dvr_SimData),            intent(inout) :: DvrData       ! Input data for initialization (intent out for getting AD WriteOutput names/units)
   type(AeroDyn_Data),           intent(inout) :: AD            ! AeroDyn data 
   real(DbKi),                   intent(inout) :: dt            ! interval
      
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None

      ! locals
   real(reKi)                                  :: theta(3)
   integer(IntKi)                              :: j, k   
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   character(*), parameter                     :: RoutineName = 'Init_AeroDyn'
                                                  
   ! local data                                
   type(AD_InitInputType)                      :: InitInData     ! Input data for initialization
   type(AD_InitOutputType)                     :: InitOutData    ! Output data from initialization
   real(ReKi)                                  :: RotAzimuth    ! Rotor azimuth -- aligned with blade 1 (deg)

   errStat = ErrID_None
   errMsg  = ''


   allocate(InitInData%rotors(1), InitOutData%rotors(1), stat=errStat) 
   if (errStat/=0) then
      call SetErrStat( ErrID_Fatal, 'Allocating rotors', errStat, errMsg, RoutineName )
      call Cleanup()
      return
   end if
      
   
   if (iCase.EQ.1) then
   
      InitInData%InputFile      = DvrData%AD_InputFile
      InitInData%RootName       = DvrData%outFileData%Root
      InitInData%rotors(1)%NumBlades      = DvrData%numBlades
      InitInData%Gravity        = 9.80665_ReKi
      InitInData%Linearize      = .false.
                        
   
         ! set initialization data:
      call AllocAry( InitInData%rotors(1)%BladeRootPosition, 3, InitInData%rotors(1)%NumBlades, 'BladeRootPosition', errStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      call AllocAry( InitInData%rotors(1)%BladeRootOrientation, 3, 3, InitInData%rotors(1)%NumBlades, 'BladeRootOrientation', errStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
      
      call eye(InitInData%rotors(1)%NacelleOrientation, ErrStat2, ErrMsg2) ! nacelle reference orientation will be identity
      
      InitInData%rotors(1)%HubPosition = (/ DvrData%Overhang * cos(DvrData%shftTilt), 0.0_ReKi, DvrData%HubHt /)
      theta(1) = 0.0_ReKi
      theta(2) = -DvrData%shftTilt
      theta(3) = 0.0_ReKi
      InitInData%rotors(1)%HubOrientation = EulerConstruct( theta )
     
   
      do k=1,InitInData%rotors(1)%numBlades
                     
         theta(1) = (k-1)*TwoPi/real(InitInData%rotors(1)%numBlades,ReKi)
         theta(2) = DvrData%precone
         theta(3) = 0.0_ReKi
         InitInData%rotors(1)%BladeRootOrientation(:,:,k) = matmul( EulerConstruct( theta ), InitInData%rotors(1)%HubOrientation )
                  
         InitInData%rotors(1)%BladeRootPosition(:,k)   = InitInData%rotors(1)%HubPosition + DvrData%hubRad * InitInData%rotors(1)%BladeRootOrientation(3,:,k)      
      
      end do
      
    
      call AD_Init(InitInData, AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, dt, InitOutData, ErrStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )


      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if   

         
      do j = 2, numInp
         call AD_CopyInput (AD%u(1),  AD%u(j),  MESH_NEWCOPY, errStat2, errMsg2)
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      end do
      
         ! move AD initOut data to AD Driver
      call move_alloc( InitOutData%rotors(1)%WriteOutputHdr, DvrData%OutFileData%WriteOutputHdr )
      call move_alloc( InitOutData%rotors(1)%WriteOutputUnt, DvrData%OutFileData%WriteOutputUnt )   
     
      DvrData%OutFileData%AD_ver = InitOutData%ver
             
      call cleanup() ! destroy init input/output data
      
   else
   
      call AD_ReInit(AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%m, dt, ErrStat2, ErrMsg2 )   
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         if (ErrStat >= AbortErrLev) return
      
   end if   
   
      
      ! we know exact values, so we're going to initialize inputs this way (instead of using the input guesses from AD_Init)
   AD%InputTime = -999
   RotAzimuth = 0.0
   DO j = 1-numInp, 0
      call Set_AD_Inputs(iCase,j,RotAzimuth,DvrData,AD,errStat2,errMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END DO              
   
      
   
contains

   subroutine cleanup()
      call AD_DestroyInitInput( InitInData, ErrStat2, ErrMsg2 )   
      call AD_DestroyInitOutput( InitOutData, ErrStat2, ErrMsg2 )      
   end subroutine cleanup
   
end subroutine Init_AeroDyn
!----------------------------------------------------------------------------------------------------------------------------------
!> this routine cycles values in the input array AD%InputTime and AD%u.
!! it then sets the inputs for nt * DvrData%Cases(iCase)%dT, which are index values 1 in the arrays.
subroutine Set_AD_Inputs(iCase,nt,RotAzimuth,DvrData,AD,errStat,errMsg)

   integer(IntKi)              , intent(in   ) :: iCase         ! case number 
   integer(IntKi)              , intent(in   ) :: nt            ! time step number

   real(ReKi)                  , intent(inout) :: RotAzimuth    ! Rotor azimuth at time nt-1 -- aligned with blade 1 (deg)
   type(Dvr_SimData),            intent(in   ) :: DvrData       ! Driver data 
   type(AeroDyn_Data),           intent(inout) :: AD            ! AeroDyn data 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None

      ! local variables
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   character(*), parameter                     :: RoutineName = 'Set_AD_Inputs'

   integer(intKi)                              :: j             ! loop counter for nodes
   integer(intKi)                              :: k             ! loop counter for blades
   integer(intKi)                              :: timeIndex     ! index for time

   real(ReKi)                                  :: z             ! height (m)
   !real(ReKi)                                  :: angle
   real(R8Ki)                                  :: theta(3)
   real(R8Ki)                                  :: position(3)
   real(R8Ki)                                  :: orientation(3,3)
   real(R8Ki)                                  :: rotateMat(3,3)
   
   
   errStat = ErrID_None
   errMsg  = ""
   
      ! note that this initialization is a little different than the general algorithm in FAST because here
      ! we can get exact values, so we are going to ignore initial guesses and not extrapolate
      
   timeIndex = min( max(1,nt+1), DvrData%Cases(iCase)%numSteps )
   
   !................
   ! shift previous calculations:
   !................
   do j = numInp-1,1,-1
      call AD_CopyInput (AD%u(j),  AD%u(j+1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
      AD%InputTime(j+1) = AD%InputTime(j)
   end do
   
   if (nt <= 0) then
         ! save the azimuth at t (not t+dt) for output to file:
         ! compare to theta(1) for calculate of HubMotion%Orientation below
      RotAzimuth = MODULO( REAL( DvrData%Cases(iCase)%dT * (nt-1) * DvrData%Cases(iCase)%RotSpeed(1), ReKi) * R2D, 360.0_ReKi )

      AD%inputTime(1) = DvrData%Cases(iCase)%time(1) + DvrData%Cases(iCase)%dT * nt ! time at nt+1
   else
      
      if (nt==1) then
         RotAzimuth = 0.0_ReKi
      else
         RotAzimuth = MODULO( RotAzimuth + REAL(DvrData%Cases(iCase)%dt * DvrData%Cases(iCase)%RotSpeed(nt), ReKi) * R2D, 360.0_ReKi ) ! add a delta angle to the previous azimuth
      end if
      
      if (nt == DvrData%Cases(iCase)%numSteps) then
         AD%inputTime(1) = DvrData%Cases(iCase)%time(timeIndex) + DvrData%Cases(iCase)%dT
      else
         AD%inputTime(1) = DvrData%Cases(iCase)%time(timeIndex)
      end if

   end if
         
   !................
   ! calculate new values
   !................
   
      ! Tower motions:
      do j=1,AD%u(1)%rotors(1)%TowerMotion%nnodes
         AD%u(1)%rotors(1)%TowerMotion%Orientation(  :,:,j) = AD%u(1)%rotors(1)%TowerMotion%RefOrientation(:,:,j) ! identity
         AD%u(1)%rotors(1)%TowerMotion%TranslationDisp(:,j) = 0.0_ReKi
         AD%u(1)%rotors(1)%TowerMotion%TranslationVel( :,j) = 0.0_ReKi
      end do !j=nnodes
      
      ! Nacelle motions:
      theta(1) = 0.0_ReKi
      theta(2) = 0.0_ReKi
      theta(3) = DvrData%Cases(iCase)%Yaw(timeIndex)
      orientation = EulerConstruct(theta)
      
      if (AD%u(1)%rotors(1)%NacelleMotion%Nnodes > 0) then
         AD%u(1)%rotors(1)%NacelleMotion%TranslationDisp = 0.0_R8Ki
         AD%u(1)%rotors(1)%NacelleMotion%TranslationVel = 0.0_R8Ki
         AD%u(1)%rotors(1)%NacelleMotion%Orientation(:,:,1) = orientation
      endif
            
      ! Hub motions:
      ! orientation set in nacelle motion calculation
      AD%u(1)%rotors(1)%HubMotion%TranslationDisp(:,1) = matmul( AD%u(1)%rotors(1)%HubMotion%Position(:,1), orientation ) - AD%u(1)%rotors(1)%HubMotion%Position(:,1) ! = matmul( transpose(orientation) - eye(3), AD%u(1)%HubMotion%Position(:,1) )

      theta(1) = RotAzimuth*D2R + DvrData%Cases(iCase)%dt * DvrData%Cases(iCase)%RotSpeed(timeIndex)  ! AD%inputTime(1) * DvrData%Cases(iCase)%RotSpeed
      theta(2) = 0.0_ReKi
      theta(3) = 0.0_ReKi
      AD%u(1)%rotors(1)%HubMotion%Orientation(  :,:,1) = matmul( AD%u(1)%rotors(1)%HubMotion%RefOrientation(:,:,1), orientation )
      orientation = EulerConstruct( theta )
      AD%u(1)%rotors(1)%HubMotion%Orientation(  :,:,1) = matmul( orientation, AD%u(1)%rotors(1)%HubMotion%Orientation(  :,:,1) )
      
      AD%u(1)%rotors(1)%HubMotion%RotationVel(    :,1) = AD%u(1)%rotors(1)%HubMotion%Orientation(1,:,1) * DvrData%Cases(iCase)%RotSpeed(timeIndex)
                  
      ! Blade motions:
      do k=1,DvrData%numBlades         
         theta(1) = (k-1)*TwoPi/real(DvrData%numBlades,ReKi)
         theta(2) =  DvrData%precone
         theta(3) = -DvrData%Cases(iCase)%pitch(timeIndex)
         orientation = EulerConstruct(theta)
         
         AD%u(1)%rotors(1)%BladeRootMotion(k)%Orientation(  :,:,1) = matmul( orientation, AD%u(1)%rotors(1)%HubMotion%Orientation(  :,:,1) )
      end do !k=numBlades
            
      ! Blade and blade root motions:
      do k=1,DvrData%numBlades
         rotateMat = transpose( AD%u(1)%rotors(1)%BladeRootMotion(k)%Orientation(  :,:,1) )
         rotateMat = matmul( rotateMat, AD%u(1)%rotors(1)%BladeRootMotion(k)%RefOrientation(  :,:,1) )
         orientation = transpose(rotateMat)
         
         rotateMat(1,1) = rotateMat(1,1) - 1.0_ReKi
         rotateMat(2,2) = rotateMat(2,2) - 1.0_ReKi
         rotateMat(3,3) = rotateMat(3,3) - 1.0_ReKi
                  

         position = AD%u(1)%rotors(1)%BladeRootMotion(k)%Position(:,1) - AD%u(1)%rotors(1)%HubMotion%Position(:,1) 
         AD%u(1)%rotors(1)%BladeRootMotion(k)%TranslationDisp(:,1) = AD%u(1)%rotors(1)%HubMotion%TranslationDisp(:,1) + matmul( rotateMat, position )

         position =  AD%u(1)%rotors(1)%BladeRootMotion(k)%Position(:,1) + AD%u(1)%rotors(1)%BladeRootMotion(k)%TranslationDisp(:,1) &
                     - AD%u(1)%rotors(1)%HubMotion%Position(:,1) - AD%u(1)%rotors(1)%HubMotion%TranslationDisp(:,1)
         AD%u(1)%rotors(1)%BladeRootMotion(k)%TranslationVel( :,1) = cross_product( AD%u(1)%rotors(1)%HubMotion%RotationVel(:,1), position )

         do j=1,AD%u(1)%rotors(1)%BladeMotion(k)%nnodes        
            position = AD%u(1)%rotors(1)%BladeMotion(k)%Position(:,j) - AD%u(1)%rotors(1)%HubMotion%Position(:,1) 
            AD%u(1)%rotors(1)%BladeMotion(k)%TranslationDisp(:,j) = AD%u(1)%rotors(1)%HubMotion%TranslationDisp(:,1) + matmul( rotateMat, position )
            
            AD%u(1)%rotors(1)%BladeMotion(k)%Orientation(  :,:,j) = matmul( AD%u(1)%rotors(1)%BladeMotion(k)%RefOrientation(:,:,j), orientation )
            
            
            position =  AD%u(1)%rotors(1)%BladeMotion(k)%Position(:,j) + AD%u(1)%rotors(1)%BladeMotion(k)%TranslationDisp(:,j) &
                      - AD%u(1)%rotors(1)%HubMotion%Position(:,1) - AD%u(1)%rotors(1)%HubMotion%TranslationDisp(:,1)
            AD%u(1)%rotors(1)%BladeMotion(k)%TranslationVel( :,j) = cross_product( AD%u(1)%rotors(1)%HubMotion%RotationVel(:,1), position )
            
            AD%u(1)%rotors(1)%BladeMotion(k)%RotationVel(:,j) = AD%u(1)%rotors(1)%HubMotion%Orientation(1,:,1) * DvrData%Cases(iCase)%RotSpeed(timeIndex) ! simplification (without pitch rate)
            AD%u(1)%rotors(1)%BladeMotion(k)%TranslationAcc(:,j) = 0.0_ReKi ! simplification
         end do !j=nnodes
                                    
      end do !k=numBlades       
      
      ! Inflow wind velocities:
      ! InflowOnBlade
      do k=1,DvrData%numBlades
         do j=1,AD%u(1)%rotors(1)%BladeMotion(k)%nnodes
            z = AD%u(1)%rotors(1)%BladeMotion(k)%Position(3,j) + AD%u(1)%rotors(1)%BladeMotion(k)%TranslationDisp(3,j)
            AD%u(1)%rotors(1)%InflowOnBlade(1,j,k) = GetU(  DvrData%Cases(iCase)%WndSpeed(timeIndex), DvrData%HubHt, DvrData%Cases(iCase)%ShearExp(timeIndex), z )
            AD%u(1)%rotors(1)%InflowOnBlade(2,j,k) = 0.0_ReKi !V
            AD%u(1)%rotors(1)%InflowOnBlade(3,j,k) = 0.0_ReKi !W      
         end do !j=nnodes
      end do !k=numBlades
      
      !InflowOnTower
      do j=1,AD%u(1)%rotors(1)%TowerMotion%nnodes
         z = AD%u(1)%rotors(1)%TowerMotion%Position(3,j) + AD%u(1)%rotors(1)%TowerMotion%TranslationDisp(3,j)
         AD%u(1)%rotors(1)%InflowOnTower(1,j) = GetU(  DvrData%Cases(iCase)%WndSpeed(timeIndex), DvrData%HubHt, DvrData%Cases(iCase)%ShearExp(timeIndex), z )
         AD%u(1)%rotors(1)%InflowOnTower(2,j) = 0.0_ReKi !V
         AD%u(1)%rotors(1)%InflowOnTower(3,j) = 0.0_ReKi !W         
      end do !j=nnodes
      
      !InflowOnNacelle
      if (AD%u(1)%rotors(1)%NacelleMotion%Committed) then
         z = AD%u(1)%rotors(1)%NacelleMotion%Position(3,1) + AD%u(1)%rotors(1)%NacelleMotion%TranslationDisp(3,1)
         AD%u(1)%rotors(1)%InflowOnNacelle(1) = GetU( DvrData%Cases(iCase)%WndSpeed(timeIndex), DvrData%HubHt, DvrData%Cases(iCase)%ShearExp(timeIndex), z )
      else
         AD%u(1)%rotors(1)%InflowOnNacelle(2) = 0.0_ReKi ! U
      end if
      AD%u(1)%rotors(1)%InflowOnNacelle(2) = 0.0_ReKi !V
      AD%u(1)%rotors(1)%InflowOnNacelle(3) = 0.0_ReKi !W
                     
end subroutine Set_AD_Inputs
!----------------------------------------------------------------------------------------------------------------------------------
function GetU( URef, ZRef, PLExp, z ) result (U)
   real(ReKi), intent(in) :: URef
   real(ReKi), intent(in) :: ZRef
   real(ReKi), intent(in) :: PLExp
   real(ReKi), intent(in) :: z
   real(ReKi)             :: U
   
   U = URef*(z/ZRef)**PLExp

end function GetU
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_ReadInputFile(fileName, DvrData, errStat, errMsg )
   ! This routine opens the gets the data from the input files.

   character(*),                  intent( in    )   :: fileName
   type(Dvr_SimData),             intent(   out )   :: DvrData
   integer,                       intent(   out )   :: errStat              ! returns a non-zero value when an error occurs  
   character(*),                  intent(   out )   :: errMsg               ! Error message if errStat /= ErrID_None
   

      ! Local variables
   character(1024)              :: PriPath
   character(1024)              :: inpVersion                               ! String containing the input-version information.
   character(1024)              :: Line                                     ! String containing a line of input.
   integer                      :: unIn, unEc
   integer                      :: ICase
   integer                      :: nt, Ind
   integer                      :: Sttus
   character( 11)               :: DateNow                                  ! Date shortly after the start of execution.
   character(  8)               :: TimeNow                                  ! Time of day shortly after the start of execution.
   
   integer, parameter           :: NumCols = 7                              ! number of columns to be read from the input file
   real(DbKi)                   :: InpCase(NumCols)                         ! Temporary array to hold combined-case input parameters. (note that we store in double precision so the time is read correctly)
   logical                      :: TabDel      
   logical                      :: echo   

   INTEGER(IntKi)               :: ErrStat2                                 ! Temporary Error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                                  ! Temporary Err msg
   CHARACTER(*), PARAMETER      :: RoutineName = 'Dvr_ReadInputFile'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ''
   UnIn = -1
   UnEc = -1
   
   ! Open the input file
   CALL GetPath( fileName, PriPath )     ! Input files will be relative to the path where the primary input file is located.

   call GetNewUnit( unIn )   
   call OpenFInpFile( unIn, fileName, errStat2, ErrMsg2 )
   call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   if ( errStat >= AbortErrLev ) then
      call cleanup()
      return
   end if

   
   call WrScr( 'Opening input file:  '//fileName )

      ! Skip a line, read the run title information.

   CALL ReadStr( UnIn, fileName, inpVersion, 'inpVersion', 'File Header: (line 1)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   CALL ReadStr( UnIn, fileName, DvrData%OutFileData%runTitle, 'runTitle', 'File Header: File Description (line 2)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   call WrScr1 ( ' '//DvrData%OutFileData%runTitle )
   
      ! Read in the title line for the input-configuration subsection.
   CALL ReadStr( UnIn, fileName, line, 'line', 'File Header: (line 3)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


      ! See if we should echo the output.     
   call ReadVar ( unIn, fileName, echo, 'Echo', 'Echo Input', errStat2, errMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   if ( echo )  then
         ! Get date and time.
      dateNow = CurDate()
      timeNow = CurTime()
      call GetNewUnit( unEc ) 
      call getroot(fileName,DvrData%OutFileData%Root)      
      call  OpenFOutFile ( unEc, trim( DvrData%OutFileData%Root )//'.ech', errStat2, errMsg2 )
         call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
         if ( errStat >= AbortErrLev ) then
            call cleanup()
            return
         end if
      
      write (unEc,'(A)')      'Echo of Input File:'
      write (unEc,'(A)')      ' "'//fileName//'"'
      write (unEc,'(A)')      'Generated on: '//trim( dateNow )//' at '//trim( timeNow )//'.'
      write (unEc,'(A)')      inpVersion
      write (unEc,'(A)')      DvrData%OutFileData%runTitle
      write (unEc,'(A)')      line
      write (unEc,Ec_LgFrmt)  echo, 'Echo', 'Echo input parameters to "rootname.ech"?'
   end if


      ! Read the rest of input-configuration section.
      
   call ReadVar ( unIn, fileName, DvrData%AD_InputFile,   'AD_InputFile',   'Name of the AeroDyn input file', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      if ( errStat >= AbortErrLev ) then
         call cleanup()
         return
      end if
   IF ( PathIsRelative( DvrData%AD_InputFile ) ) DvrData%AD_InputFile = TRIM(PriPath)//TRIM(DvrData%AD_InputFile)

   
      ! Read the turbine-data section.

   call ReadCom ( unIn, fileName, 'the turbine-data subtitle', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%numBlades,'NumBlades','Number of blades', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%HubRad,   'HubRad',   'Hub radius (m)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%HubHt,    'HubHt',    'Hub height (m)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%Overhang, 'Overhang',  'Overhang (m)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%ShftTilt, 'ShftTilt',  'Shaft tilt (deg)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      DvrData%ShftTilt = DvrData%ShftTilt*D2R
   call ReadVar ( unIn, fileName, DvrData%precone, 'Precone',  'Precone (deg)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      DvrData%precone = DvrData%precone*D2R
            
      if ( errStat >= AbortErrLev ) then
         call cleanup()
         return
      end if           


      ! Read the I/O-configuration section.

   call ReadCom ( unIn, fileName, 'the I/O-configuration subtitle', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%OutFileData%Root, 'OutFileRoot', 'Root name for any output files', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      IF ( PathIsRelative( DvrData%OutFileData%Root ) ) DvrData%OutFileData%Root = TRIM(PriPath)//TRIM(DvrData%OutFileData%Root)
   if (len_trim(DvrData%OutFileData%Root) == 0) then
      call getroot(fileName,DvrData%OutFileData%Root)
   end if
   
   call ReadVar ( unIn, fileName, TabDel,   'TabDel',   'Make output tab-delimited (fixed-width otherwise)?', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      if (TabDel) then
         DvrData%OutFileData%delim = TAB
      else
         DvrData%OutFileData%delim = " "
      end if
               
      ! OutFmt - Format used for text tabular output (except time).  Resulting field should be 10 characters. (-):
   call ReadVar( UnIn, fileName, DvrData%OutFileData%OutFmt, "OutFmt", "Format used for text tabular output (except time).  Resulting field should be 10 characters. (-)", ErrStat2, ErrMsg2, UnEc)
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName ) !bjj: this is a global variable in NWTC_Library            
   call ReadVar ( unIn, fileName, Beep,  'Beep',     'Beep on exit?', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName ) !bjj: this is a global variable in NWTC_Library
      if ( errStat >= AbortErrLev ) then
         call cleanup()
         return
      end if



      ! Read the combined-case section.

   call ReadCom  ( unIn, fileName, 'the combined-case subtitle', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar  ( unIn, fileName, DvrData%NumCases, 'NumCases', 'Number of cases to run', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadCom  ( unIn, fileName, 'the combined-case-block header (names)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadCom  ( unIn, fileName, 'the combined-case-block header (units)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )

      if ( errStat >= AbortErrLev ) then
         call cleanup()
         return
      end if
      
   if ( DvrData%NumCases < 1 )  then
      call setErrStat( ErrID_Fatal,'Variable "NumCases" must be > 0.' ,errstat,errmsg,routinename)
      call cleanup()
      return
   end if
   
   allocate ( DvrData%Cases(DvrData%NumCases) , STAT=Sttus )
   if ( Sttus /= 0 )  then
      call setErrStat( ErrID_Fatal,'Error allocating memory for the Cases array.',errstat,errmsg,routinename)
      call cleanup()
      return
   end if

   do ICase=1,DvrData%NumCases

      CALL ReadStr( UnIn, fileName, line, 'line', 'Input data for case #'//trim( Int2LStr( ICase ) ), ErrStat2, ErrMsg2, UnEc )
         CALL setErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if ( errStat >= AbortErrLev ) then
            call cleanup()
            return
         end if

      Line = ADJUSTL( Line ) ! remove leading spaces
      
      if (Line(1:1) == '@') then
         Line = Line(2:) ! remove leading character
         call ReadTimeHistoryFile(Line, PriPath, DvrData%Cases(iCase), ErrStat2, ErrMsg2, UnEc )
            call setErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if ( errStat >= AbortErrLev ) then
               call cleanup()
               return
            end if
      else
         
         READ (Line,*,IOSTAT=Sttus)  InpCase ! read whole array (hopefully!)
         
            ! check errors
            CALL CheckIOS ( Sttus, fileName, 'InpCase', NumType, ErrStat2, ErrMsg2 )
               call setErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

            DO Ind=1,NumCols
               CALL CheckRealVar( InpCase(Ind), 'InpCase', ErrStat2, ErrMsg2)
               CALL setErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            END DO
            
            if (ErrStat>=AbortErrLev) then
               call cleanup()
               return
            end if
            
            IF ( UnEc > 0 ) THEN
               WRITE( UnEc, Ec_ReAryFrmt ) TRIM( 'InpCase' ), 'Parameters for Case #'//trim( Int2LStr( ICase ) ), InpCase
            END IF
            
         ! set data
         DvrData%Cases(iCase)%numSteps = ceiling( InpCase( 7) / InpCase( 6) )
         call AllocateCase(DvrData%Cases(iCase), ErrStat2, ErrMsg2) ! needs %numSteps set prior to call
            if (ErrStat2>=AbortErrLev) then
               call SetErrStat( ErrStat2, ErrMsg2 , ErrStat, ErrMsg, RoutineName )
               call Cleanup()
               return
            end if

         DvrData%Cases(iCase)%WndSpeed        = InpCase( 1)
         DvrData%Cases(ICase)%ShearExp        = InpCase( 2)
         DvrData%Cases(ICase)%RotSpeed        = InpCase( 3)*RPM2RPS
         DvrData%Cases(ICase)%Pitch           = InpCase( 4)*D2R
         DvrData%Cases(ICase)%Yaw             = InpCase( 5)*D2R
         DvrData%Cases(iCase)%dT              = InpCase( 6)
        !DvrData%Cases(iCase)%Tmax            = InpCase( 7)
      
         do nt = 1,DvrData%Cases(iCase)%numSteps
            DvrData%Cases(iCase)%time(nt) = (nt-1) * DvrData%Cases(iCase)%dT
         end do
      
      end if
      
   end do ! ICase
   
   call cleanup ( )


   RETURN
contains
   subroutine cleanup()
      if (UnIn>0) close(UnIn)
      if (UnEc>0) close(UnEc)
   end subroutine cleanup
end subroutine Dvr_ReadInputFile
!----------------------------------------------------------------------------------------------------------------------------------
subroutine ReadTimeHistoryFile(FileName, PriPath, CaseData, ErrStat, ErrMsg, UnEc )
   character(*),                  intent(inout) :: FileName
   character(*),                  intent(in   ) :: PriPath
   type(Dvr_Case),                intent(inout) :: CaseData
   integer,                       intent(  out) :: ErrStat           ! returns a non-zero value when an error occurs  
   character(*),                  intent(  out) :: ErrMsg            ! Error message if errStat /= ErrID_None
   integer,                       intent(in   ) :: UnEc

   integer                                      :: UnIn
   integer                                      :: i
   integer, parameter                           :: NumHeaderLines = 2
   
   character(*), parameter                      :: RoutineName = 'AllocateCase'
   integer                                      :: ErrStat2
   character(ErrMsgLen)                         :: ErrMsg2
   
   integer, parameter                           :: NumCols = 6                              ! number of columns to be read from the input file
   real(DbKi)                                   :: InpCase(NumCols)                         ! Temporary array to hold combined-case input parameters. (note that we store in double precision so the time is read correctly)

   
   
   ! Open the input file
   IF ( PathIsRelative( FileName ) ) FileName = TRIM(PriPath)//TRIM(FileName)

   
   call GetNewUnit( UnIn )
   call OpenFInpFile( UnIn, FileName, errStat2, ErrMsg2 )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      if ( errStat >= AbortErrLev ) then
         call cleanup()
         return
      end if

   
   DO I=1,NumHeaderLines
      call ReadCom(UnIn, FileName, 'Header', ErrStat2, ErrMsg2, UnEc)
         call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
         IF (ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
   END DO
   
   
         ! find out how many rows there are to the end of the file
   CaseData%NumSteps   = -1
   ErrStat2 = 0
   DO WHILE ( ErrStat2 == 0 )
      
     CaseData%NumSteps = CaseData%NumSteps + 1
     READ(UnIn, *, IOSTAT=ErrStat2) InpCase(1)
      
   END DO
      
   CALL WrScr( '   Found '//TRIM(Num2LStr(CaseData%NumSteps))//' lines of time-series data.' )
   
   IF (CaseData%NumSteps < 2) THEN
      CALL SetErrStat(ErrID_Fatal, 'The user time-series input file must contain at least 2 rows of time data.', ErrStat, ErrMsg, RoutineName)
      CALL Cleanup()
      RETURN
   END IF
   
   call AllocateCase(CaseData, ErrStat2, ErrMsg2)
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      IF (ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF 
   
      ! now rewind and skip the first few lines. 
   REWIND( UnIn, IOSTAT=ErrStat2 )
      IF (ErrStat2 /= 0_IntKi ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Error rewinding file "'//TRIM(FileName)//'".', ErrStat, ErrMsg, RoutineName)
         CALL Cleanup()
      END IF 

      !IMPORTANT: any changes to the number of lines in the header must be reflected in NumHeaderLines
   DO I=1,NumHeaderLines
      call ReadCom(UnIn, FileName, 'Header', ErrStat2, ErrMsg2, UnEc) ! I'm going to ignore this error because we should have caught any issues the first time we read the file. 
   END DO

   
   DO i=1,CaseData%NumSteps
   
      call ReadAry ( unIn, fileName, InpCase,  NumCols, 'InpCase',  'parameters for Case', errStat2, errMsg2, UnEc )
         call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
            
         if (ErrStat>=AbortErrLev) then
            call Cleanup()
            return
         end if

         CaseData%time(i)        = InpCase( 1)
         CaseData%WndSpeed(i)    = InpCase( 2)
         CaseData%ShearExp(i)    = InpCase( 3)
         CaseData%RotSpeed(i)    = InpCase( 4)*RPM2RPS
         CaseData%Pitch(i)       = InpCase( 5)*D2R
         CaseData%Yaw(i)         = InpCase( 6)*D2R

   END DO
   
   CaseData%dT = CaseData%time(2) - CaseData%time(1)
   
   do i=3,CaseData%NumSteps
      if (.not. EqualRealNos( CaseData%time(i), CaseData%time(i-1) + CaseData%dT ) ) then
         call SetErrStat(ErrID_Fatal,'Time history file must contain time constant deltas in the time channel.', ErrStat, ErrMsg, RoutineName)
         call cleanup()
         return
      end if
   end do
      
   call cleanup()
   
contains
   subroutine cleanup
      close(UnIn)
   end subroutine cleanup
end subroutine ReadTimeHistoryFile
!----------------------------------------------------------------------------------------------------------------------------------
subroutine AllocateCase(CaseData, ErrStat, ErrMsg)
   type(Dvr_Case),                intent(inout) :: CaseData
   integer,                       intent(  out) :: errStat           ! returns a non-zero value when an error occurs  
   character(*),                  intent(  out) :: errMsg            ! Error message if errStat /= ErrID_None
   
   character(*), parameter                      :: routineName = 'AllocateCase'
   integer                                      :: ErrStat2
   character(ErrMsgLen)                         :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   CaseData%numSteps = max(1, CaseData%numSteps)
   
   call AllocAry( CaseData%time,     CaseData%numSteps, 'time',     ErrStat2,ErrMsg2); call setErrStat( ErrStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call AllocAry( CaseData%WndSpeed, CaseData%numSteps, 'WndSpeed', ErrStat2,ErrMsg2); call setErrStat( ErrStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call AllocAry( CaseData%ShearExp, CaseData%numSteps, 'ShearExp', ErrStat2,ErrMsg2); call setErrStat( ErrStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call AllocAry( CaseData%RotSpeed, CaseData%numSteps, 'RotSpeed', ErrStat2,ErrMsg2); call setErrStat( ErrStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call AllocAry( CaseData%Pitch,    CaseData%numSteps, 'Pitch',    ErrStat2,ErrMsg2); call setErrStat( ErrStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call AllocAry( CaseData%Yaw,      CaseData%numSteps, 'Yaw',      ErrStat2,ErrMsg2); call setErrStat( ErrStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )

end subroutine AllocateCase
!----------------------------------------------------------------------------------------------------------------------------------
subroutine ValidateInputs(DvrData, errStat, errMsg)

   type(Dvr_SimData),             intent(inout) :: DvrData           ! intent(out) only so that we can save FmtWidth in DvrData%OutFileData%ActualChanLen
   integer,                       intent(  out) :: errStat           ! returns a non-zero value when an error occurs  
   character(*),                  intent(  out) :: errMsg            ! Error message if errStat /= ErrID_None

   ! local variables:
   integer(intKi)                               :: i
   integer(intKi)                               :: FmtWidth          ! number of characters in string produced by DvrData%OutFmt
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'ValidateInputs'
   
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
      ! Turbine Data:
   if ( DvrData%numBlades < 1 ) call SetErrStat( ErrID_Fatal, "There must be at least 1 blade (numBlades).", ErrStat, ErrMsg, RoutineName)
!   if ( DvrData%numBlades > 3 ) call SetErrStat( ErrID_Fatal, "There can be no more than 3 blades (numBlades).", ErrStat, ErrMsg, RoutineName)
   if ( DvrData%HubRad < 0.0_ReKi .or. EqualRealNos(DvrData%HubRad, 0.0_ReKi) ) call SetErrStat( ErrID_Fatal, "HubRad must be a positive number.", ErrStat, ErrMsg, RoutineName)
   if ( DvrData%HubHt < DvrData%HubRad ) call SetErrStat( ErrID_Fatal, "HubHt must be at least HubRad.", ErrStat, ErrMsg, RoutineName)
   
      
      ! I-O Settings:
      ! Check that DvrData%OutFileData%OutFmt is a valid format specifier and will fit over the column headings
   call ChkRealFmtStr( DvrData%OutFileData%OutFmt, 'OutFmt', FmtWidth, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   !if ( FmtWidth /= ChanLen ) call SetErrStat( ErrID_Warn, 'OutFmt produces a column width of '// &
   !   TRIM(Num2LStr(FmtWidth))//' instead of '//TRIM(Num2LStr(ChanLen))//' characters.', ErrStat, ErrMsg, RoutineName )
   if ( FmtWidth < MinChanLen ) call SetErrStat( ErrID_Warn, 'OutFmt produces a column less than '//trim(num2lstr(MinChanLen))//' characters wide ('// &
      TRIM(Num2LStr(FmtWidth))//'), which may be too small.', ErrStat, ErrMsg, RoutineName )
   DvrData%OutFileData%ActualChanLen = FmtWidth
   
      ! Combined-Case Analysis:
   do i=1,DvrData%NumCases
   
      if (DvrData%Cases(i)%DT < epsilon(0.0_ReKi) ) call SetErrStat(ErrID_Fatal,'dT must be larger than 0 in case '//trim(num2lstr(i))//'.',ErrStat, ErrMsg,RoutineName)
      
   end do
   
   
   
end subroutine ValidateInputs
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_WriteOutputLine(OutFileData, nt, RotAzimuth, output, CaseData, iCase, errStat, errMsg)

   integer(IntKi)         ,  intent(in   )   :: nt                   ! simulation time step (-)
   integer(IntKi)         ,  intent(in   )   :: iCase                ! case # to write to file
   type(Dvr_OutputFile)   ,  intent(in   )   :: OutFileData
   type(Dvr_Case),           intent(in   )   :: CaseData
   real(ReKi)             ,  intent(in   )   :: RotAzimuth           ! Rotor azimuth -- aligned with blade 1 (deg)
   real(ReKi)             ,  intent(in   )   :: output(:)            ! array of requested outputs
   integer(IntKi)         ,  intent(inout)   :: errStat              ! Status of error message
   character(*)           ,  intent(inout)   :: errMsg               ! Error message if ErrStat /= ErrID_None
      
   ! Local variables.

   character(ChanLen)                    :: tmpStr                                    ! temporary string to print the time output as text
   
   errStat = ErrID_None
   errMsg  = ''
  
      ! time
   write( tmpStr, OutFileData%Fmt_t ) CaseData%time(nt)  ! '(F15.4)'
   call WrFileNR( OutFileData%unOutFile, tmpStr(1:OutFileData%ActualChanLen) )
   
   call WrNumAryFileNR ( OutFileData%unOutFile, (/iCase/),  OutFileData%Fmt_i, errStat, errMsg )
   call WrNumAryFileNR ( OutFileData%unOutFile, (/CaseData%WNDSPEED(nt), CaseData%SHEAREXP(nt), RotAzimuth, CaseData%Yaw(nt)*R2D/),  OutFileData%Fmt_a, errStat, errMsg )
   call WrNumAryFileNR ( OutFileData%unOutFile, output,  OutFileData%Fmt_a, errStat, errMsg )
   if ( errStat >= AbortErrLev ) return

     ! write a new line (advance to the next line)
   write (OutFileData%unOutFile,'()')
      
end subroutine Dvr_WriteOutputLine
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_InitializeOutputFile(numBlades, iCase, CaseData, OutFileData, errStat, errMsg)      !TODO:ADP -- how do we tell this routine that we are creating the summary file right now
      integer(IntKi),           intent(in   )   :: numBlades              ! driver data.  neeeded for number of blades
      type(Dvr_OutputFile),     intent(inout)   :: OutFileData 
      
      integer(IntKi)         ,  intent(in   )   :: iCase                ! case number (to write in file description line and use for file name)
      type(Dvr_Case),           intent(in   )   :: CaseData
      
      integer(IntKi)         ,  intent(  out)   :: errStat              ! Status of error message
      character(*)           ,  intent(  out)   :: errMsg               ! Error message if ErrStat /= ErrID_None

         ! locals
      integer(IntKi)                            ::  i      
      integer(IntKi)                            :: numSpaces
      integer(IntKi)                            :: numOuts
      character(ChanLen)                        :: colTxt
      character(ChanLen)                        :: caseTxt
      
      
      
      call GetNewUnit( OutFileData%unOutFile, ErrStat, ErrMsg )
         if ( ErrStat >= AbortErrLev ) then
            OutFileData%unOutFile = -1
            return
         end if
         
      numOuts = size(OutFileData%WriteOutputHdr)

      ! compute the width of the column output
      numSpaces = OutFileData%ActualChanLen ! the size of column produced by OutFmt
      OutFileData%ActualChanLen = max( OutFileData%ActualChanLen, MinChanLen ) ! set this to at least MinChanLen , or the size of the column produced by OutFmt
      do i=1,NumOuts
         OutFileData%ActualChanLen = max(OutFileData%ActualChanLen, LEN_TRIM(OutFileData%WriteOutputHdr(i)))
         OutFileData%ActualChanLen = max(OutFileData%ActualChanLen, LEN_TRIM(OutFileData%WriteOutputUnt(i)))
      end do
      
      ! create format statements for time and the array outputs:
      OutFileData%Fmt_t = '(F'//trim(num2lstr(OutFileData%ActualChanLen))//'.4)'
      OutFileData%Fmt_i = '(I'//trim(num2lstr(OutFileData%ActualChanLen))//')'
      OutFileData%Fmt_a = '"'//OutFileData%delim//'"'//trim(OutFileData%outFmt)      ! format for array elements from individual modules
      numSpaces = OutFileData%ActualChanLen - numSpaces  ! the difference between the size of the headers and what is produced by OutFmt
      if (numSpaces > 0) then
         OutFileData%Fmt_a = trim(OutFileData%Fmt_a)//','//trim(num2lstr(numSpaces))//'x'
      end if
         
      
!      call OpenFOutFile ( OutFileData%unOutFile, trim(outFileData%Root)//'.'//trim(num2lstr(iCase))//'.out', ErrStat, ErrMsg )
      call OpenFOutFile ( OutFileData%unOutFile, trim(outFileData%Root)//'.out', ErrStat, ErrMsg )
         if ( ErrStat >= AbortErrLev ) return
         
      write (OutFileData%unOutFile,'(/,A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//trim( version%Name )
      write (OutFileData%unOutFile,'(1X,A)') trim(GetNVD(OutFileData%AD_ver))
      write (OutFileData%unOutFile,'()' )    !print a blank line
      write (OutFileData%unOutFile,'()' )    !print a blank line
      
      write (OutFileData%unOutFile,'()' )    !print a blank line
              

         !......................................................
         ! Write the names of the output parameters on one line:
         !......................................................

      colTxt = 'Time'
      call WrFileNR ( OutFileData%unOutFile, colTxt(1:OutFileData%ActualChanLen))
      
      colTxt = 'Case'
      call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//colTxt(1:OutFileData%ActualChanLen))
      
      colTxt = 'WindSpeed'
      call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//colTxt(1:OutFileData%ActualChanLen) )
      
      colTxt = 'ShearExp'
      call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//colTxt(1:OutFileData%ActualChanLen) )

      colTxt = 'RotAzimuth'
      call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//colTxt(1:OutFileData%ActualChanLen) )

      colTxt = 'Yaw'
      call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//colTxt(1:OutFileData%ActualChanLen) )

      do i=1,NumOuts
         call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//OutFileData%WriteOutputHdr(i)(1:OutFileData%ActualChanLen) )
      end do ! i

      write (OutFileData%unOutFile,'()')

         !......................................................
         ! Write the units of the output parameters on one line:
         !......................................................
      colTxt = '(s)'
      call WrFileNR ( OutFileData%unOutFile, colTxt(1:OutFileData%ActualChanLen))
      
      colTxt = '(-)'
      call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//colTxt(1:OutFileData%ActualChanLen) )

      colTxt = '(m/s)'
      call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//colTxt(1:OutFileData%ActualChanLen) )
      
      colTxt = '(-)'
      call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//colTxt(1:OutFileData%ActualChanLen) )

      colTxt = '(deg)'
      call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//colTxt(1:OutFileData%ActualChanLen) )

      colTxt = '(deg)'     ! Yaw
      call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//colTxt(1:OutFileData%ActualChanLen) )

      do i=1,NumOuts
         call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//OutFileData%WriteOutputUnt(i)(1:OutFileData%ActualChanLen) )
      end do ! i

      write (OutFileData%unOutFile,'()')

end subroutine Dvr_InitializeOutputFile
!----------------------------------------------------------------------------------------------------------------------------------
end module AeroDyn_Driver_Subs
