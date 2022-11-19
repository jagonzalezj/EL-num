!====================================================================
! NUMODIS - FORTRAN coupler
!--------------------------------------------------------------------
! Authors: Peter Råback,Laurent Dupuy and Javier Gonzalez 2019-2020
!====================================================================

SUBROUTINE ExportStress( Model,Solver,dt,TransientSimulation )

  USE DefUtils
  USE ElementDescription
  use global
  use libnumodis
  use libstressexchange

  IMPLICIT NONE

  TYPE(Solver_t) :: Solver 
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  REAL(KIND=dp), ALLOCATABLE :: InterpStress(:)
  LOGICAL::Found

  ! Numodisvariables
  integer :: isnapshot
  double precision :: actualdmax,actualdtime
  double precision :: vmax
  integer :: nGPointsPerSegment = 1
  integer :: nGPoints
  integer :: nGrains
  integer, dimension(:), allocatable :: shift
  real(kind=dp), dimension(:), allocatable :: GPoints, testDE

  ! nucleation variables
  integer :: WaitTime, Pubellator, Point, StepFlag,i, l, option, burguer(3), plane(3), flagAdd, flagRem
  integer :: arrElem = 1, NucStep = 1
  REAL(KIND=dp) :: MgOmain(320,18), MgOopt(96,18), DEinit(320), DEinit1(320),  DEfly(320), SiteWtime(3)=0.0, nKMC(1), dt_KMC(1)
  REAL(KIND=dp) :: center(3), relative(3), StressV, dE, Radii, FEMatpoint(6), coordmin(3), coordmax(3), NucPoitAllStress(320)
  REAL(KIND=dp) :: LocalStressPerPoint(6), MeanStressAtNucPoints, LocCoord(3), MeanStress, choise
  REAL(KIND=dp) :: plasticstrain(6),t_DDD, t_KMC, dt_DDD
  integer :: NucleatedPoints(320,2), PassTrouve, NotNeedPoint,ll


#ifdef USE_MEDCOUPLING
  integer :: nMeshNodes
  character(4) :: cisnapshot
  double precision, pointer  :: meshNodes(:)
  double precision, dimension(:), allocatable :: meshStress
  character(len=255) :: directory
#endif
       

  SAVE t_KMC, t_DDD, dt_DDD, DEinit, DEfly, NucleatedPoints, flagAdd 
  !------------------------------------------------------------------------------
  ! Nucleate dislocation on-fly V3 (Michel)
  !------------------------------------------------------------------------------
  
  IF (istep == 0) THEN
      MgOmain = ReadNuclData(1)
      MgOopt = ReadNuclData(2)
      !WaitTime = GetConstReal(Model%Constants, 'WaitTime', Found)
      NucleatedPoints = 0
      flagAdd = 0
      t_DDD = 0.0D0
      t_KMC = 0.0D0
      dt_DDD = numodis%getDTime()
  END IF

  ! modifications here -----------------------------------------------------------
  ! Get box dimension
  call numodis%exportXminXmax(coordmin,coordmax)

  DO  l = 1, size(DEfly)
      LocCoord = GetRalCoord((/MgOmain(l,1), MgOmain(l,2), MgOmain(l,3)/), coordmin, coordmax)
      LocalStressPerPoint = GetoneFEMStress(LocCoord)
      NucPoitAllStress(l) = abs(LocalStressPerPoint(3))
      if (NucPoitAllStress(l) == 0.0D0) then  ! correction if any stress value is = 0.0
         NucPoitAllStress(l) = dabs(AppStressLoad)/1E6
      end if
  END DO
  NucPoitAllStress = (dabs(NucPoitAllStress)/2)*1e-3
  !print*,NucPoitAllStress
  ! modifications here -----------------------------------------------------------

  ! update step
  !print*,'ISTEP', istep, time, dtime

  ! retrieve the loading stress
  StressV = (dabs(AppStressLoad)/2)*1e-9
  !print*, 'StressV',StressV





  ! Evaluate dE array
  DEinit = EvaldEArray(MgOmain(:,12), MgOmain(:,13), MgOmain(:,14), StressV)* 1.60218e-19
  ! modifications here -----------------------------------------------------------
  !DEinit1 = EvaldEArrayV2(MgOmain(:,12), MgOmain(:,13), MgOmain(:,14), NucPoitAllStress)* 1.60218e-19


  DEfly = DEinit

  !print*, "Sum arrays", SUM(DEfly)


  NotNeedPoint = 0
  ll = 0
  IF (NucleatedPoints(1,1) .NE. 0) THEN
  !print*, 'I M HEEERE 2', NotNeedPoint
   DO WHILE (NotNeedPoint == 0)
   !   print*, 'I M HEEERE'
      ll=ll+1
   !   print*, NucleatedPoints(ll,1)
      DEfly(NucleatedPoints(ll,1))=20
      
      IF (NucleatedPoints(ll+1,1) == 0) THEN
   !      print*, "UPDATED"
        NotNeedPoint = 1
      END IF
   END DO
  ENDIF
  ! modifications here -----------------------------------------------------------

  ! calculate dt_kmc as Michel
  dt_KMC = NuclMath(DEinit,1) * 1E9
  !print*, 'dt_KMC', dt_KMC

  print*, "t_KMC                      t_DDD "
  print*, t_KMC+dt_KMC(1), t_DDD, dt_DDD

  !print*, NucleatedPoints 

  IF (t_KMC+dt_KMC(1) .GT. t_DDD+dt_DDD) THEN
     ! DDD step case (nothing to do)
     t_DDD=t_DDD+dt_DDD  

     ! resset list of nuelcated points after a full stress decrease
     IF (t_KMC+dt_KMC(1) .GT. t_DDD*2)  THEN
         NucleatedPoints = 0
     END IF

     ! test code to handle dislocation that nucleates a dissapear


  ELSE
     ! NUCLEATION CASE

     ! select a random nucleation point from all available  
     !Point = SelectNuclMode(DEinit) !DEinit DEfly! modifications here -----------------


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!***********************************************************************************************
      ! SELCTION OF NUCLEATION POINT accounting for already nucleated sites
      !----------------------------- 
      Print*, "Several sites nucleated..... searching... "

      PassTrouve = 0
 10   DO WHILE (PassTrouve == 0)

         ! select a random nucleation point from all available  
         Point = SelectNuclMode(DEfly)

         ! get its corrdinates
         center = GetRalCoord((/MgOmain(Point,1), MgOmain(Point,2), MgOmain(Point,3)/), coordmin, coordmax)

         IF (flagAdd == 0) then ! case of first nucleation, so.., nucleate it! :)
            NucleatedPoints = addremove(NucleatedPoints, Point, 1, 1)
            PassTrouve = 1
         ELSE
            ! check if slected point is near an already nucleated point
            !----------------------------- ----------------------------
            DO l = 1, size(NucleatedPoints,1)
               IF (NucleatedPoints(l,1) /= 0) THEN

                  ! get the relative coordinates to estimate distance
                  relative = GetRalCoord((/MgOmain(NucleatedPoints(l,1),1), MgOmain(NucleatedPoints(l,1), &
                             2), MgOmain(NucleatedPoints(l,1),3)/), coordmin, coordmax)

      !print*, DistPoints(center,relative)
                  ! Calculate distance
                  IF (DistPoints(center,relative)<30) THEN
      !print*, "NO useful point...try next"
                     !PassTrouve = 1
                     goto 10 ! porque esta muy cerca  :)
                  END IF
               END IF 
            END DO
         ! haaa..., ahora no esta cerca de nadie :) nucleate entonces pendejo!   
         PassTrouve = 1
         exit                
         END IF


      ! update t_KMC
      t_KMC = t_KMC + dt_KMC(1)
      END DO


      ! Check if selected point has double glide system available
      !----------------------------- ----------------------------
      IF (MgOmain(Point,10)  == 0.0D0 ) THEN
            burguer = ([int(MgOmain(point,4)*6),int(MgOmain(point,5)*6),int(MgOmain(point,6)*6)])
            plane = ([int(MgOmain(point,7)),int(MgOmain(point,8)),int(MgOmain(point,9))])
            radii = EvalFunctions(MgOmain(Point,16), MgOmain(Point,17), 0.0D0, StressV, 2)
      ELSE
         Call random_number(choise)
         option = nint(choise) 
         IF (option == 0) THEN
            ! deja loaded  :)
         ELSE
            ! load values form optional repeated GS database Mg0opt.txt
            point =  findloc(MgOopt(:,11), MgOmain(Point, 11))
            burguer = ([int(MgOopt(point,4)*6),int(MgOopt(point,5)*6),int(MgOopt(point,6)*6)])
            plane = ([int(MgOopt(point,7)),int(MgOopt(point,8)),int(MgOopt(point,9))])
            radii = EvalFunctions(MgOopt(Point,16), MgOopt(Point,17), 0.0D0, StressV, 2)    
         END IF   
      END IF
      !!!!!!***********************************************************************************************
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



     ! Get box dimension
     !call numodis%exportXminXmax(coordmin,coordmax)

     ! get its PARAMETERS
     ! center = GetRalCoord((/MgOmain(Point,1), MgOmain(Point,2), MgOmain(Point,3)/), coordmin, coordmax)
     ! burguer = ([int(MgOmain(Point,4)*6),int(MgOmain(Point,5)*6),int(MgOmain(Point,6)*6)])
     ! plane = ([int(MgOmain(Point,7)),int(MgOmain(Point,8)),int(MgOmain(Point,9))])
     ! radii = EvalFunctions(MgOmain(Point,16), MgOmain(Point,17), 0.0D0, StressV, 2)

     print*, "point", point
     print*, center
     print*, burguer
     print*, plane
     print*, radii
     ! Nucleate
     call numodis%nucleateLoop(center,radii, plane, burguer)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!**************************************************************
     flagAdd = flagAdd + 1
      !Turn off the recently nucleated  site
     DEfly(point) = 10.0D0

      ! update the register of nucleated points (NucleatedPoints) to extimate distance nexts 
     NucleatedPoints = addremove(NucleatedPoints, point,  flagAdd, 1)
     !!**************************************************************
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! update t_KMC
     !t_KMC = t_KMC + dt_KMC(1)


  END IF
 
  !call SaveNuclData(istep, t_DDD, t_KMC, dt_KMC(1), dt_DDD, StressV, StressV*2, t_DDD, t_DDD)

  !------------------------------------------------------------------------------
  ! Interpolate the values of Stress at dislocation_gaus points
  !------------------------------------------------------------------------------

  ! get the number of Gauss points
  call numodis%exportNumberDislocationGaussPoints(nGPointsPerSegment,nGPoints,nGrains)

  ! construct the Gauss points containers
  allocate(shift(nGrains))
  allocate(GPoints(3*nGPoints))

  ! get the Gauss points on the dislocations
  call numodis%exportDislocationGaussPoints(nGPointsPerSegment,nGrains,shift,nGpoints,GPoints)

  ! Interpolate FEM stress at NumGPoints
  InterpStress = GetFemStress(GPoints)

  !-------------------------
  ! Continue to Numodis flow
  !-------------------------

  ! reset the nodal forces
  call numodis%ResetNodalForces()

  ! compute mirror forces using Weygang method
  call numodis%computeMirrorForces()

  ! compute the elastic forces
  call numodis%computeInternalForces()

  ! compute the core energy forces
  call numodis%computeCoreForces()

  ! compute forces from finite-elements
  call numodis%importStressFieldOnDislocations(nGPointsPerSegment,shift,InterpStress)

  ! add lattice resistance to nodal forces
  call numodis%computeFrictionForces()

  !------------------------------
  ! compute the nodal velocities
  !------------------------------
  vmax = 0.0
  call numodis%computeVelocities(vmax)

  !---------------------
  ! update current time
  !---------------------
  actualdtime = dtime
  if( vmax > 0.0 ) then
     actualdtime = min( dmax/vmax, dtime )    
  end if
  time = time + actualdtime

  call numodis%setTime(time)
  
  !---------------------------
  ! maximum node displacement
  !---------------------------
  actualdmax = min( dmax, vmax*actualdtime )

  !----------------
  ! move the nodes
  !----------------
  call numodis%moveNodes(actualdtime,actualdmax)

  !---------------------
  ! line discretization
  !---------------------
  call numodis%remesh()

  !--------------------------------------------------------
  ! pass the strain increment to the stress/strain control
  !--------------------------------------------------------
  !call numodis%PlasticStrainIncrement()

  !--------------
  ! clean graphs
  !--------------
  call numodis%cleanGraphs()

  !-----------------------------------------------------------
  ! set new tags on all the nodes and lines of the simulation
  !-----------------------------------------------------------
  call numodis%renumber()

  ! update step
  istep = istep + 1
  print*, "Time and fsnapshot:",time, fsnapshot, mod(istep,fsnapshot)

  !--------------------------------------
  ! save current configuration if needed
  !--------------------------------------
  if( mod(istep,fsnapshot) == 0 ) then        
     isnapshot = istep / fsnapshot        
     !call numodis % Save(isnapshot)
     call numodis % Save(istep)

     
     !--------------
     ! save data
     !--------------
     call numodis % PrintData(istep)

     !--------------
     ! save SIGEPS
     !--------------
     !call numodis % PrintSigeps(istep)

     !---------------------------
     ! save mirror into vtk file
     !---------------------------
     IF (MirrorON) then
      call numodis%printMirrorDislocations2VTK(istep)
     END IF

     !--------------------------
     ! save total stress fields
     !--------------------------
#ifdef USE_MEDCOUPLING
     if( numodis%saveStress() ) then
        
        ! get mesh nodes
        nMeshNodes = stressExchange % getNumberOfNodes() ! get number of nodes
        meshNodes => stressExchange % getCoordinates() ! get pointer to nodes coordinates
!!$        do i=1,nMeshNodes
!!$           print*,"node#",i,meshNodes(3*i-2)," ",meshNodes(3*i-1)," ",meshNodes(3*i)
!!$        end do

        ! allocate memory for stress field
!        allocate(meshStress(6*nMeshNodes)) 
        
        ! compute elmer finite-element stress on the mesh        
        meshStress = GetFemStress(MeshNodes) ! get stress of nodes
        call stressExchange % pushField('finite-element',meshStress ) ! push stress to stressExchange

        ! compute numodis internal stress on the mesh
        meshStress = 0.0D0
        call numodis % exportInternalStressFieldOnNodes(nMeshNodes,meshNodes,meshStress)
        call stressExchange % pushField('internal',meshStress)

        ! compute numodis mirror stress on the mesh
        meshStress = 0.0D0
        call numodis % addMirrorStressesOnNodes(nMeshNodes,meshNodes,meshStress)
        call stressExchange % pushField('mirror',meshStress)

        ! compute total stress
        call stressExchange % addFields('total')

        ! save stress field
        write(cisnapshot,'(I0.4)') isnapshot
        directory = numodis % getSaveDirectory()
        print *,'local directory to save stuff (directory) = ',directory
        call stressExchange % writeVTK(trim(directory)//'/STRESS'//cisnapshot)

        ! flush stress exchange fields and deallocate
        call stressExchange % flushFields()
        deallocate(meshStress)
        
     endif
#endif

  end if

  write(*,*) 'time = ', time, ' (ns) time step = ', actualdtime, ' (ns)'
  write(*,*)

  ! TotalStressBenchmark calling at first timestep remove when not needed anymore!
  !if (time-actualdtime == 0.0 ) then
  !  INQUIRE(FILE="pos.txt", EXIST=THERE)
  !  if (THERE) CALL TotalStressBenchmark()
  !end if


  deallocate(shift)
  deallocate(GPoints)     
  deallocate(InterpStress)



CONTAINS

  !=============================================================
  ! function GetOneFEMStress
  !-------------------------------------------------------------
  !> return stress at a given position
  !-------------------------------------------------------------
  !! 
  !! @param position position
  !! @param FEMstress stress field (NUMODIS format)
  !!
  !! Note: the FEM stress field is exported in NUMODIS format
  !!       meaning xx, yy, zz xy, xz, yz in MPa
  !!
  !!============================================================
  FUNCTION GetOneFEMStress(position) RESULT(FEMStress)

    implicit none

    double precision, dimension(3), intent(in) :: position
    double precision, dimension(6) :: FEMStress

    TYPE(Variable_t), POINTER ::  Stress
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t), POINTER :: Element 
    INTEGER::  k, i, n, m, nNodes
    integer, pointer :: ElementNodeIndices(:) 
    REAL(KIND=dp):: u, v, w
    double precision, allocatable :: NodalStress(:,:)
    REAL(KIND=dp), POINTER :: StressValues(:)
    INTEGER, POINTER :: StressPerm(:)

    ! Inquire the Stress matrix from Elmer
    Stress => VariableGet( Solver % Mesh % Variables, 'Stress' )        
    StressValues => Stress % Values
    StressPerm => Stress % Perm
    
    ! default value, in case we don't find an element containing "position"
    FEMstress = 0.0D0
    
    ! Loop over all elements in the Elmer FEM domain
    DO i=1, Solver % Mesh % NumberOfBulkElements

       ! get the active element
       Element => GetActiveElement(i)

       ! get the number of nodes of the element(Element type dependent)  
       nNodes = GetElementNOFNodes(Element)

       ! get the nodes of the element
       CALL GetElementNodes(ElementNodes, Element, Solver)

       !------------------------------------------------------------------------
       ! Check if respective Gauss point is inside the 3D Element
       ! use BrickInside if Hexahedrom mesh; use TetraInside if tetrahedron mesh
       !------------------------------------------------------------------------
       IF ( BrickInside(ElementNodes%x, ElementNodes%y, ElementNodes%z, position(1), position(2), position(3) ) ) THEN
          !IF ( TetraInside(ElementNodes%x, ElementNodes%y, ElementNodes%z, x, y, z) ) THEN !to use tetrahedron mesh 

          ! Get the indices of the nodes of the element
          ElementNodeIndices => Element % NodeIndexes
          
          ! Get Stress component at each node
          ! size depends of the mesh type 8 to 48 for Hexahedrom
          !print *, "======================================================================="

          allocate(NodalStress(size(ElementNodeIndices),6))
          DO k=1, size(ElementNodeIndices)

             ! index of the current node
             n = StressPerm(ElementNodeIndices(k))

             ! Adapted to numodis stress array format:
             ! Elmer  : XX YY ZZ XY YZ XZ
             ! Numodis: XX YY ZZ XY XZ YZ
             NodalStress(k,1) = StressValues(6*n-5)
             NodalStress(k,2) = StressValues(6*n-4)
             NodalStress(k,3) = StressValues(6*n-3)
             NodalStress(k,4) = StressValues(6*n-2)
             NodalStress(k,5) = StressValues(6*n)
             NodalStress(k,6) = StressValues(6*n-1)

             !print*, NodalStress(6*k-5), NodalStress(6*k-4), NodalStress(6*k-3), NodalStress(6*k-2), &
             !NodalStress(6*k-1), NodalStress(6*k)
          END DO
          !print *, "======================================================================="
          
          ! Convert global to local coordinates
          CALL GlobalToLocal(u,v,w,position(1),position(2),position(3),Element,ElementNodes)

          !------------------------------------------------------------
          ! Here create the on-fly array of stress component one per m
          ! m=1 => Sxx; m=2 => Syy; m=3 => Szz; m=4 => Sxy.......etc :) 
          ! Have to deallocate for every Sxx, Syy, Szz, ....etc
          !------------------------------------------------------------
          
          ! Looping over the components of the variable          
          DO m = 1, 6 ! components of the stress tensor 
             FEMStress(m) = 1.0D-6 * InterpolateInElement( Element, NodalStress(:,m), u, v, w )  ! [Pa] => [MPa]
          END DO

          ! deallocate stuff
          deallocate(NodalStress)
          
          ! exit do loop if the element containing the position was found
          exit
          
       END IF

    END DO

    ! clear memory
    if( ASSOCIATED( ElementNodes % x ) ) then
       deallocate(ElementNodes % x, ElementNodes % y, ElementNodes % z) 
    endif
           
  end FUNCTION GetOneFEMStress
  !-------------------------------------------------------



  !-------------------------------------------------------
  ! create arrays of node indexes at each BC
  ! l = number of nodes
  ! store them in Indx
  !-------------------------------------------------------
  FUNCTION GetFemStress(GPoints) RESULT(InterpStress)
    IMPLICIT NONE
    REAL(KIND=dp), ALLOCATABLE :: InterpStress(:)
    real(kind=dp), dimension(:) :: GPoints
    
    integer::h    
    double precision, dimension(3) :: point
    double precision, dimension(6) :: stress
    
    ! Give size to array to save interpolated stress values: (6 * Gpoints) => 2 = 3 / 6
    ALLOCATE(InterpStress( 2 * SIZE(GPoints) ) )

    !------------------------------------------------------------------------------
    ! Loop over all dislocation_gaus points
    !------------------------------------------------------------------------------
    DO h = 1, SIZE(GPoints)/3

       ! Select coordinates of respective Guss point(h)
       point(1:3) = GPoints(3*h-2:3*h)

       stress = GetOneFEMStress(point)

       InterpStress( 6*(h-1)+1:6*h ) = stress(1:6)

    END DO

  END FUNCTION GetFemStress
  !-------------------------------------------------------




  !-------------------------------------------------------
  ! Read and select the parameter for nucleate dislocations
  ! on fly. 
  !-------------------------------------------------------
  FUNCTION SelectNuclMode(DE) RESULT(Ncase)
    IMPLICIT NONE
    REAL(KIND=dp) :: dt, FlagSum, RSel
    REAL(KIND=dp), ALLOCATABLE :: DN(:), Nu(:), Cint(:), NewDN(:)
    REAL(KIND=dp), INTENT(IN)::DE(:)
!    LOGICAL :: Found
    INTEGER :: h, i, k, j, Ncase
    
    !----------------------------------------
    ! Read the Activ. Energ. & store in DE 
    !----------------------------------------
    !ALLOCATE(DE(0))
    !DE = ReadNuclData(0) * 1.60218e-19 ! read the DE from database (mode = 0)

    !----------------------------------------
    ! Call nucleation equations & get DN
    !---------------------------------------- 
    ALLOCATE(DN(size(DE)), NewDN(size(DE)+1), Cint(size(DE)+2))
    DN = NuclMath(DE,5)


    !----------------------------------------
    ! add the "Do Nothing" into DN --> 1-sum(DN) 
    !----------------------------------------  
    NewDN = [DN, 1-sum(DN)]   

    !----------------------------------------
    ! Make the Markov Chain:  Cint 
    !---------------------------------------- 
    FlagSum = 0.0
    Cint(1) = FlagSum
    DO h = 1, size(NewDN)
      FlagSum = NewDN(h)+FlagSum
      Cint(h+1)=FlagSum
    END DO
    
    !----------------------------------------
    ! generate a random number [0,1] for posible event
    !----------------------------------------
    call random_number(Rsel)
    !print*, 'random_number', Rsel
    !--------------------------------------------------------
    ! Select between "nucleation events" and "Do nothing" cases
    !--------------------------------------------------------
    DO i =1, size(Cint)-1 ! loop on all step of Cint; Nucleate (or do nothing) and break if found

      IF (i == size(Cint)-1) THEN  ! this is the "Do nothing" case [last part of Markov]
          PRINT*, "DO nothing case....."
          Ncase = i
          exit
      END IF

      IF (Rsel<Cint(i+1) .AND. Rsel>Cint(i)) THEN ! this is the "nucleation" case
        ! PRINT*, 'Nucleating event ',i, ' with DN = ', DN(i), "and DE", DE(i)/1.60218e-19
        ! print*, "Nucleating with ",i!, ReadNuclData(i) 
        ! call NucleateEvent(i)
        Ncase = i
        exit
      END IF

    END DO

    DEALLOCATE(DN, NewDN, Cint)


  END FUNCTION SelectNuclMode
  !-------------------------------------------------------




  !-------------------------------------------------------
  ! Reading Jonathan nucleation database V2
  !-------------------------------------------------------
  FUNCTION ReadNuclData(mode) RESULT(salida)
    IMPLICIT NONE
    INTEGER :: OutUnit, InUnit, i, ios, mode, numfiles
    REAL*8 :: x, y, z, bh, bk, bl, ph,pk,pl, rep, ref
    REAL*8 :: Cp1, Cp2, Cp3, CRsq, lm, ln, lRsq
    REAL(kind=dp), ALLOCATABLE:: salida(:,:)

    IF (mode == 1) then
        OPEN(NEWUNIT=InUnit, FILE='MgOmain1.txt', action='read', status='old')
        numfiles = 320
    END IF

    IF (mode == 2) then
        OPEN(NEWUNIT=InUnit, FILE='MgOopt.txt', action='read', status='old')
        numfiles = 98
    END IF

    ALLOCATE(salida(numfiles,18))
    DO i = 1,numfiles
        READ(InUnit,*, IOSTAT = ios) x, y, z, bh, bk, bl, ph,pk,pl, rep, ref, Cp1, Cp2, Cp3, CRsq, lm, ln, lRsq
        salida(i,:) = (/x, y, z, bh, bk, bl, ph,pk,pl, rep, ref, Cp1, Cp2, Cp3, CRsq, lm, ln, lRsq/)
    END DO

    close(InUnit)

  END FUNCTION ReadNuclData
  !-------------------------------------------------------



  !-------------------------------------------------------
  ! Evaluate dE array
  !-------------------------------------------------------
   FUNCTION EvaldEArray(par1,par2,par3, StressValue) result(outArr)
    IMPLICIT NONE
    INTEGER:: mode
    REAL(kind=dp):: StressValue
    REAL(kind=dp), INTENT(IN):: par1(:), par2(:), par3(:)!, StressValue(:)
    REAL(kind=dp), ALLOCATABLE ::outArr(:)

    ALLOCATE(outArr(size(par1)))
    !cooks eval ---> check with matlab
    ! DO i = 1,size(par1)
    !     outArr(i) = par1(i)*(1-(StressValue(i)/par2(i)))**par3(i)
    !     If (isnan(outArr(i))) outArr(i) = par1(i)*(1-(valormedio(StressValue(:))/par2(i)))**par3(i)
    ! END DO

    DO i = 1,size(par1)
        outArr(i) = par1(i)*(1-(StressValue/par2(i)))**par3(i)
    END DO
   END FUNCTION EvaldEArray
  !-------------------------------------------------------


  !-------------------------------------------------------
  ! Evaluate dE array V2
  !-------------------------------------------------------
   FUNCTION EvaldEArrayV2(par1,par2,par3, StressValue) result(outArr)
    IMPLICIT NONE
    INTEGER:: mode
    !REAL(kind=dp):: StressValue
    REAL(kind=dp), INTENT(IN):: par1(:), par2(:), par3(:), StressValue(:)
    REAL(kind=dp), ALLOCATABLE ::outArr(:)

    ALLOCATE(outArr(size(par1)))

     DO i = 1,size(par1)
         outArr(i) = par1(i)*(1-(StressValue(i)/par2(i)))**par3(i)
         If (isnan(outArr(i))) outArr(i) = par1(i)*(1-(valormedio(StressValue(:))/par2(i)))**par3(i)
     END DO

    !DO i = 1,size(par1)
    !    outArr(i) = par1(i)*(1-(StressValue/par2(i)))**par3(i)
    !END DO
   END FUNCTION EvaldEArrayV2
  !-------------------------------------------------------




  !-------------------------------------------------------
  ! Evaluate single dE or Radii depending on mode
  !-------------------------------------------------------
   FUNCTION EvalFunctions(par1,par2,par3, StressValue, mode) result(outArr)
    IMPLICIT NONE
    INTEGER:: mode
    REAL(kind=dp):: StressValue
    REAL(kind=dp):: par1, par2, par3 
    REAL(kind=dp), ALLOCATABLE ::outArr

    !cooks eval ---> check with matlab
    IF (mode == 1) THEN
        outArr = par1*(1-(StressValue/par2))**par3
    END IF

    IF (mode == 2) THEN
        outArr = par1*StressValue + par2
    END IF

   END FUNCTION EvalFunctions
  !-------------------------------------------------------


  !-------------------------------------------------------
  ! Nucleate dislocation with information from Jonathan
  ! Database and Box dimension
  ! With option to random plane 
  !-------------------------------------------------------
  SUBROUTINE NucleateEvent(case)
    IMPLICIT NONE
    REAL(KIND=dp) :: NuclPar(3), FacesArr(5), xmax, xmin, ymax, ymin, zmax, zmin
    REAL(KIND=dp) :: dcut, planecut, radii, xcord, ycord, zcord, NucFace
    INTEGER :: case, i, LevelNumber, StageNumber, burgers(3), plane(3)
    REAL(KIND=dp) :: coordmin(3)
    REAL(KIND=dp) :: coordmax(3)
    REAL(KIND=dp) :: center(3), Mfact

  ! stress shit to delete
    REAL(KIND=dp), dimension(6):: FEMatpoint, SELFatpoint

    !database multiplicative factor (bizzare dcut and planecut stuff)
    Mfact = 1.0
    !--------------------------------------------------------
    ! Read the "nucleation" parameters from Jonatan database
    !--------------------------------------------------------
    !NuclPar = ReadNuclData(case)
    dcut = NuclPar(1)
    planecut = NuclPar(2)
    radii = NuclPar(3)*1
    print*, "NucleateEvent: ReadNuclData case:", case, NuclPar 

    !--------------------------------------------------------
    ! Get the box sizes min(xmin,ymin,zmin), max(xmax,ymax,zmax)
    !--------------------------------------------------------
    call numodis%exportXminXmax(coordmin,coordmax)

    !print*, coordmin
    !print*, coordmax
    !print*, "From Nuclete event: ", ReadNuclData(case)

    !--------------------------------------------------------
    ! Selecte one of the four equivalent lateral faces of the cube
    !--------------------------------------------------------
    !Temporal for half loop
    If (case == 3)THEN
      xcord = coordmax(1)
      ycord = (coordmax(2)-coordmin(2))*(planecut*Mfact)
      ycord = coordmin(2) + ycord 
    ELSE
      xcord = coordmin(1)
      ycord = (coordmax(2)-coordmin(2))*(planecut*Mfact)
      ycord = coordmin(2) + ycord 
    END IF

    !--------------------------------------------------------
    ! The height
    !--------------------------------------------------------
    zcord = (coordmax(3)-coordmin(3))*(dcut * Mfact) 
    zcord = coordmax(3) - zcord

    center(1) = xcord
    center(2) = ycord 
    center(3) = zcord 
    burgers = (/-3,0,3/) ! FCC
    plane = (/1,1,1/)    ! FCC

    !--------------------------------------------------------
    ! test to retieve the stresses at center (delete it)
    !--------------------------------------------------------
    !FEMatpoint = GetOneFEMStress(center)
    !print*, "FEMSTRES", FEMatpoint
    !call numodis % exportInternalStressFieldOnNodes(1,center,SELFatpoint)
    !print*, "SELFSTRES", SELFatpoint

    !print*, "Nucleating at :", center
    !print*, "With radii :", radii, "at plane case :", i
    !print*, "Box lengths :", (coordmax(1)-coordmin(1)), (coordmax(2)-coordmin(2)), (coordmax(3)-coordmin(3))

    !--------------------------------------------------------
    ! Finally Nucleate
    !--------------------------------------------------------
    call numodis%nucleateLoop(center,radii,plane,burgers)

    !call numodis%nucleateLoop(center,length,[1,1,1],[-3,0,3])
    !call numodis%nucleateDislocation((/5000.D0,2000.D0,5000.D0/) , 2000.D0) 
    !call numodis%nucleateDislocation(center,length)

  END SUBROUTINE NucleateEvent
  !----------------------------------------------------------

  


  !----------------------------------------------------------
  ! The equations uses in the nucleation algorithm
  ! Read the constants values for nucleation in the .sif file
  ! Calculate nucleation rate (mode  = 0)
  ! Calculate nucleation probability per cite (mode  = #) 
  !----------------------------------------------------------
  FUNCTION NuclMath(DE, mode) RESULT(MathOut)
    IMPLICIT NONE
    REAL(KIND=dp) :: Temp, TempM, KB, Alpha, Nu_0, dt, Ncalls
    REAL(KIND=dp), ALLOCATABLE :: MathOut(:), Nu(:), DN(:)
    REAL(KIND=dp), INTENT(IN) :: DE(:)
    INTEGER :: mode
    LOGICAL :: Found

    ALLOCATE(NU(SIZE(DE)), DN(SIZE(DE)))

    !----------------------------------------
    ! Reading the constants from case.sif
    !----------------------------------------
     Temp = GetConstReal(Model%Constants, 'Temp', Found)
    TempM = GetConstReal(Model%Constants, 'TempM', Found)
     Nu_0 = GetConstReal(Model%Constants, 'Nu_0', Found)
       KB = GetConstReal(Model%Constants, 'KB', Found)
    alpha = GetConstReal(Model%Constants, 'alpha', Found)

    !----------------------------------------    
    ! EQUATIONS
    !----------------------------------------
    Nu = Size(DE) * Nu_0 * exp(-(DE *(1-Temp/TempM))/(KB*Temp)) 
    !Nu =  Nu_0 * exp(-(DE *(1-Temp/TempM))/(KB*Temp)) 
    dt = 1/(alpha*sum(Nu))  
    DN = Nu * dt 
    Ncalls = dtime/(dt*1.0E9)  ! nKMC, converted dt to (ns)
    !----------------------------------------
    ! Select mode of output
    !----------------------------------------
    IF (mode  == 0) THEN
      ALLOCATE(MathOut(1))
      !MathOut = dt *1e9 ! converted to ns (numodis unit)
      !MathOut(1) = dt
      MathOut(1) = Ncalls
    ELSE IF (mode  == 1) THEN
      ALLOCATE(MathOut(1))
      MathOut(1) = dt

    ELSE
      ALLOCATE(MathOut(SIZE(DE)))
      MathOut = DN
    END IF

  END FUNCTION NuclMath
  !----------------------------------------------------------




  !----------------------------------------------------------
  ! Retrieve the FEM stress 
  ! at a given set of point stored in file pos.txt
  !----------------------------------------------------------
  SUBROUTINE TotalStressBenchmark() 
    INTEGER :: InUnit, OutUnit, total
    REAL(KIND=dp), allocatable :: BenchStress(:), TotalStress(:)

    OPEN(NEWUNIT=InUnit, FILE='pos.txt', STATUS='old', ACTION='read')
    OPEN(NEWUNIT=OutUnit, FILE='TotalStress.txt')

    READ(InUnit,*) total          ! read the numbers of elements

    ALLOCATE(BenchStress(total))
    ALLOCATE(TotalStress(total*2))

    READ(InUnit,*) BenchStress    ! read the rest of values

    ! inquiring to numodis the stress values in MPa
    TotalStress = GetFemStress(BenchStress) 
    TotalStress = TotalStress * 1.0D6
    WRITE(OutUnit,*) TotalStress

    deallocate(TotalStress, BenchStress)

    CLOSE(InUnit)
    CLOSE(OutUnit)
  END SUBROUTINE TOTALSTRESSBENCHMARK
  !-------------------------------------------------------



  !----------------------------------------------------------
  !  Calculate the mean value of an array
  !----------------------------------------------------------
  function valormedio(Inarray) result(outVM)
      implicit none
      REAL(KIND=dp), INTENT(IN) :: Inarray(:) 
      REAL(KIND=dp) :: outVM 
      integer :: i

      outVM = 0.0D0

      do i=1,size(Inarray)
         outVM = outVM + Inarray(i)  
      end do

      outVM = outVM/size(Inarray)

  end function
!-------------------------------------------------------

!----------------------------------------------------------
!  Transform normalize coordinates to local coordinates
!----------------------------------------------------------
function GetRalCoord(Inpoint, cordmin, cordmax) result(outcord)
   implicit none
   REAL(KIND=dp) :: Inpoint(3), cordmin(3), cordmax(3), outcord(3)
  
   outcord(1) = cordmin(1) + (cordmax(1) - cordmin(1))*Inpoint(1) 
   outcord(2) = cordmin(2) + (cordmax(2) - cordmin(2))*Inpoint(2) 
   outcord(3) = cordmin(3) + (cordmax(3) - cordmin(3))*Inpoint(3) 
end function
!-------------------------------------------------------


!----------------------------------------------------------
!  Distance between points
!----------------------------------------------------------
Function DistPoints(point1, point2) result(outdist)
   implicit none
   REAL(KIND=dp) :: outdist, point1(3), point2(3)
   outdist = sqrt( (point2(1)-point1(1))**2 + (point2(2)-point1(2))**2 + (point2(3)-point1(3))**2)
end function
!-------------------------------------------------------


!----------------------------------------------------------
!  find location of value iside a 1D array
!----------------------------------------------------------
function findloc(Inp1DArray, invalue) result(position)
   implicit none
   REAL(kind=dp):: Inp1DArray(:), invalue
   INTEGER :: position

   position = 0
   do i = 1, size(Inp1DArray)
      if (Inp1DArray(i) == invalue) then
         position = i
         exit
      else 
         cycle        
      end if
   enddo
end function findloc


!----------------------------------------------------------
!  
!----------------------------------------------------------
function addremove(MainArray, element, flag, mode) result(outArray)
   implicit none
   INTEGER, intent(in) :: MainArray(:,:)
   INTEGER, allocatable :: MarrayCopy(:,:), outArray(:,:)
   INTEGER :: element, mode, flag

   ! case of size MainArray == 0
   if (size(MainArray) == 0) then
      !print*, "case of size MainArray == 0"
      allocate(MarrayCopy(1,2))
      MarrayCopy(1,1) = Element
      MarrayCopy(1,2) = 0
      outArray = MarrayCopy
      !print*, outArray
      !print*,  "---------------------------"
   endif


   ! include an element
   if (mode == 1) then
      outArray  = MainArray
      outArray(flag,1) = element
      outArray(flag,2) = 0
      !print*, outArray, 'outArray after add'
   endif

   ! remove an element
   if (mode == 2) then
      outArray  = MainArray
      outArray(flag,1) = 0
      outArray(flag,2) = 0
      !print*, outArray, 'outArray after remove'      
   endif


end function addremove



  !----------------------------------------------------------
  !  
  !----------------------------------------------------------
  SUBROUTINE SaveNuclData(par0, par1, par2, par3, par4, par5, par6, par7, par8)
    !IMPLICIT NONE
    INTEGER :: OutUnit, par0 
    REAL(Kind=dp) :: par1, par2, par3, par4, par5,par6, par7, par8

    IF (istep==0) THEN
      OPEN(NEWUNIT=OutUnit, FILE='NUCLEATION.txt', action='write', status='replace')
      WRITE(OutUnit,*) "       step         t_DDD               t_KMC                 ", & 
                               "dt_KMC                    dt_DDD                  ShearStress",&
          "                     NuclStresss                   t_DDD             t_DDD"

      WRITE(OutUnit,*)par0, par1, par2, par3, par4, par5, par6, par7, par8
    ELSE  
      OPEN(NEWUNIT=OutUnit, FILE='NUCLEATION.txt', action='write', position='append')
      WRITE(OutUnit,*)par0, par1, par2, par3, par4, par5, par6, par7, par8
    END IF
    CLOSE(OutUnit)

  END SUBROUTINE SaveNuclData
  !----------------------------------------------------------


END SUBROUTINE ExportStress






!   !------------------------------------------------------------------------------
!   ! Nucleate dislocation on-fly V1
!   !------------------------------------------------------------------------------
!   IF (istep == 0) THEN
! nucleation variables
!  integer :: WaitTime, Pubellator, StepFlag,i
!  integer :: arrElem = 1, NucStep = 1
!  REAL(KIND=dp) :: DEinit(3), DEfly(3), SiteWtime(3)=0.0, nKMC(1)

!   SAVE WaitTime, DEinit, DEfly, nKMC, StepFlag
!   !------------------------------------------------------------------------------
!   ! Nucleate dislocation on-fly V1
!   !------------------------------------------------------------------------------
!   IF (istep == 0) THEN
!      WaitTime = GetConstReal(Model%Constants, 'WaitTime', Found)
!      Print*,WaitTime, "WAITING TIME"
!      DEinit = ReadNuclData(0) * 1.60218e-19 ! read the DE from database 
!      DEfly = DEinit
!      nKMC = NuclMath(DEinit,0)
!   END IF
! 
!   !SiteWtime = UpdateWaitTime()
!   DO i=1,size(DEfly)
!     IF (SiteWtime(i) .GT. 0.0) THEN
!        DEfly(i) = 20 * 1.60218e-19
!        SiteWtime(i) = SiteWtime(i)-1
!        IF (SiteWtime(i) == 0.0) THEN
!           DEfly(i) = DEinit(i)
!        END IF
!     END IF
!   END DO
! 
!   print*,"current step",istep+1
!   print*,'nextcall',WaitTime*NucStep 
!   print*,'nKMC=',nKMC
!   print*,'stepflag=',stepflag
!   print*,'SiteWtime',SiteWtime
!   print*,"defly",DEfly/1.60218e-19
! 
!     ! ! NUCLEATION STUFF FOR TOY PAPER (Javier Version)
!   IF (istep+1 == WaitTime*NucStep .OR. istep+1 == StepFlag) THEN
! 
!       ! Update waiting time calling without stepflag for low T
!       IF (istep+1 == WaitTime*NucStep) THEN 
!         NucStep = NucStep + 1
!       END IF
! 
!       IF (nKMC(1) .LT. 1) THEN  ! Case when T low and nKMC < 1 (call nucleation only one time)
! 
!           ! nucleate 
!           Pubellator = SelectNuclMode(DEfly)
!           ! Update DE
!           DEfly(Pubellator) = 20 * 1.60218e-19
!           ! get new nKMC
!           nKMC = NuclMath(DEfly,0)
!           ! update SiteWtime
!           SiteWtime(Pubellator) = WaitTime
!           ! update the step flag for calling again before waitingtime
!           StepFlag = FLOOR(1/nKMC(1)) + istep+1
! 
!       ELSE  ! Case when T High and nKMC > 1 (call nucleation several times)
!           ! Do STUFFS
!           DO WHILE (nKMC(1) .GT. 1)
!             ! nucleate
!             Pubellator = SelectNuclMode(DEfly)
!             ! update dE
!             DEfly(Pubellator) = 20 * 1.60218e-19
!             ! get new nKMC
!             nKMC = NuclMath(DEfly,0)
! 
!           END DO
!           !reset DE and nKMC
!           DEfly = DEinit
!           nKMC = NuclMath(DEfly,0)
!           ! ! counter of nucleation step
!           ! NucStep = NucStep + 1
!       END IF
! 
!   END IF
!------------------------------------------------------------------------------






!   SAVE WaitTime, NucleatedPoints, DEinit, DEfly , flagAdd!, nKMC, StepFlag
!   !------------------------------------------------------------------------------
!   ! Nucleate dislocation on-fly V2.1
!   !------------------------------------------------------------------------------
  
!   ! Read the databases
!   ! 1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17  18
!   ! X  Y  Z  bh bk bl h  k  l REP  REF P1  P2  P3  R^2  m  n   R^2
!   IF (istep == 0) THEN
!       MgOmain = ReadNuclData(1)
!       MgOopt = ReadNuclData(2)
!       !WaitTime = GetConstReal(Model%Constants, 'WaitTime', Found)
!       NucleatedPoints = 0
!       flagAdd = 0
!   END IF

!   ! Get box dimension
!   call numodis%exportXminXmax(coordmin,coordmax)

!   ! retrieve the loading stress
!   StressV = dabs(AppStressLoad)/2*1e-9

!   ! Evaluate dE array
!   DEinit = EvaldEArray(MgOmain(:,12), MgOmain(:,13), MgOmain(:,14), StressV)* 1.60218e-19
!   DEfly = DEinit

!   ! Calculate nKMC
!   nKMC = NuclMath(DEfly,0)
!   print*, "nKMC = ", nKMC

!   ! CALL NUCLEATION ALGO IF nKMC > 1
!   !---------------------------------
!   IF (nKMC(1) .GT. 1.0) THEN

!       ! get the average Stress at  all the nucleation points
!       DO  l = 1, size(DEfly)
!          LocCoord = GetRalCoord((/MgOmain(l,1), MgOmain(l,3), MgOmain(l,3)/), coordmin, coordmax)
!          LocalStressPerPoint = GetoneFEMStress(LocCoord)
!          NucPoitAllStress(l) = abs(LocalStressPerPoint(3))
!       END DO
!       MeanStress = valormedio(NucPoitAllStress)

!       ! filter the values with larger or lower values of stress (works fine [0.5 -- 0.8])
!       DO  l = 1, size(DEfly)
!          IF (NucPoitAllStress(l) > MeanStress+MeanStress*0.4 .OR. NucPoitAllStress(l) < MeanStress-MeanStress*0.4) THEN
!             DEfly(l) = 20.0
!          END IF
!       END DO

!       ! SELCTION OF NUCLEATION POINT
!       !----------------------------- 
!       PassTrouve = 0
!  10   DO WHILE (PassTrouve == 0)

!          ! select a random nucleation point from all available  
!          Point = SelectNuclMode(DEfly)

!          ! get its corrdinates
!          center = GetRalCoord((/MgOmain(Point,1), MgOmain(Point,2), MgOmain(Point,3)/), coordmin, coordmax)

!          IF (flagAdd == 0) then ! case of first nucleation, so.., nucleate it! :)
!             NucleatedPoints = addremove(NucleatedPoints, Point, 1, 1)
!             PassTrouve = 1
!          ELSE
!             ! check if slected point is near an already nucleated point
!             !----------------------------- ----------------------------
!             DO l = 1, size(NucleatedPoints,1)
!                IF (NucleatedPoints(l,1) /= 0) THEN

!                   ! get the relative coordinates to estimate distance
!                   relative = GetRalCoord((/MgOmain(NucleatedPoints(l,1),1), MgOmain(NucleatedPoints(l,1), &
!                              2), MgOmain(NucleatedPoints(l,1),3)/), coordmin, coordmax)

! print*, DistPoints(center,relative)
!                   ! Calculate distance
!                   IF (DistPoints(center,relative)<10) THEN
! print*, "NO este no,,, epmezamos de nuevo"
!                      !PassTrouve = 1
!                      goto 10 ! porque esta muy cerca  :)
!                   END IF
!                END IF 
!             END DO
!          ! haaa..., ahora no esta cerca de nadie :) nucleate entonces pendejo!   
!          PassTrouve = 1
!          exit                
!          END IF
!       END DO
         
   
!       ! Check if selected point has double glide system available
!       !----------------------------- ----------------------------
!       IF (MgOmain(Point,10)  == 0.0D0 ) THEN
!             burguer = ([int(MgOmain(point,4)*6),int(MgOmain(point,5)*6),int(MgOmain(point,6)*6)])
!             plane = ([int(MgOmain(point,7)),int(MgOmain(point,8)),int(MgOmain(point,9))])
!             radii = EvalFunctions(MgOmain(Point,16), MgOmain(Point,17), 0.0D0, StressV, 2)
!       ELSE
!          Call random_number(choise)
!          option = nint(choise) 
!          IF (option == 0) THEN
!             ! deja loaded  :)
!          ELSE
!             ! load values form optional repeated GS database Mg0opt.txt
!             point =  findloc(MgOopt(:,11), MgOmain(Point, 11))
!             burguer = ([int(MgOopt(point,4)*6),int(MgOopt(point,5)*6),int(MgOopt(point,6)*6)])
!             plane = ([int(MgOopt(point,7)),int(MgOopt(point,8)),int(MgOopt(point,9))])
!             radii = EvalFunctions(MgOopt(Point,16), MgOopt(Point,17), 0.0D0, StressV, 2)    
!          END IF   
!       END IF

!       ! Nucleate
!       call numodis%nucleateLoop(center,radii, plane, burguer)
!       flagAdd = flagAdd + 1

!       !Turn off the recently nculeated  site
!       DEfly(point) = 10.0D0

!       ! update the register of nucleated points (NucleatedPoints) to extimate distance nexts 
!       NucleatedPoints = addremove(NucleatedPoints, point,  flagAdd, 1)

!       !Update nKMC
!       nKMC = NuclMath(DEfly,0)
!   END IF



!   SAVE WaitTime, NucleatedPoints, DEinit, DEfly , flagAdd!, nKMC, StepFlag
!   !------------------------------------------------------------------------------
!   ! Nucleate dislocation on-fly V2.1
!   !------------------------------------------------------------------------------
  
!   ! Read the databases
!   ! 1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17  18
!   ! X  Y  Z  bh bk bl h  k  l REP  REF P1  P2  P3  R^2  m  n   R^2
!   IF (istep == 0) THEN
!       MgOmain = ReadNuclData(1)
!       MgOopt = ReadNuclData(2)
!       WaitTime = GetConstReal(Model%Constants, 'WaitTime', Found)
!       NucleatedPoints = 0
!       flagAdd = 0
!   END IF

!   ! Get box dimension
!   call numodis%exportXminXmax(coordmin,coordmax)

!   ! get the average Stress at  all the nucleation points
!   DO  l = 1, size(DEfly)
!       LocCoord = GetRalCoord((/MgOmain(l,1), MgOmain(l,2), MgOmain(l,3)/), coordmin, coordmax)
!       LocalStressPerPoint = GetoneFEMStress(LocCoord)
!       NucPoitAllStress(l) = abs(LocalStressPerPoint(3))
!       if (NucPoitAllStress(l) == 0.0D0) then  ! correction if any stress value is = 0.0
!          NucPoitAllStress(l) = dabs(AppStressLoad)/1E6
!       end if
!   END DO
!   print*,NucPoitAllStress
  
!   ! retrieve the loading stress
!   StressV = dabs(AppStressLoad)!/2*1e-9
!   print*,StressV

!   ! Evaluate dE array
!   !DEinit = EvaldEArray(MgOmain(:,12), MgOmain(:,13), MgOmain(:,14), (NucPoitAllStress/1000)/2)* 1.60218e-19
!   DEinit = EvaldEArray(MgOmain(:,12), MgOmain(:,13), MgOmain(:,14),(StressV/100)/2) * 1.60218e-19
!   DEfly = DEinit
  

!   ! Calculate nKMC
!   nKMC = NuclMath(DEfly,0)
!   print*, "nKMC = ", nKMC

!   ! CALL NUCLEATION ALGO IF nKMC > 1
!   !---------------------------------
!   !IF (nKMC(1) .GT. 1.0) THEN
!   IF (WaitTime == 0.0) THEN

!       WaitTime = GetConstReal(Model%Constants, 'WaitTime', Found)


!       ! get the average Stress at  all the nucleation points
!       DO  l = 1, size(DEfly)
!          LocCoord = GetRalCoord((/MgOmain(l,1), MgOmain(l,3), MgOmain(l,3)/), coordmin, coordmax)
!          LocalStressPerPoint = GetoneFEMStress(LocCoord)
!          NucPoitAllStress(l) = abs(LocalStressPerPoint(3))
!       END DO
!       MeanStress = valormedio(NucPoitAllStress)

!       ! filter the values with larger or lower values of stress (works fine [0.5 -- 0.8])
!       DO  l = 1, size(DEfly)
!          IF (NucPoitAllStress(l) > MeanStress+MeanStress*0.6 .OR. NucPoitAllStress(l) < MeanStress-MeanStress*0.6) THEN
!             DEfly(l) = 20.0
!          END IF
!       END DO

!       ! SELCTION OF NUCLEATION POINT
!       !----------------------------- 
!       PassTrouve = 0
!  !10   DO WHILE (PassTrouve == 0)
!       do i=1,50
    
!          ! select a random nucleation point from all available  
!   10     Point = SelectNuclMode(DEfly)

!          ! get its corrdinates
!          center = GetRalCoord((/MgOmain(Point,1), MgOmain(Point,2), MgOmain(Point,3)/), coordmin, coordmax)

!          IF (flagAdd == 0) then ! case of first nucleation, so.., nucleate it! :)
!             NucleatedPoints = addremove(NucleatedPoints, Point, 1, 1)
!             PassTrouve = 1
!          ELSE
!             ! check if slected point is near an already nucleated point
!             !----------------------------- ----------------------------
!             DO l = 1, size(NucleatedPoints,1)
!                IF (NucleatedPoints(l,1) /= 0) THEN

!                   ! get the relative coordinates to estimate distance
!                   relative = GetRalCoord((/MgOmain(NucleatedPoints(l,1),1), MgOmain(NucleatedPoints(l,1), &
!                              2), MgOmain(NucleatedPoints(l,1),3)/), coordmin, coordmax)

! print*, DistPoints(center,relative)
!                   ! Calculate distance
!                   IF (DistPoints(center,relative)<10) THEN
! print*, "NO este no,,, epmezamos de nuevo"
!                      !PassTrouve = 1
!                      !goto 10 ! porque esta muy cerca  :)
!                   END IF
!                END IF 
!             END DO
!          ! haaa..., ahora no esta cerca de nadie :) nucleate entonces pendejo!   
!          PassTrouve = 1
!          exit                
!          END IF
!       END DO
         
   
!       ! Check if selected point has double glide system available
!       !----------------------------- ----------------------------
!       IF (MgOmain(Point,10)  == 0.0D0 ) THEN
!             burguer = ([int(MgOmain(point,4)*6),int(MgOmain(point,5)*6),int(MgOmain(point,6)*6)])
!             plane = ([int(MgOmain(point,7)),int(MgOmain(point,8)),int(MgOmain(point,9))])
!             radii = EvalFunctions(MgOmain(Point,16), MgOmain(Point,17), 0.0D0, StressV, 2)
!       ELSE
!          Call random_number(choise)
!          option = nint(choise) 
!          IF (option == 0) THEN
!             ! deja loaded  :)
!          ELSE
!             ! load values form optional repeated GS database Mg0opt.txt
!             point =  findloc(MgOopt(:,11), MgOmain(Point, 11))
!             burguer = ([int(MgOopt(point,4)*6),int(MgOopt(point,5)*6),int(MgOopt(point,6)*6)])
!             plane = ([int(MgOopt(point,7)),int(MgOopt(point,8)),int(MgOopt(point,9))])
!             radii = EvalFunctions(MgOopt(Point,16), MgOopt(Point,17), 0.0D0, StressV, 2)    
!          END IF   
!       END IF

!       ! Nucleate
!       call numodis%nucleateLoop(center,radii, plane, burguer)
!       flagAdd = flagAdd + 1

!       !Turn off the recently nculeated  site
!       DEfly(point) = 10.0D0

!       ! update the register of nucleated points (NucleatedPoints) to extimate distance nexts 
!       NucleatedPoints = addremove(NucleatedPoints, point,  flagAdd, 1)

!       !Update nKMC
!       nKMC = NuclMath(DEfly,0)

!   END IF
!       WaitTime = WaitTime -1









 ! SAVE t_KMC, t_DDD, dt_DDD, DEinit 
 !  !------------------------------------------------------------------------------
 !  ! Nucleate dislocation on-fly V3 (Michel)
 !  !------------------------------------------------------------------------------
  
 !  IF (istep == 0) THEN
 !      MgOmain = ReadNuclData(1)
 !      MgOopt = ReadNuclData(2)
 !      !WaitTime = GetConstReal(Model%Constants, 'WaitTime', Found)
 !      !NucleatedPoints = 0
 !      !flagAdd = 0
 !      t_DDD = 0.0D0
 !      t_KMC = 0.0D0
 !      dt_DDD = numodis%getDTime()
 !  END IF

 !  ! update step
 !  !->
 !  print*,'ISTEP', istep, time, dtime
 !  ! retrieve the loading stress
 !  StressV = (dabs(AppStressLoad)/2)*1e-9
 !  print*, 'StressV',StressV

 !  ! Evaluate dE array
 !  DEinit = EvaldEArray(MgOmain(:,12), MgOmain(:,13), MgOmain(:,14), StressV)* 1.60218e-19

 !  ! calculate dt_kmc as Michel
 !  dt_KMC = NuclMath(DEinit,1) * 1E9
 !  print*, 'dt_KMC', dt_KMC

 !  IF (t_KMC+dt_KMC(1) .GT. t_DDD+dt_DDD) THEN
 !     ! DDD step case (nothing to do)
 !     t_DDD=t_DDD+dt_DDD  
 !  ELSE
 !     ! Nucleation case

 !     ! select a random nucleation point from all available  
 !     Point = SelectNuclMode(DEinit)

 !     ! Get box dimension
 !     call numodis%exportXminXmax(coordmin,coordmax)

 !     ! get its PARAMETERS
 !     center = GetRalCoord((/MgOmain(Point,1), MgOmain(Point,2), MgOmain(Point,3)/), coordmin, coordmax)
 !     burguer = ([int(MgOmain(Point,4)*6),int(MgOmain(Point,5)*6),int(MgOmain(Point,6)*6)])
 !     plane = ([int(MgOmain(Point,7)),int(MgOmain(Point,8)),int(MgOmain(Point,9))])
 !     radii = EvalFunctions(MgOmain(Point,16), MgOmain(Point,17), 0.0D0, StressV, 2)

 !     print*, "point", point
 !     print*, center
 !     print*, burguer
 !     print*, plane
 !     print*, radii
 !     ! Nucleate
 !     call numodis%nucleateLoop(center,radii, plane, burguer)

 !     t_KMC = t_KMC + dt_KMC(1)
 !  END IF
 
 !  !call SaveNuclData(istep, t_DDD, t_KMC, dt_KMC(1), dt_DDD, StressV, StressV*2, t_DDD, t_DDD)