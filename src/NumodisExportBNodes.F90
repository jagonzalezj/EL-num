!====================================================================
! NUMODIS - FORTRAN coupler
!--------------------------------------------------------------------
! Authors: Peter Råback,Laurent Dupuy and Javier Gonzalez 2019-2020
!====================================================================
! The code export boundary nodes coordinates and id to be used for 
! implementing the superposition method coupling between DDD<->FEM.
!
!  Assume a brick geometry (6 faces)
!
!  Traction  = boudary where 0 stress is applied (traction free BC)
!  Neumann   = boundary where non 0 stress is aplied.
!  Dirichlet = boudary where most of the time displacement is aplied
!  Strainrate = boudary where non 0 displcement BC is applied
!
!====================================================================
SUBROUTINE NumodisExportBNodes( Model,Solver,dt,Transient)

  use DefUtils
  use global
  use libnumodis
  use libstressexchange

  IMPLICIT NONE

  TYPE(Solver_t) :: Solver 
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient

  !-------------------------------------------------------
  LOGICAL :: NumodisInitialized = .FALSE. ! true after first step
  
  TYPE(Mesh_t), POINTER :: Mesh
  INTEGER, POINTER :: SurfPerm(:)
  INTEGER :: i,SurfNodes, k, l
  !REAL(KIND=dp), ALLOCATABLE,target :: Dirichlet(:), Neumann(:), &
  !     TractionFree(:), Strainrate(:)
  LOGICAL :: Visited = .FALSE.
  ! numodis variables
  INTEGER :: argc
  character(kind=c_char),dimension(:),allocatable,target :: argv
  character mycommand*16
  LOGICAL :: EXIST_Dirichelt, EXIST_Neumann, EXIST_Traction, EXIST_Strainrate, Found !THERE
  REAL(KIND=dp):: ControlLoad

#ifdef USE_MEDCOUPLING
  double precision :: xmin(3),xmax(3)
  integer :: npts(3)
#endif

  ! Weygand related stuff
  double precision :: MirrorCutoff = 0.D0 ! cutoff for the mirror / Weygand method
  
  ! peter stuff
!  INTEGER, POINTER :: StrPerm(:)
!  REAL(KIND=dp), POINTER :: StrValues(:)
!  INTEGER :: h,n
!  LOGICAL :: Visitado = .FALSE.
!  REAL(KIND=dp), ALLOCATABLE :: NumodisStrValues(:)

!  SAVE :: StrValues

  ! static variables
  SAVE MirrorCutoff
  SAVE Mesh, SurfPerm
  SAVE NumodisInitialized

  EXIST_Dirichelt  = .TRUE.
  EXIST_Neumann = .TRUE.
  EXIST_Traction = .TRUE.
  EXIST_Strainrate = .TRUE.

  CALL Info('NumodisExport','Export Nodes to Numodis')

  CALL CheckAllocations()

  !-----------------------------------------------------------------------
  ! create the numodis object if this function is call for the first time
  !-----------------------------------------------------------------------
  if(.NOT. NumodisInitialized) then
     call NumodisInitialization()
     MirrorCutoff = GetConstReal(Model%Constants, 'Weygand cutoff', MirrorON)
  end if

  ! compute Mirror dislocations once per step (hopefully)
  if (MirrorON) then 
     call numodis % computeMirrorDislocations(MirrorCutoff)
     print*, " computing Mirror dislocations with cutoff = ", MirrorCutoff
  endif        

  !---------------------------------------------------------------------------
  ! Assembly process; detect BC from the .sif file and list the boundary nodes 
  !---------------------------------------------------------------------------
  DO k = 1, 4   ! number of kind of BC surfaces

     IF(.NOT. Visited ) THEN
        Mesh => GetMesh()

        SELECT CASE ( k ) ! List the nodes that are related to a Numodis coupler

        CASE( 1 )   ! Dirichlet (Displ applied constant or rate [fix bottom mostly])
           !print*, "EXPORTBNODES SELCT CASE 1"
           CALL MakePermUsingMask( Model, Solver, Mesh,'Numodis Dirichlet BC',&
                .FALSE.,SurfPerm, SurfNodes )
           IF( SurfNodes == 0 ) THEN
              EXIST_Dirichelt = .FALSE. 
              CALL Info ('NumodisExport','Exporting Dirichlet nodes!')
           END IF

        CASE( 2 )   ! Neumann (Applied Stress constant or rate)
           !print*, "EXPORTBNODES SELCT CASE 2"
           CALL MakePermUsingMask( Model, Solver, Mesh,'Numodis Stress BC',&
                .FALSE.,SurfPerm, SurfNodes )
           IF( SurfNodes == 0 ) THEN
              EXIST_Neumann = .FALSE. 
              CALL Info ('NumodisExport','Exporting Neumann nodes!')
           END IF

        CASE( 3 )   ! Traction_Free (Lateral surfaces)
           !print*, "EXPORTBNODES SELCT CASE 3"
           CALL MakePermUsingMask( Model, Solver, Mesh,'Numodis Traction BC',&
                .FALSE.,SurfPerm, SurfNodes )
           IF( SurfNodes == 0 ) THEN
              EXIST_Traction = .FALSE. 
              CALL Info ('NumodisExport','Exporting Traction_Free nodes!')
           END IF

        CASE( 4 )   ! Strainrate (Displ applied/indenter constant or rate)
           !print*, "EXPORTBNODES SELCT CASE 4"
           CALL MakePermUsingMask( Model, Solver, Mesh,'Numodis Strainrate BC',&
                .FALSE.,SurfPerm, SurfNodes )
           IF( SurfNodes == 0 ) THEN
              EXIST_Strainrate = .FALSE. 
              CALL Info ('NumodisExport','Exporting Strainrate nodes!')
           END IF


        END SELECT

        Visited = .TRUE.

     END IF


     !-----------------------------
     ! Saving the nodes as 1D array
     !-----------------------------
     ! count the number of nodes in the respective BC
     l = COUNT_NODES()

     IF ( k == 1 .AND. EXIST_Dirichelt) then   ! allocate Dirichelet nodes & their indexes
        !print*, "EXPORTBNODES K == 1"
!        ALLOCATE(Dirichlet( l * 3 ) );   ALLOCATE(NodIdDir( l ) )
        Dirichlet = CREATE_NODELISTS(l); NodIdDir = CREATE_NODES_ID(l)
        allocate( BCdisplacement( size( Dirichlet ) ) )
        call numodis % exportDisplacementFieldOnNodes(size( Dirichlet ) / 3, Dirichlet, BCdisplacement)
        Current_bottom = DETECT_MAXMIN(Dirichlet, 3, 0)
        IF (istep == 0) THEN
           InitBot = DETECT_MAXMIN(Dirichlet, 3, 0)
        END IF
        deallocate( Dirichlet )
     END IF

     IF ( k == 2 .AND. EXIST_Neumann) then ! allocate Neumann nodes & their indexes
        !print*, "EXPORTBNODES K == 2"
!        ALLOCATE(Neumann( l * 3 ) );    ALLOCATE(NodIdNeu( l ) )
        Neumann = CREATE_NODELISTS(l);  NodIdNeu = CREATE_NODES_ID(l)
        allocate( BCstress( 2 * size( Neumann ) ) )
        call numodis%exportInternalStressFieldOnNodes(size( Neumann ) / 3, Neumann, BCstress)
        IF (MirrorON) THEN
           Print*, "WeygandMirror (k=2)"
           call numodis % addMirrorStressesOnNodes(size( Neumann ) / 3, Neumann, BCstress)
        END IF
        BCstress = BCstress * 1.0D6 

        !This is the function to save data at boundaries
        ! IF( Visited ) THEN
        !    print*, "Saving SelfStress  at top"
        !    ! Create Permutation from the node id list
        !    n = Mesh % NumberOfNodes
        !    ALLOCATE( StrPerm(n) ) 
        !    StrPerm = 0 
        !    DO h=1,SIZE(NodIdNeu)
        !      StrPerm(NodIdNeu(h)) = h    
        !    END DO
        !    ! Can you see this? 
        !    StrValues => BCstress
        !    CALL VariableAddVector( Mesh % Variables, Mesh, Solver, &
        !   'Str[Str_xx:1 Str_yy:1 Str_zz:1 Str_xy:1 Str_xz:1 Str_yz:1]', &
        !    6, StrValues, StrPerm ) 
        ! END IF

        Current_top = DETECT_MAXMIN(Neumann, 3, 1)
        IF (istep == 0) THEN
           InitTop = DETECT_MAXMIN(Neumann, 3, 1)
        END IF
        deallocate( Neumann )
     END IF

     IF ( k == 3 .AND. EXIST_Traction) then ! allocate Traction free nodes & their indexes
        !print*, "EXPORTBNODES K == 3"
!        ALLOCATE(TractionFree( l * 3 ) );    ALLOCATE(NodIdTract( l ) )
        TractionFree = CREATE_NODELISTS(l);  NodIdTract = CREATE_NODES_ID(l)
        allocate( BCtraction( 2 * size( TractionFree ) ) )
        call numodis % exportInternalStressFieldOnNodes(size( TractionFree ) / 3, TractionFree, BCtraction)
        IF (MirrorON) THEN
           Print*, "WeygandMirror (k=3)"
           call numodis % addMirrorStressesOnNodes(size( TractionFree ) / 3, TractionFree, BCtraction)
        END IF
        BCtraction = BCtraction * 1e6
        deallocate( TractionFree )
     END IF

     IF ( k == 4 .AND. EXIST_Strainrate) then ! allocate StrainRate nodes & their indexes
        !print*, "EXPORTBNODES K == 4"
        ! find the kind of control schema (Stress controlled or displacement controlled)
        ControlLoad = GetConstReal(Model%Constants, 'StrainControlled', Found)


        IF (Found) THEN
           !print*, "EXPORTBNODES K == 4.11111111 Strain controlled"
 !          ALLOCATE(Strainrate( l * 3 ) );    ALLOCATE(NodIdStr( l ) )
           Strainrate = CREATE_NODELISTS(l);  NodIdStr = CREATE_NODES_ID(l)
           allocate( BCstrain( size( Strainrate ) ) )
           call numodis % exportDisplacementFieldOnNodes(size( Strainrate ) / 3, Strainrate, BCstrain)
           Current_top = DETECT_MAXMIN(Strainrate, 3, 1)
           IF (istep == 0) THEN
              InitTop = DETECT_MAXMIN(Strainrate, 3, 1)
           END IF
           deallocate( Strainrate )
        END IF

        IF (.NOT. Found) then
           !print*, "EXPORTBNODES K == 4.222222222 stress Controlled" 
!           ALLOCATE(Neumann( l * 3 ) );    ALLOCATE(NodIdNeu( l ) )
           Neumann = CREATE_NODELISTS(l);  NodIdNeu = CREATE_NODES_ID(l)
           allocate( BCstress( 2 * size( Neumann ) ) )
           call numodis % exportInternalStressFieldOnNodes(size( Neumann ) / 3, Neumann, BCstress)
           IF (MirrorON) THEN
              Print*, "WeygandMirror (k=4)"
              call numodis%addMirrorStressesOnNodes(size( Neumann ) / 3, Neumann, BCstress)
           END IF
           BCstress = BCstress * 1.0D6 
           Current_top = DETECT_MAXMIN(Neumann, 3, 1)
           IF (istep == 0) THEN
              InitTop = DETECT_MAXMIN(Neumann, 3, 1)
           END IF
           deallocate( Neumann )
        END IF
     END IF
     !-------------------------------------

     Visited = .FALSE.

  END DO

  CALL Info('NumodisExport','Finished exporting nodal coordinates',Level=8)

  ! for the benchmarks and test cases remove when not needed
  !if (time == 0.0 ) then
  !  INQUIRE(FILE="pos.txt", EXIST=THERE)
  !  if (THERE) CALL SelfStressBenchmark()
  !end if

CONTAINS

  !-----------------------------------------------------------------------
  ! create the numodis object at the first step and other initializations
  !-----------------------------------------------------------------------  
  SUBROUTINE NumodisInitialization()

    ! set initalized to TRUE
    NumodisInitialized = .TRUE.

    ! create "fake" argument input files
    argc = 3
    allocate(argv(80*argc))
    do i=1,240
       argv(i)=' '
    end do

    mycommand="ElmerSolver"
    do i=1,11
       argv(i) = mycommand(i:i)
    end do
    mycommand="-f"  
    do i=1,2
       argv(80+i) = mycommand(i:i)
    end do
    mycommand="datafiles.xml"
    do i=1,13
       argv(160+i) = mycommand(i:i)
    end do

    numodis = Fnumodis(argc,argv)

    call numodis%initialize()

    fsnapshot = numodis%getSaveFrequency()

    time = numodis%getTime()

    dtime = numodis%getDTime()

    dmax = numodis%getDmax()

    deallocate(argv)

#ifdef USE_MEDCOUPLING     
    call numodis % exportStressGrid(xmin,xmax,npts)
    !print *,"xmin=",xmin, " xmax=",xmax," npts= ",npts
    stressExchange = FstressExchange(xmin,xmax,npts)
#endif

  END SUBROUTINE NumodisInitialization
  !----------------------------------------------------------
  
  !----------------------------------------------------------
  !check allocations of Numodis BC arrays and nodeidexes
  !----------------------------------------------------------
  SUBROUTINE CheckAllocations()
    If (allocated(NodIdDir)) then
       Deallocate(NodIdDir)
       Deallocate(BCdisplacement)
    End If

    If (allocated(NodIdNeu)) then
       Deallocate(NodIdNeu)
       Deallocate(BCstress)
    End If

    If (allocated(NodIdTract)) then
       Deallocate(NodIdTract)
       Deallocate(BCtraction)
    End If

    If (allocated(NodIdStr)) then
       Deallocate(NodIdStr)
       Deallocate(BCstrain)
    End If
  END SUBROUTINE CheckAllocations
  !---------------------------------------------------------


  !----------------------------------------------------------
  ! count the number of nodes in each boundary and store in l
  !----------------------------------------------------------
  FUNCTION COUNT_NODES() RESULT(l)
    IMPLICIT NONE
    INTEGER:: l, i, j

    l= 0  ! to allocate array size
    DO i = 1, Mesh % NumberOfNodes 
       j = SurfPerm( i )
       IF( j == 0 ) CYCLE ! Not a node on the boundary (do not remove!)
       l = l + 1
    END DO
  END FUNCTION COUNT_NODES
  !-------------------------------------------------------


  !-------------------------------------------------------
  ! create arrays of node indexes at each BC
  ! l = number of nodes
  ! store them in Indx
  !-------------------------------------------------------
  FUNCTION CREATE_NODES_ID(l) RESULT(Indx)
    IMPLICIT NONE
    INTEGER:: m, j, i,l
    INTEGER, ALLOCATABLE :: Indx(:)

    IF( l == 0 ) CALL Fatal('NumodisExport','No number of nodes passed to CREATE_NODES_ID function')

    ALLOCATE(Indx( l ) )
    m = 0
    DO i = 1, Mesh % NumberOfNodes 
       j = SurfPerm( i )
       IF( j == 0 ) CYCLE ! Not a node on the boundary (do not remove!)
       m = m + 1
       Indx( m ) = i
    END DO
  END FUNCTION CREATE_NODES_ID
  !-------------------------------------------------------



  !-------------------------------------------------------
  ! create arrays of node positions at each BC
  ! l = number of nodes
  ! store them in Poss
  !-------------------------------------------------------
  FUNCTION CREATE_NODELISTS(l) RESULT(Poss)
    IMPLICIT NONE
    INTEGER:: m, j, i, l
    REAL(KIND=dp), ALLOCATABLE,target :: Poss(:)

    IF( l == 0 ) CALL Fatal('NumodisExport','No number of nodes passed to CREATE_NODELISTS function')

    ALLOCATE(Poss( l * 3 ) )
    m = 0
    DO i = 1, Mesh % NumberOfNodes 
       j = SurfPerm( i )
       IF( j == 0 ) CYCLE ! Not a node on the boundary (do not remove!)
       m = m + 1
       Poss( m * 3 - 2 ) = Mesh % Nodes % x(i)
       Poss( m * 3 - 1 ) = Mesh % Nodes % y(i)
       Poss( m * 3 ) = Mesh % Nodes % z(i)
    END DO
  END FUNCTION CREATE_NODELISTS
  !-------------------------------------------------------



  !---------------------------------------------------------
  ! Fucntion to detect the maximum and minimum  of the box 
  ! dimention (to be used for feedback loop)
  !   - Impvector: List of nodal coordinates
  !   - Corrd:  X, Y, or Z  ===>   (1, 2 or 3)
  !   - MaxOrMin:  "Max"==1   otherwise "Min"==0 
  !---------------------------------------------------------
  FUNCTION DETECT_MAXMIN(Impvector, Coord, MaxOrMin) RESULT(res)
    IMPLICIT NONE
    INTEGER :: Coord, MaxOrMin
    REAL(KIND=dp) :: res, Impvector(:)
    REAL(KIND=dp), ALLOCATABLE :: vect(:)

    ALLOCATE(vect(size(Impvector)/3))

    ! insolate the respective component

    Do i = 1, size(Impvector)/3
       vect(i) = Impvector((i*3) - (3-Coord) ) 
    end do

    ! get the max (1) or min (0)
    IF (MaxOrMin == 1) THEN
       res = maxval(vect)
    ELSE IF (MaxOrMin == 0) THEN
       res = minval(vect)
    ELSE
       CALL Fatal('DETECT_MAXMIN','Non Max or Min defined..')
    END IF

    DEALLOCATE(vect)
    
  END FUNCTION DETECT_MAXMIN
  !-------------------------------------------------------



  !----------------------------------------------------------
  ! Retrive the self-stress of the dislocation microsturctured
  ! at a given set of point stored in file pos.txt
  !----------------------------------------------------------
  SUBROUTINE SelfStressBenchmark() 
    INTEGER :: InUnit, OutUnit1, OutUnit2, total
    REAL(KIND=dp), allocatable :: NodePos(:), SelfStress(:), WeygandStress(:)


    print *, "SELF STRESS BENCHMARK CALLED"


    ! file with node positions
    OPEN(NEWUNIT=InUnit, FILE='pos.txt', STATUS='old', ACTION='read')

    !file to save the dislocations self stress
    OPEN(NEWUNIT=OutUnit1, FILE='SelfStress.txt')

    !file to save the dislocations self stress
    OPEN(NEWUNIT=OutUnit2, FILE='WeygandStress.txt')

    READ(InUnit,*) total          ! read the numbers of elements

    ! alocations 
    ALLOCATE(NodePos(total))   ! node position array
    ALLOCATE(SelfStress(total*2))  ! self stress array
    ALLOCATE(WeygandStress(total*2))  ! Weygand stress array

    READ(InUnit,*) NodePos    ! read the rest of values

    ! inquiring to numodis the Self stress values
    call numodis % exportInternalStressFieldOnNodes(size(NodePos)/3, NodePos, SelfStress)
    WRITE(OutUnit1,*) SelfStress

    ! inquiring to numodis the Weygand stress values
    call numodis % exportInternalStressFieldOnNodes(size(NodePos)/3, NodePos, WeygandStress)
    call numodis % addMirrorStressesOnNodes(size(NodePos)/3, NodePos, WeygandStress)
    WRITE(OutUnit2,*) WeygandStress

    ! dealocations
    deallocate(SelfStress, WeygandStress, NodePos)

    CLOSE(InUnit)
    CLOSE(OutUnit1)
    CLOSE(OutUnit2)

  END SUBROUTINE SelfStressBenchmark
  !-------------------------------------------------------


END SUBROUTINE NumodisExportBnodes


