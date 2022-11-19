!====================================================================
! NUMODIS - FORTRAN coupler
!--------------------------------------------------------------------
! Authors: Peter Råback,Laurent Dupuy and Javier Gonzalez 2019-2020
!====================================================================

SUBROUTINE NumodisImport( Model,Solver,dt,Transient)

  USE DefUtils
  use global

  IMPLICIT NONE

  TYPE(Solver_t) :: Solver 
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
  !------------------------------------------------------------------
  TYPE(Variable_t), POINTER :: DispVar
  TYPE(Solver_t), POINTER :: DispSolver
  TYPE(Mesh_t), POINTER :: Mesh
  INTEGER, POINTER :: DispSurfPerm(:),TractionSurfPerm(:),StressSurfPerm(:), &
                      StrainSurfPerm(:)
  INTEGER :: i,j,dim,DispSurfNodes,TractionSurfNodes,StressSurfNodes, &
             StrainSurfNodes
  REAL(KIND=dp), POINTER :: Rhs(:)
  LOGICAL :: Visited = .FALSE., Found
!  CHARACTER(LEN=MAX_NAME_LEN) :: FileName
  TYPE(ValueList_t), POINTER :: Params
  !REAL(KIND=dp) :: AppliedDisp(3), AppliedStress(6)
  double precision :: AppliedDisp(3), AppliedStress(6) 
  !------------------------------------------------------------------

  REAL(KIND=dp) ::  ControlLoad
  double precision :: LoadValue



  SAVE DispVar, DispSolver, Mesh, Rhs, dim, &
       DispSurfPerm, DispSurfNodes, StressSurfPerm, StressSurfNodes, &
       TractionSurfPerm, TractionSurfNodes, StrainSurfPerm, StrainSurfNodes


        

  CALL Info('NumodisImport','Setting coupling values on rhs of Displacement')

  Params => GetSolverParams()

  IF(.NOT. Visited ) THEN
     Mesh => GetMesh()

     !-------------------------------------------------------------------
     ! Get the nodes of each BC save then to be used in assembly  
     !-------------------------------------------------------------------

     ! List the nodes that are related to a Numodis coupler Dirichlet nodes
     CALL MakePermUsingMask( Model, Solver, Mesh,'Numodis Dirichlet BC',&
          .FALSE.,DispSurfPerm, DispSurfNodes )
     CALL Info('NumorisImport','Number of coupling nodes for displacement: '&
          //TRIM(I2S(DispSurfNodes)))

     ! List the nodes that are related to a Numodis coupler, Stress nodes 
     CALL MakePermUsingMask( Model, Solver, Mesh,'Numodis Stress BC',&
          .FALSE.,StressSurfPerm, StressSurfNodes )
     CALL Info('NumorisImport','Number of coupling nodes for stress: '&
          //TRIM(I2S(StressSurfNodes)))

     ! List the nodes that are related to a Numodis coupler, Traction nodes 
     CALL MakePermUsingMask( Model, Solver, Mesh,'Numodis Traction BC',&
          .FALSE.,TractionSurfPerm, TractionSurfNodes )
     CALL Info('NumorisImport','Number of coupling nodes for stress: '&
          //TRIM(I2S(TractionSurfNodes)))

     ! List the nodes that are related to a Numodis coupler, Strainrate nodes 
     CALL MakePermUsingMask( Model, Solver, Mesh,'Numodis Strainrate BC',&
          .FALSE.,StrainSurfPerm, StrainSurfNodes )
     CALL Info('NumorisImport','Number of coupling nodes for Indenter : '&
          //TRIM(I2S(StrainSurfNodes)))


     DispVar => VariableGet( Mesh % Variables,'Displacement')
     IF(.NOT. ASSOCIATED(DispVar) ) THEN
        CALL Fatal('NumodisImport','Could not find variable "Displacement"')
     END IF
     dim = DispVar % Dofs

     DispSolver => DispVar % Solver
     IF(.NOT. ASSOCIATED( DispSolver ) ) THEN
        CALL Fatal('NumodisImport','Could not find solver for "Displacement"')
     END IF
     Rhs => DispSolver % Matrix % Rhs

     Visited = .TRUE.
  END IF



  !----------------------------------------------------------------------
  ! set the loading mode; control, and call the assembly functions 
  !----------------------------------------------------------------------

  IF( DispSurfNodes > 0 ) THEN        !!!! Bottom Substrate
!print*, "!!!! Bottom Substrate"
     AppliedDisp = (/0.0d0, 0.0d0, 0.0d0/)
     CALL SetDirichletNumodis(NodIdDir, BCdisplacement, AppliedDisp, DispSurfPerm)
  END IF

  IF( TractionSurfNodes > 0 ) THEN    !!!! Lateral Surfaces
!print *, "!!!! Lateral Surfaces"
     AppliedStress = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
     AppliedStress = AppliedStress * 1.0D6        ! Pa
     !CALL SetStressNumodis(NodIdNeu, BCstress, AppliedStress, StressSurfPerm, StressSurfNodes, 'Stress')
     CALL SetTractionNumodis()
  END IF

  IF( StressSurfNodes > 0 ) THEN      !!!! Constant stress at TOP
!print*, "!!!! Constant stress at TOP"
     LoadValue = SetLoad('StressLoad','ConstantStress')
     AppliedStress = (/0.0d0, 0.0d0, LoadValue, 0.0d0, 0.0d0, 0.0d0/)
     AppliedStress = AppliedStress * 1.0D6        ! Pa
     CALL SetStressNumodis(NodIdNeu, BCstress, AppliedStress, StressSurfPerm, StressSurfNodes, 'Stress')
  END IF


  ControlLoad = GetConstReal(Model%Constants, 'StressControlled', Found)
  IF( Found .AND. StrainSurfNodes > 0 ) THEN
!Print*, 'Strain rate load----- Strees controlled'
     LoadValue = SetLoad('StrainRate','StressControlled')
     AppliedStress = (/0.0d0, 0.0d0, LoadValue, 0.0d0, 0.0d0, 0.0d0/)
     !AppliedStress = AppliedStress * 1.0D6        ! Pa
     AppStressLoad = AppliedStress(3)
     CALL SetStressNumodis(NodIdNeu, BCstress, AppliedStress, StrainSurfPerm, StrainSurfNodes, 'Strain')
  END IF


  ControlLoad = GetConstReal(Model%Constants, 'StrainControlled', Found)
  IF( Found .AND. StrainSurfNodes > 0 ) THEN
!Print*, 'Strain rate load----- Displacement controlled' 
     LoadValue = SetLoad('StrainRate','StrainControlled')
     AppliedDisp = (/0.0d0, 0.0d0, LoadValue/)
     CALL SetDirichletNumodis(NodIdStr, BCstrain, AppliedDisp, StrainSurfPerm)
  END IF
        

  CALL Info('NumodisImport','Finished setting rhs',Level=8)

!!!!!-----------------------------------------------------------------------------------------------------------

CONTAINS
  !====================================================================
  ! Subroutines and functions
  !--------------------------------------------------------------------
  !  SetLoad
  !  SaveLoadData
  !  SetDirichletNumodis
  !  SetStressNumodis
  !  SetTractionNumodis
  !====================================================================



    !----------------------------------------------------------
  ! Loading Feed-Back LOOP
  !----------------------------------------------------------
  FUNCTION SetLoad(app, ContLoad) RESULT(lola)
  !   FEEDBACK LOOP (loading correction plastic + elastic)
  !      Stress corrected == 1
  !      Strain corrected == 2
  ! assuming brick where bottom is fix and load is applied 
  ! only at top
  ! strainrate experiment will have two controll: Stress controlled or displacement controlled
  ! in FEM coupling it directly refer to the kind of BC or LOAD.
  ! contrary to DDD where always app stress is done even if Strainrate is defined.  
  !-------------------------------------------------------
    !IMPLICIT NONE
    REAL*8:: lola
    Real(Kind=dp)::ExtLoad, erate, InitStress, InitStrain, ElmerStrain
    CHARACTER(LEN=10):: app
    CHARACTER(LEN=16):: ContLoad

    REAL(KIND=dp) :: Young_mod, SigmaNC, SigmaC, DispC, DispNC, StrainC, StrainNC
    REAL(KIND=dp) :: plasticstrain(6)

    ! read the young modulus from .sif constants definitions
    Young_mod = GetConstReal(Model%Constants, 'Youngs modulus', Found)
    IF (.NOT. Found) call Fatal('NumodisImport/ SetLoad',& 
                              'Non Young modulus on .sif file')   

    ! get plastic strain
    CALL numodis%exportPlasticStrain(plasticstrain)

      !---------------------------------------------------------------------   
      ! check and assign "CONSTANT STRESS LOAD" if not strainrate experiment
      !---------------------------------------------------------------------  
      ExtLoad = GetConstReal(Model%Constants, 'StressValue', Found)
      IF (Found) THEN
        lola = ExtLoad
        !                 time, Sig33NC,  Sig33C, Str33NC, Str33C, Plstrain(33), Elmerstrain
        call SaveLoadData(istep+1, time, & 
                          0.0_dp, &
                          lola*1e6, &
                          lola*1e6/Young_mod, &
                          0.0_dp, &
                          plasticstrain(3), &
                          abs((InitTop-InitBot)-(Current_top-Current_bottom))/(InitTop-InitBot),&
                          ((InitTop-InitBot)-(Current_top-Current_bottom))/(InitTop-InitBot)*Young_mod/1e6 )
        !                                  Str33NC                                          plasticstrain(33)
      END IF

      !--------------------------------------------------------------------------------
      ! check for "STRAIN RATE LOAD" and error masage if not LOAD defined or misspelled  
      !--------------------------------------------------------------------------------
      IF (.NOT. Found) THEN
        erate = GetConstReal(Model%Constants, 'StrainRate', Found)
        IF (.NOT. Found) CALL Fatal('NumodisImport/ SetLoad',& 
                              'Non StressValue or StrainRate defined...(Check misspelling 1)')  

        !--------------------------------------------------------------------------------
        !  "STRAIN RATE LOAD"  -------->  STRESS CONTROL   
        !--------------------------------------------------------------------------------
        IF (ContLoad == 'StressControlled') THEN  
          ! Set an intial stress for saving time in the elastic regimen
          InitStress = GetConstReal(Model%Constants, 'Initial Stress', Found)
          IF (.NOT. Found) InitStress = 0.0_dp

          ! ! non corrected stress 
           SigmaNC = Young_mod*(erate*time)
          ! ! corrected stress 
          SigmaC = Young_mod*((erate*time)-plasticstrain(3))

          ! ! update output value for Load
          lola = SigmaC+InitStress*1e6

          IF( mod(istep,fsnapshot) == 0 ) THEN ! save the loading data
          !                   step, time, Sig33NC,  Sig33C, Str33NC, Str33C, Plstrain(33), ElmerRealstrain 
            CALL SaveLoadData(istep, time, &
                              SigmaNC/1e6 +InitStress, &
                              SigmaC/1e6 +InitStress, &
                              (erate*time)+(InitStress*1e6/Young_mod), &
                              ((erate*time)-plasticstrain(3))+(InitStress*1e6/Young_mod), &
                              plasticstrain(3), & 
                              ((InitTop-InitBot)-(Current_top-Current_bottom))/(InitTop-InitBot), &
                              ((InitTop-InitBot)-(Current_top-Current_bottom))/(InitTop-InitBot)*Young_mod/1e6 )
          END IF
          
          print*, "Initrial Stress (MPa) -->  ", InitStress
          print*, "******      Step              Time(ns)          PlasticStrain"
          print*, "******   SigmaNC(MPa)        SigmaC(MPa)       ElmerStress(MPa)"
          Print*, "******    StrainNC            StrainC            ElmerStrain" 
          print*,"-------------------------------------------------"
          print*, istep+1,'              ',time, plasticstrain(3)
          print*, Young_mod*(erate * time)/1e6 +InitStress, Young_mod*((erate*time)-plasticstrain(3))/1e6 +InitStress, &
                  ((InitTop-InitBot)-(Current_top-Current_bottom))/(InitTop-InitBot)*Young_mod/1e6
          print*, (erate*time)+(InitStress*1e6/Young_mod), ((erate*time)-plasticstrain(3))+(InitStress*1e6/Young_mod), &
                  ((InitTop-InitBot)-(Current_top-Current_bottom))/(InitTop-InitBot)
          print*,"-------------------------------------------------"


        !--------------------------------------------------------------------------------
        !  "STRAIN RATE LOAD"  -------->  STRAIN CONTROL   
        !--------------------------------------------------------------------------------        
        ELSE IF (ContLoad == 'StrainControlled') THEN  ! change to displacement control
          ! set initial strain for saving time in the elastic regimen
          InitStrain = GetConstReal(Model%Constants, 'Initial Strain', Found)
          IF (.NOT. Found) InitStrain = 0.0_dp

          ! load in displacement control for FEM
          DispNC = (InitStrain+(erate*time))*(InitTop-InitBot)
          DispC = (InitStrain+(erate*time)-plasticstrain(3))*(InitTop-InitBot)
          lola = DispC

          ! save the loading data
          IF( mod(istep,fsnapshot) == 0 ) THEN   
            !                 step, time, Sig33NC,  Sig33C, Str33NC, Str33C, Plstrain(33), ElmerRealstrain                     
            CALL SaveLoadData(istep, time, &
                              Young_mod*(InitStrain+(erate * time))/1e6, &
                              Young_mod*(InitStrain+(erate * time) - plasticstrain(3))/1e6, & 
                              (InitStrain+(erate * time)), &
                              (InitStrain+(erate * time) - plasticstrain(3)), &
                              plasticstrain(3), &
                              ((InitTop-InitBot)-(Current_top-Current_bottom))/(InitTop-InitBot), &
                              ((InitTop-InitBot)-(Current_top-Current_bottom))/(InitTop-InitBot)*Young_mod/1e6 )                 
          END IF
          print*, "Initial Strain  --> ",InitStrain 
          print*, "******      Step              Time(ns)          PlasticStrain"
          print*, "******   SigmaNC(MPa)        SigmaC(MPa)       ElmerStress(MPa)"
          Print*, "******    StrainNC            StrainC            ElmerStrain" 
          print*,"-------------------------------------------------"
          print*, istep+1,'              ',time, plasticstrain(3)
          print*, Young_mod*(InitStrain+(erate * time))/1e6, Young_mod*(InitStrain+(erate* time)-plasticstrain(3))/1e6, &
                  ((InitTop-InitBot)-(Current_top-Current_bottom))/(InitTop-InitBot)*Young_mod/1e6
          print*, InitStrain+(erate * time), InitStrain+((erate * time) - plasticstrain(3)), &
                  ((InitTop-InitBot)-(Current_top-Current_bottom))/(InitTop-InitBot)
          print*,"-------------------------------------------------"

        ELSE
          CALL Fatal('NumodisImport/ SetLoad','Non Controll algorithm defined for StrainRate Load')
        END IF

      END IF
  
  END FUNCTION
  !-------------------------------------------------------


  !----------------------------------------------------------
  !  
  !----------------------------------------------------------
  SUBROUTINE SaveLoadData(par0, par1, par2, par3, par4, par5, par6, par7, par8)
    !IMPLICIT NONE
    INTEGER :: OutUnit, par0 
    REAL(Kind=dp) :: par1, par2, par3, par4, par5,par6, par7, par8

    IF (istep==0) THEN
      OPEN(NEWUNIT=OutUnit, FILE='res/FEMLOAD.txt', action='write', status='replace')
      WRITE(OutUnit,*) "       step         time                    Sigma33NC                 ", & 
                               "Sigma33C                    Str33NC                  Str33C",&
          "                      Pstr                   ElmerStrain             ElmerStress"

      WRITE(OutUnit,*)par0, par1, par2, par3, par4, par5, par6, par7, par8
    ELSE  
      OPEN(NEWUNIT=OutUnit, FILE='res/FEMLOAD.txt', action='write', position='append')
      WRITE(OutUnit,*)par0, par1, par2, par3, par4, par5, par6, par7, par8
    END IF
    CLOSE(OutUnit)

  END SUBROUTINE SaveLoadData
  !----------------------------------------------------------



  !-------------------------------------------------------------------
  ! Set BCs when Dirichlet conditions  
  !-------------------------------------------------------------------
  SUBROUTINE SetDirichletNumodis(ID_node, Num_disp, Ext_disp, PermDisp)

    INTEGER :: i,j,k, m, ID_node(:), PermDisp(:)
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp) :: NodeDisp(3), Ext_disp(3), Num_disp(:)
   
    CALL Info('==> NumodisImport:SetDirichletNumodis','Setting Dirichlet conditions',Level=1)

    A => DispSolver % Matrix

    ! First time these are not allocated as we do our own Dirichlet conditions
    IF(.NOT.ALLOCATED(A % ConstrainedDOF)) THEN
       ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
       A % ConstrainedDOF = .FALSE.
    END IF

    IF(.NOT.ALLOCATED(A % Dvalues)) THEN
       ALLOCATE(A % Dvalues(A % NumberOfRows))
       A % Dvalues = 0._dp
    END IF

    m = 0
    DO i=1,Mesh % NumberOfNodes 
       j = PermDisp(i)
       IF( j == 0 ) CYCLE 

       k = 0  ! is this necessary?
   
       ! reading from passed arrays
       m = m + 1
       k = ID_node( m )

       NodeDisp = Ext_disp - Num_disp(3*m-2:3*m)  ! per node
      
       IF( k /= i ) CALL Fatal('NumodisImport','Node Dirichlet indexing is inconsistent!')

       k = DispVar % Perm(i)            
       A % DValues(dim*(k-1)+1:dim*k) = NodeDisp(1:dim)
       A % ConstrainedDOF(dim*(k-1)+1:dim*k) = .TRUE.
    END DO

  END SUBROUTINE SetDirichletNumodis
  !---------------------------------


  !------------------------------------------------------------
  ! Set BCs when stress components are given 
  !------------------------------------------------------------
  SUBROUTINE SetStressNumodis(ID_node, Num_stress, Ext_stress, PermStress, SurfNodes,reference)

    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t), TARGET :: IP
    INTEGER ::k, e, t, n, m, ElemStart, ElemFin, ID_node(:), PermStress(:), SurfNodes
    INTEGER, POINTER :: Indexes(:)
    REAL(KIND=dp) :: detJ, Weight,Normal(3)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:)
    LOGICAL :: Stat
    TYPE(ValueList_t), POINTER :: BC

    REAL(KIND=dp), ALLOCATABLE :: SurfaceNormals(:,:), SurfaceWeights(:)
    REAL(KIND=dp) :: NodeStress(6),StressT(3,3),NodeForce(3)
    LOGICAL :: WeightsComputed = .FALSE.

    CHARACTER(len=6)::reference
    REAL(KIND=dp):: Num_stress(:), Ext_stress(6)


    SAVE WeightsComputed, SurfaceNormals, SurfaceWeights


    CALL Info('==> NumodisImport:SetStressNumodis','Setting Stress BC', Level=1)


    IF( .NOT. WeightsComputed ) THEN

       ElemStart = Mesh % NumberOfBulkElements + 1
       ElemFin = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

       n = Mesh % MaxElementNodes
       ALLOCATE(Basis(n))

       ALLOCATE(SurfaceNormals(SurfNodes,3),SurfaceWeights(SurfNodes))
       SurfaceNormals = 0.0_dp
       SurfaceWeights = 0.0_dp

       DO e=ElemStart,ElemFin

          Element => Mesh % Elements( e )
          Model % CurrentElement => Element

          BC => GetBC( Element )

          ! here is alowed the function to be used depending of the BC specified
          IF (reference == 'Stress') THEN
            IF(.NOT. GetLogical( BC,'Numodis Stress BC',Stat ) ) CYCLE
          ELSE IF (reference == 'Strain') THEN
            IF(.NOT. GetLogical( BC,'Numodis Strainrate BC',Stat ) ) CYCLE
          ELSE
            CALL FATAL('NumodisImport/SetStressNumodis:','Non BC control specified calling the fucntion to set stress BC')
          END IF

          n = Element % TYPE % NumberOfNodes
          Indexes => Element % NodeIndexes

          CALL GetElementNodes( Nodes, Element )      

          IP = GaussPoints( Element )

          DO t=1,IP % n  

             stat = ElementInfo( Element, Nodes, &
                  IP % U(t), IP % V(t), IP % W(t), detJ, Basis )
             Weight = detJ * IP % s(t)

             Normal = NormalVector( Element, Nodes, IP % U(t), IP % V(t), .TRUE. ) 

             SurfaceWeights(PermStress(Indexes)) = &
                  SurfaceWeights(PermStress(Indexes)) + Weight * Basis(1:n)

             DO k=1,dim
                SurfaceNormals(PermStress(Indexes),k) = &
                     SurfaceNormals(PermStress(Indexes),k) + &
                     Weight * Basis(1:n) * Normal(k)
             END DO
          END DO
       END DO

       DO k=1,dim
          SurfaceNormals(:,k) = SurfaceNormals(:,k) / SurfaceWeights

       END DO

       WeightsComputed = .TRUE.
    END IF

    m = 0
    DO i=1,Mesh % NumberOfNodes 
       j = PermStress(i)
       IF( j == 0 ) CYCLE 
       


       ! Reading from memory array
       m = m + 1
       k = ID_node( m )
       NodeStress = Num_stress(6*m-5:6*m)
       
       ! Adapted to numodis stress array format:
       !  Elmer : XX YY ZZ XY YZ XZ
       ! Numodis: XX YY ZZ XY XZ YZ
       NodeStress = (/NodeStress(1), NodeStress(2), NodeStress(3), NodeStress(4), NodeStress(6), NodeStress(5)/)
       
       ! summing up the BC Applied Stress
       NodeStress = Ext_stress - NodeStress
      

       IF( k /= i ) CALL Fatal('NumodisImport','Node stress indexing is inconsistent!')

       Normal = SurfaceNormals(j,:) 

       StressT(1,1) = NodeStress(1)
       StressT(2,2) = NodeStress(2)
       StressT(3,3) = NodeStress(3)

       StressT(1,2) = NodeStress(4)
       StressT(2,3) = NodeStress(5)
       StressT(1,3) = NodeStress(6)

       StressT(2,1) = StressT(1,2)
       StressT(3,2) = StressT(2,3)
       StressT(3,1) = StressT(1,3)

       NodeForce = MATMUL( StressT, Normal )

       !print*, NodeStress, Normal, SurfaceWeights(j), NodeForce(1:dim)*SurfaceWeights(j)
       !print*, StressT
       !print*, Normal
       !print*, SurfaceWeights(j)
       !print*, NodeForce(1:dim)

       k = DispVar % Perm(i)

       !print*, Rhs(dim*(k-1)+1:dim*k)
       !print*, "--------------------------------------------------------"


       Rhs(dim*(k-1)+1:dim*k) = Rhs(dim*(k-1)+1:dim*k) + &
            SurfaceWeights(j) * NodeForce(1:dim)
    END DO

  END SUBROUTINE SetStressNumodis
  !------------------------------------------------------------


  

  !------------------------------------------------------------
  ! Set traction free stress components 
  !------------------------------------------------------------
  SUBROUTINE SetTractionNumodis()

    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t), TARGET :: IP
    INTEGER ::k, e, t, n, m, ElemStart, ElemFin
    INTEGER, POINTER :: Indexes(:)
    REAL(KIND=dp) :: detJ, Weight,Normal(3)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:)
    LOGICAL :: Stat
    TYPE(ValueList_t), POINTER :: BC

    REAL(KIND=dp), ALLOCATABLE :: SurfaceNormals(:,:), SurfaceWeights(:)
    REAL(KIND=dp) :: NodeStress(6),StressT(3,3),NodeForce(3)
    LOGICAL :: WeightsComputed = .FALSE.

    REAL(KIND=dp):: AppliedStress(6)

    SAVE WeightsComputed, SurfaceNormals, SurfaceWeights

    CALL Info('==> NumodisImport:SetTractionNumodis','Setting Traction Free BC', Level=1)

     ! Applied stress (Load at the surface)
    AppliedStress =(/0.0, 0.0, 0.0, 0.0, 0.0, 0.0/) ! MPa
    AppliedStress = AppliedStress * 1.0D6        ! Pa


    IF( .NOT. WeightsComputed ) THEN


       ElemStart = Mesh % NumberOfBulkElements + 1
       ElemFin = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

       n = Mesh % MaxElementNodes
       ALLOCATE(Basis(n))

       ALLOCATE(SurfaceNormals(TractionSurfNodes,3),SurfaceWeights(TractionSurfNodes))
       SurfaceNormals = 0.0_dp
       SurfaceWeights = 0.0_dp


       DO e=ElemStart,ElemFin

          Element => Mesh % Elements( e )
          Model % CurrentElement => Element

          BC => GetBC( Element )

          IF(.NOT. GetLogical( BC,'Numodis Traction BC',Stat ) ) CYCLE

          n = Element % TYPE % NumberOfNodes
          Indexes => Element % NodeIndexes

          CALL GetElementNodes( Nodes, Element )      

          IP = GaussPoints( Element )

          DO t=1,IP % n  

             stat = ElementInfo( Element, Nodes, &
                  IP % U(t), IP % V(t), IP % W(t), detJ, Basis )
             Weight = detJ * IP % s(t)

             Normal = NormalVector( Element, Nodes, IP % U(t), IP % V(t), .TRUE. ) 

             SurfaceWeights(TractionSurfPerm(Indexes)) = &
                  SurfaceWeights(TractionSurfPerm(Indexes)) + Weight * Basis(1:n)

             DO k=1,dim
                SurfaceNormals(TractionSurfPerm(Indexes),k) = &
                     SurfaceNormals(TractionSurfPerm(Indexes),k) + &
                     Weight * Basis(1:n) * Normal(k)
             END DO
          END DO
       END DO

       DO k=1,dim
          SurfaceNormals(:,k) = SurfaceNormals(:,k) / SurfaceWeights

       END DO

       WeightsComputed = .TRUE.
    END IF

    ! open the ASCII file to read
    !OPEN(NEWUNIT=InUnit, FILE=FileName, STATUS='old', ACTION='read')

    m = 0
    DO i=1,Mesh % NumberOfNodes 
       j = TractionSurfPerm(i)
       IF( j == 0 ) CYCLE 

       ! Reading from memory array
       m = m + 1
       k = NodIdTract( m )
       NodeStress = BCtraction(6*m-5:6*m)
       
       ! Adapted to numodis stress array format:
       !  Elmer : XX YY ZZ XY YZ XZ
       ! Numodis: XX YY ZZ XY XZ YZ
       NodeStress = (/NodeStress(1), NodeStress(2), NodeStress(3), NodeStress(4), NodeStress(6), NodeStress(5)/)
      
       ! summing up the BC Applied Stress
       NodeStress = AppliedStress - NodeStress

       IF( k /= i ) CALL Fatal('NumodisImport','Node indexing is inconsistent!')

       Normal = SurfaceNormals(j,:) 

       StressT(1,1) = NodeStress(1)
       StressT(2,2) = NodeStress(2)
       StressT(3,3) = NodeStress(3)

       StressT(1,2) = NodeStress(4)
       StressT(2,3) = NodeStress(5)
       StressT(1,3) = NodeStress(6)

       StressT(2,1) = StressT(1,2)
       StressT(3,2) = StressT(2,3)
       StressT(3,1) = StressT(1,3)

       NodeForce = MATMUL( StressT, Normal )

       k = DispVar % Perm(i)


       Rhs(dim*(k-1)+1:dim*k) = Rhs(dim*(k-1)+1:dim*k) + &
            SurfaceWeights(j) * NodeForce(1:dim)
    END DO

    ! close the ASCII file
    !CLOSE(InUnit)

    !deallocate(NodIdNeu)
    !deallocate(BCstress)

!!$    DEALLOCATE(basis)
!!$    DEALLOCATE(SurfaceNormals)
!!$    DEALLOCATE(SurfaceWeights)

  END SUBROUTINE SetTractionNumodis
  !------------------------------



END SUBROUTINE NumodisImport



