module libnumodis

  use iso_c_binding


  private

  public :: Fnumodis

  include "numodis_cc2fortran.bind"


  type Fnumodis

     private

     type(c_ptr) :: ptr

   contains

     final :: destroy_numodis

     procedure :: initialize => numodis_initialize
     
!     procedure :: getNstep => numodis_getNstep

     procedure :: saveStress => numodis_saveStress

     procedure :: getSaveFrequency => numodis_getSaveFrequency

     procedure :: getSaveDirectory => numodis_getSaveDirectory

     procedure :: getTime => numodis_getTime

     procedure :: getDmax => numodis_getDmax

     procedure :: getDTime => numodis_getDTime

     procedure :: setTime => numodis_setTime

     procedure :: compute => numodis_compute

     procedure :: Save => numodis_Save

     procedure :: computeAppliedForces => numodis_computeAppliedForces

     procedure :: computeCoreForces => numodis_computeCoreForces

     procedure :: computeFrictionForces => numodis_computeFrictionForces

     procedure :: computeVelocities => numodis_computeVelocities

     procedure :: cleanGraphs => numodis_cleanGraphs

     procedure :: renumber => numodis_renumber

     procedure :: computeInternalForces => numodis_computeInternalForces

     procedure :: moveNodes => numodis_moveNodes

     procedure :: remesh => numodis_remesh

!     procedure :: oPlasticStrainIncrement => numodis_aPlasticStrainIncrement

     procedure :: ResetNodalForces => numodis_ResetNodalForces

     procedure :: exportInternalStressFieldOnNodes => numodis_exportInternalStressFieldOnNodes

     procedure :: exportDisplacementFieldOnNodes => numodis_exportDisplacementFieldOnNodes

     procedure :: exportNumberDislocationGaussPoints => numodis_exportNumberDislocationGaussPoints

     procedure :: exportDislocationGaussPoints => numodis_exportDislocationGaussPoints

     procedure :: importStressFieldOnDislocations => numodis_importStressFieldOnDislocations

!!$     procedure :: OpenSIGEPS => numodis_OpenSIGEPS
!!$
!!$     procedure :: PrintSIGEPS => numodis_PrintSIGEPS
!!$
!!$     procedure :: CloseSIGEPS => numodis_CloseSIGEPS

     procedure :: OpenData => numodis_OpenData

     procedure :: PrintData => numodis_PrintData

     procedure :: CloseData => numodis_CloseData

     procedure :: nucleateDislocation => numodis_nucleateDislocation

     procedure :: nucleateLoop => numodis_nucleateLoop

     procedure :: exportXminXmax => numodis_exportXminXmax

     procedure :: exportStressGrid => numodis_exportStressGrid

     procedure :: exportPlasticStrain => numodis_exportPlasticStrain

     procedure :: computeMirrorDislocations => numodis_computeMirrorDislocations

     procedure :: printMirrorDislocations2VTK => numodis_printMirrorDislocations2VTK

     procedure :: computeMirrorForces => numodis_computeMirrorForces

     procedure :: addMirrorStressesOnNodes => numodis_addMirrorStressesOnNodes
     
  end type Fnumodis


  interface Fnumodis
     procedure init_numodis
  end interface Fnumodis



contains



  function init_numodis(argc,argv)
    implicit none

    type(Fnumodis) :: init_numodis
    integer, intent(in) :: argc
    character(kind=c_char), dimension(:), allocatable, target :: argv

    type(c_ptr), dimension(argc) :: argv_to_c
    integer :: ii

    do ii = 1, argc
       argv_to_c(ii) = c_loc(argv((ii-1)*80+1))
    end do
    init_numodis%ptr = init_numodis_cc(argc,argv_to_c)

  end function init_numodis

  subroutine destroy_numodis(this)
    implicit none
    type(Fnumodis) :: this
    call destroy_numodis_cc(this%ptr)
  end subroutine destroy_numodis

  subroutine numodis_initialize(this)
    implicit none
    class(Fnumodis) :: this
    call numodis_initialize_cc(this%ptr)
  end subroutine numodis_initialize
    
!!$  integer function numodis_getNstep(this)
!!$    implicit none
!!$    class(Fnumodis) :: this
!!$    numodis_getNstep = numodis_getNstep_cc(this%ptr)
!!$  end function numodis_getNstep

  logical function numodis_saveStress(this)
    implicit none
    class(Fnumodis) :: this
    numodis_saveStress = numodis_saveStress_cc(this%ptr)
  end function numodis_saveStress
  
  integer function numodis_getSaveFrequency(this)
    implicit none
    class(Fnumodis) :: this
    numodis_getSaveFrequency = numodis_getSaveFrequency_cc(this%ptr)
  end function numodis_getSaveFrequency

  function numodis_getSaveDirectory(this) result(directory)
    implicit none
    class(Fnumodis) :: this
    character(len=255) :: directory
    character(len=255) :: directory1
    integer :: fsize
    call numodis_getSaveDirectory_cc(this%ptr,fsize,directory1)
    directory = directory1(1:fsize)
  end function numodis_getSaveDirectory
  
  double precision function numodis_getTime(this)
    implicit none
    class(Fnumodis) :: this
    numodis_getTime = numodis_getTime_cc(this%ptr)
  end function numodis_getTime

  double precision function numodis_getDmax(this)
    implicit none
    class(Fnumodis) :: this
    numodis_getDmax = numodis_getDmax_cc(this%ptr)
  end function numodis_getDmax

  double precision function numodis_getDTime(this)
    implicit none
    class(Fnumodis) :: this
    numodis_getDTime = numodis_getDTime_cc(this%ptr)
  end function numodis_getDTime

  subroutine numodis_setTime(this,time)
    implicit none
    class(Fnumodis) :: this
    double precision, intent(in) :: time
    call numodis_setTime_cc(this%ptr,time)
  end subroutine numodis_setTime

  subroutine numodis_compute(this,step)
    implicit none
    class(Fnumodis) :: this
    integer, intent(in) :: step
    call numodis_compute_cc(this%ptr,step)
  end subroutine numodis_compute

  subroutine numodis_Save(this, insave)
    implicit none
    class(Fnumodis) :: this
    integer, intent(in) :: insave
    call numodis_Save_cc(this%ptr,insave)
  end subroutine numodis_Save

  subroutine numodis_computeAppliedForces(this, step)
    implicit none
    class(Fnumodis) :: this
    integer, intent(in) :: step
    call numodis_computeAppliedForces_cc(this%ptr,step)
  end subroutine numodis_computeAppliedForces

  subroutine numodis_computeCoreForces(this)
    implicit none
    class(Fnumodis) :: this
    call numodis_computeCoreForces_cc(this%ptr)
  end subroutine numodis_computeCoreForces

  subroutine numodis_computeFrictionForces(this)
    implicit none
    class(Fnumodis) :: this
    call numodis_computeFrictionForces_cc(this%ptr)
  end subroutine numodis_computeFrictionForces

  subroutine numodis_computeVelocities(this, vmax)
    implicit none
    class(Fnumodis) :: this
    real*8, intent(out) :: vmax
    call numodis_computeVelocities_cc(this%ptr,vmax)
  end subroutine numodis_computeVelocities

  subroutine numodis_cleanGraphs(this)
    implicit none
    class(Fnumodis) :: this
    call numodis_cleanGraphs_cc(this%ptr)
  end subroutine numodis_cleanGraphs

  subroutine numodis_renumber(this)
    implicit none
    class(Fnumodis) :: this
    call numodis_renumber_cc(this%ptr)
  end subroutine numodis_renumber

  subroutine numodis_computeInternalForces(this)
    implicit none
    class(Fnumodis) :: this
    call numodis_computeInternalForces_cc(this%ptr)
  end subroutine numodis_computeInternalForces

  subroutine numodis_moveNodes(this,dtime,dmax)
    implicit none
    class(Fnumodis) :: this
    real*8, intent(in) :: dtime
    real*8, intent(in) :: dmax
    call numodis_moveNodes_cc(this%ptr,dtime,dmax)
  end subroutine numodis_moveNodes

  subroutine numodis_remesh(this)
    implicit none
    class(Fnumodis) :: this
    call numodis_remesh_cc(this%ptr)
  end subroutine numodis_remesh

!  subroutine numodis_aPlasticStrainIncrement(this)
!    implicit none
!    class(Fnumodis) :: this
!    call numodis_assignPlasticStrainInc_cc(this%ptr)
!  end subroutine numodis_aPlasticStrainIncrement  

  subroutine numodis_ResetNodalForces(this)
    implicit none
    class(Fnumodis) :: this
    call numodis_ResetNodalForces_cc(this%ptr)
  end subroutine numodis_ResetNodalForces

  subroutine numodis_exportInternalStressFieldOnNodes(this, npts, pts, stress)
    implicit none
    class(Fnumodis) :: this
    integer, value, intent(in) :: npts
    double precision, dimension(0:), intent(in) :: pts
    double precision, dimension(:), intent(out) :: stress
    call numodis_exportInternalStressFieldOnNodes_cc(this%ptr, npts, pts, stress)
  end subroutine numodis_exportInternalStressFieldOnNodes

  subroutine numodis_exportDisplacementFieldOnNodes(this, npts, pts, displ)
    implicit none
    class(Fnumodis) :: this
    integer, value, intent(in) :: npts
    double precision, dimension(:), intent(in) :: pts
    double precision, dimension(:), intent(out) :: displ
    call numodis_exportDisplacementFieldOnNodes_cc(this%ptr, npts, pts, displ)
  end subroutine numodis_exportDisplacementFieldOnNodes

  subroutine numodis_exportNumberDislocationGaussPoints(this, nGauss, npts, ngrains)
    implicit none
    class(Fnumodis) :: this
    integer, intent(in) :: nGauss
    integer, intent(out) :: npts
    integer, intent(out) :: ngrains
    call numodis_exportNumberDislocationGaussPoints_cc(this%ptr, nGauss, npts, ngrains)
  end subroutine numodis_exportNumberDislocationGaussPoints
  
  subroutine numodis_exportDislocationGaussPoints(this, nGauss, nshift, shift, npts, pts)
    implicit none
    class(Fnumodis) :: this
    integer, intent(in) :: nGauss
    integer, intent(in) :: nshift
    integer, dimension(:), intent(out) :: shift
    integer, intent(in) :: npts
    double precision, dimension(:), intent(out) :: pts
    call numodis_exportDislocationGaussPoints_cc(this%ptr, nGauss, nshift, shift, npts, pts) 
  end subroutine numodis_exportDislocationGaussPoints

  subroutine numodis_importStressFieldOnDislocations(this, nGauss, shift, stress)
    implicit none
    class(Fnumodis) :: this
    integer, value, intent(in) :: nGauss
    integer, intent(in) :: shift(*)
    double precision, intent(in) :: stress(*)
    call numodis_importStressFieldOnDislocations_cc(this%ptr, nGauss, shift, stress)
  end subroutine numodis_importStressFieldOnDislocations

!!$  subroutine numodis_OpenSIGEPS(this)
!!$    implicit none
!!$    class(Fnumodis) :: this
!!$    call numodis_OpenSIGEPS_cc(this%ptr)
!!$  end subroutine numodis_OpenSIGEPS
!!$
!!$  subroutine numodis_PrintSIGEPS(this,step)
!!$    implicit none
!!$    class(Fnumodis) :: this
!!$    integer, intent(in) :: step
!!$    call numodis_PrintSIGEPS_cc(this%ptr,step)
!!$  end subroutine numodis_PrintSIGEPS
!!$
!!$  subroutine numodis_CloseSIGEPS(this)
!!$    implicit none
!!$    class(Fnumodis) :: this
!!$    call numodis_CloseSIGEPS_cc(this%ptr)
!!$  end subroutine numodis_CloseSIGEPS

  subroutine numodis_OpenData(this)
    implicit none
    class(Fnumodis) :: this
    call numodis_OpenData_cc(this%ptr)
  end subroutine numodis_OpenData

  subroutine numodis_PrintData(this,step)
    implicit none
    class(Fnumodis) :: this
    integer, intent(in) :: step
    call numodis_PrintData_cc(this%ptr,step)
  end subroutine numodis_PrintData

  subroutine numodis_CloseData(this)
    implicit none
    class(Fnumodis) :: this
    call numodis_CloseData_cc(this%ptr)
  end subroutine numodis_CloseData

  subroutine numodis_nucleateDislocation(this, center, length )
    implicit none
    class(Fnumodis) :: this
    double precision, intent(in) :: center(*)
    double precision, intent(in) :: length
    write(*,*) 'From numodis_module.f90::numodis_nucleateDislocation: length = ',length
    call numodis_nucleateDislocation_cc(this%ptr,center,length)
  end subroutine numodis_nucleateDislocation

  subroutine numodis_nucleateLoop(this, center, radius, plane, burgers )
    implicit none
    class(Fnumodis) :: this
    double precision, intent(in) :: center(*)
    double precision, intent(in) :: radius
    integer, intent(in)          :: plane(*)
    integer, intent(in)          :: burgers(*)
    write(*,*) 'From numodis_module.f90::numodis_nucleateLoop: radius = ',radius
    call numodis_nucleateLoop_cc(this%ptr,center,radius,plane,burgers)
  end subroutine numodis_nucleateLoop

  subroutine numodis_exportXminXmax(this,xmin,xmax)
    implicit none
    class(Fnumodis) :: this
    double precision, intent(out) :: xmin(*)
    double precision, intent(out) :: xmax(*)
    call numodis_exportXminXmax_cc(this%ptr,xmin,xmax)
  end subroutine numodis_exportXminXmax

  subroutine numodis_exportStressGrid(this,xmin,xmax,npts)
    implicit none
    class(Fnumodis) :: this
    double precision, intent(out) :: xmin(*)
    double precision, intent(out) :: xmax(*)
    integer, intent(out) :: npts(*)
    call numodis_exportStressGrid_cc(this%ptr,xmin,xmax,npts)
  end subroutine numodis_exportStressGrid

  subroutine numodis_exportPlasticStrain(this,strain)
    implicit none
    class(Fnumodis) :: this
    double precision, intent(out) :: strain(*)
    call numodis_exportPlasticStrain_cc(this%ptr,strain)
  end subroutine numodis_exportPlasticStrain

  subroutine numodis_computeMirrorDislocations(this, cutoff)
    implicit none
    class(Fnumodis) :: this
    double precision, intent(in) :: cutoff
    call numodis_computeMirrorDislocations_cc(this%ptr,cutoff)
  end subroutine numodis_computeMirrorDislocations

  subroutine numodis_printMirrorDislocations2VTK(this, insave)
    implicit none
    class(Fnumodis) :: this
    integer, intent(in) :: insave
    call numodis_printMirrorDislocations2VTK_cc(this%ptr,insave)
  end subroutine numodis_printMirrorDislocations2VTK

  subroutine numodis_computeMirrorForces(this)
    implicit none
    class(Fnumodis) :: this
    call numodis_computeMirrorForces_cc(this%ptr)
  end subroutine numodis_computeMirrorForces

  subroutine numodis_addMirrorStressesOnNodes(this, npts, pts, stress)
    implicit none
    class(Fnumodis) :: this
    integer, value, intent(in) :: npts
    double precision, dimension(0:), intent(in) :: pts
    double precision, dimension(:), intent(out) :: stress
    call numodis_addMirrorStressesOnNodes_cc(this%ptr, npts, pts, stress)
  end subroutine numodis_addMirrorStressesOnNodes

end module libnumodis
