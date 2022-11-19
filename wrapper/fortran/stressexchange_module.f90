module libstressexchange

  use iso_c_binding


  private

  public :: FstressExchange

  include "stressexchange_cc2fortran.bind"


  type FstressExchange

     private

     type(c_ptr) :: ptr

   contains

     final :: destroy_stressExchange

     procedure :: getNumberOfNodes => stressExchange_getNumberOfNodes

     procedure :: getCoordinates => stressExchange_getCoordinates

     procedure :: writeVTK => stressExchange_writeVTK
     
     procedure :: flushFields => stressExchange_flush

     procedure :: pushField => stressExchange_pushField

     procedure :: addFields => stressExchange_addFields
     
  end type FstressExchange


  interface FstressExchange
     procedure construct_stressExchange
  end interface FstressExchange



contains


  !====================================================
  ! construct_stressExchange
  !----------------------------------------------------
  !> "Constructor" of the StressExchange container
  !----------------------------------------------------
  !! @param xmin bottom-left corner of the grid
  !! @param xmax upper-right corner of the grid
  !! @param npts number of nodes along the 3 axes
  !!===================================================
  function construct_stressExchange(xmin,xmax,npts)
    implicit none
    type(FstressExchange) :: construct_stressExchange
    real(kind=c_double), intent(in) :: xmin(*)
    real(kind=c_double), intent(in) :: xmax(*)
    integer(c_int), intent(in)      :: npts(*)    
    construct_stressExchange%ptr = construct_stressExchange_cc(xmin,xmax,npts)
  end function construct_stressExchange

  !====================================================
  ! subroutine destroy_stressExchange
  !----------------------------------------------------
  !> "Destructor" of the StressExchange container
  !----------------------------------------------------
  !! @param this our StressExchange container
  !!===================================================
  subroutine destroy_stressExchange(this)
    implicit none
    type(FstressExchange) :: this
    call destroy_stressExchange_cc(this%ptr)
  end subroutine destroy_stressExchange

  !==========================================================
  ! function stressExchange_getNumberOfNodes
  !----------------------------------------------------------
  !> get the number of nodes of the StressExchange container
  !----------------------------------------------------------
  !! @param this our StressExchange container
  !! @return number of nodes
  !!=========================================================
  integer function stressExchange_getNumberOfNodes(this)
    implicit none
    class(FstressExchange) :: this
    stressExchange_getNumberOfNodes = stressExchange_getNumberOfNodes_cc(this%ptr)
  end function stressExchange_getNumberOfNodes
  
  !==========================================================
  ! function stressExchange_getCoordinates
  !----------------------------------------------------------
  !> return a pointer to the array containing the mesh nodes
  !----------------------------------------------------------
  !! @param this our StressExchange container
  !! @return pointer to the corresponding array
  !==========================================================
  function stressExchange_getCoordinates(this)
    use iso_c_binding
    implicit none
    class(FstressExchange) :: this
    real(kind=c_double), pointer :: stressExchange_getCoordinates(:)

    integer(c_int) :: nNodes
    type(c_ptr) :: cCoordinates
    real(kind=c_double), pointer :: fCoordinates(:)

    ! get number of nodes of the mesh
    nNodes = stressExchange_getNumberOfNodes_cc(this%ptr)

    ! get C pointer to the coordinates of the mesh
    call stressExchange_getCoordinates_cc(this%ptr,cCoordinates)

    ! convert this C pointer to a proper fortran pointer to an array of positions
    call c_f_pointer(cCoordinates,fCoordinates,(/3*nNodes/))

    ! return this pointer
    stressExchange_getCoordinates => fCoordinates
    
  end function stressExchange_getCoordinates

  !================================================
  ! subroutine stressExchange_writeVTK
  !------------------------------------------------
  !>  Save all fields to a VTK file
  !------------------------------------------------
  !! @param this our StressExchange class
  !! @param filename VTK filename
  !================================================
  subroutine stressExchange_writeVTK(this,filename)
    use iso_c_binding
    implicit none
    class(FstressExchange) :: this
    character(kind=c_char,len=*), intent(in) :: filename
    call stressExchange_writeVTK_cc(this%ptr,filename//C_NULL_CHAR) ! add C_NULL_CHAR means "end of character array"
  end subroutine stressExchange_writeVTK

  !================================================
  ! subroutine stressExchange_flush
  !------------------------------------------------
  !> Flush all fields contained in stressExchange
  !------------------------------------------------
  !! @param this our StressExchange class
  !================================================  
  subroutine stressExchange_flush(this)
    implicit none
    class(FstressExchange) :: this
    call stressExchange_flush_cc(this%ptr)
  end subroutine stressExchange_flush

  !=========================================================
  ! subroutine stressExchange_pushField
  !---------------------------------------------------------
  !> Store a stress field in stressExchange
  !---------------------------------------------------------
  !! @param this our StressExchange class
  !! @fieldname name of the field in the class
  !! @stress stress field
  !=========================================================
  subroutine stressExchange_pushField(this,fieldname,stress)
    implicit none
    class(FstressExchange) :: this
    character(kind=c_char,len=*), intent(in) :: fieldname
    double precision, dimension(:), intent(out) :: stress
    call stressExchange_pushField_cc(this%ptr,fieldname//C_NULL_CHAR,stress)
  end subroutine stressExchange_pushField

  !========================================================
  ! subroutine stressExchange_addFields
  !--------------------------------------------------------
  !> Compute and store the total stress field
  !--------------------------------------------------------
  !! @param this our StressExchange class
  !! @fieldname name of the field in the class
  !========================================================
  subroutine stressExchange_addFields(this,fieldname)
    implicit none
    class(FstressExchange) :: this
    character(kind=c_char,len=*), intent(in) :: fieldname
    call stressExchange_addFields_cc(this%ptr,fieldname//C_NULL_CHAR)
  end subroutine stressExchange_addFields
  
end module libstressexchange
