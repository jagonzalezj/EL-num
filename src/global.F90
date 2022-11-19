module global

  use DefUtils
  use libnumodis
  use libstressexchange
  
  implicit none
  
!  integer :: numodisStarted = 0 ! flag set to 1 after NUMODIS initialization
  
  type(Fnumodis) :: numodis ! our NUMODIS object wrapper
  type(FstressExchange) :: stressExchange ! our wrapper to exchange stress fields

  integer, dimension(:), allocatable :: NodIdDir
  integer, dimension(:), allocatable :: NodIdNeu
  integer, dimension(:), allocatable :: NodIdTract
  integer, dimension(:), allocatable :: NodIdStr


  REAL(KIND=dp), ALLOCATABLE,target :: Dirichlet(:), Neumann(:), &
       TractionFree(:), Strainrate(:)
  
  real(kind=dp),dimension(:),allocatable,target :: BCdisplacement
  !real(kind=dp),dimension(:),allocatable,target :: BCstress
  double precision,dimension(:),allocatable,target :: BCstress
  !real(kind=dp),dimension(:),allocatable,target :: BCtraction
  double precision,dimension(:),allocatable,target :: BCtraction
  real(kind=dp),dimension(:),allocatable,target :: BCstrain

  REAL(KIND=dp) :: Current_bottom, Current_top, InitBot, InitTop, AppStressLoad

  LOGICAL :: MirrorON = .FALSE. ! mirror / Weygand method activated or not
  integer :: istep = 0 ! current time step
  integer :: fsnapshot ! numodis VTK saving period
  double precision :: time  ! current time
  double precision :: dtime ! theoretical time step
  double precision :: dmax  ! theoretical maximum displacement for a node during one time step

  double precision :: Pstress = 0.D0 !
  double precision :: Pstrain = 0.D0 !

end module global
