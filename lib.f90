function IsInRange(i, imin, imax)
  implicit none
  logical :: isInRange
  integer(kind=4), intent(in) :: i, imin, imax
  isInRange = .false.
  if(i<=imax .and. i>=imin) isInRange = .true.
end function IsInRange

function Rmsd(n, a1, a2)
  use precision_m
  implicit none
  integer(kind=4), intent(in) :: n
  real(kind=fp_kind), intent(in) :: a1(n), a2(n)
  real(kind=fp_kind) :: rmsd
  integer(kind=4) :: i
  if(n==0)then
    rmsd = 0.d0
    return
  end if
  rmsd = 0.d0
  rmsd = sqrt(sum((a1-a2)**2)/n)
end function Rmsd

subroutine Mean_and_StandardDev(n,x,mean,stdDev)
  use precision_m
  implicit none
  integer(kind=4), intent(in) :: n
  real(kind=fp_kind), intent(in) :: x(n)
  real(kind=fp_kind), intent(out) :: mean, stdDev
  mean = sum(x)/n
  stdDev = sqrt(sum((x-mean)**2)/n)
end subroutine Mean_and_StandardDev

subroutine Mean_and_StandardErr(n,x,mean,stdErr)
  use precision_m
  implicit none
  integer(kind=4), intent(in) :: n
  real(kind=fp_kind), intent(in) :: x(n)
  real(kind=fp_kind), intent(out) :: mean, stdErr
  mean = sum(x)/n
  stdErr = sqrt(sum((x-mean)**2))/n
end subroutine Mean_and_StandardErr

function MyGaussianRand() result (y)
! polar form of the Box-Muller transformation
  use precision_m
  implicit none
  real(kind=fp_kind) :: x(2)
  real(kind=fp_kind) :: w
  real(kind=fp_kind) :: y(2)
  real(kind=fp_kind) :: MyUniformRand
  do while (.true.)
    x(1) = 2.0 * MyUniformRand() - 1.0
    x(2) = 2.0 * MyUniformRand() - 1.0
    w = x(1) * x(1) + x(2) * x(2)
    if(w < 1.0d0) exit
  end do
  w = sqrt( ( -2.0 * log(w) ) / w )
  y(1) = x(1) * w
  y(2) = x(2) * w
end function MyGaussianRand

function MyUniformRand()
  use precision_m
  implicit none
  integer(kind=4), save :: initialized = 0
  real(kind=fp_kind) :: MyUniformRand
  real(kind=fp_kind) :: r
  if( initialized .eq. 0 )then
    CALL Init_Random_Seed()         ! see example of RANDOM_SEED
    initialized = 1
  end if
  CALL RANDOM_NUMBER(r)
  myUniformRand = r
end function MyUniformRand

subroutine Init_Random_Seed()
  use iso_fortran_env, only: int64
  USE IFPORT, only: getpid
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t

  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
     read(un) seed
     close(un)
  else
     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
     call system_clock(t)
     if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24_int64 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
     end if
     pid = getpid()
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
        seed(i) = lcg(t)
     end do
  end if
  call random_seed(put=seed)
contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
end subroutine Init_Random_Seed

function IsOdd( i )
  implicit none
  integer( kind = 4 ) :: i
  logical :: IsOdd
  IsOdd = .true.
  if( i / 2 * 2 == i ) IsOdd = .false.
end function IsOdd

function IsEven( i ) 
  implicit none
  integer( kind = 4 ) :: i
  logical :: IsEven
  IsEven = .false.
  if( i / 2 * 2 == i ) IsEven = .true.
end function IsEven

!interface swap
!  subroutine swap_int( i1, i2 )
!    implicit none
!    integer( kind = 4 ), intent( in out ) :: i1, i2
!  end subroutine swap_int 
!
!  subroutine swap_real( f1, f2 )
!    implicit none
!    real( kind = 4 ), intent( in out ) :: f1, f2
!  end subroutine swap_real
!
!  subroutine swap_double( d1, d2 )
!    implicit none
!    real( kind = 8 ), intent( in out ) :: d1, d2
!  end subroutine swap_double
!
!end interface swap

subroutine swap_int( i1, i2 )
  implicit none
  integer( kind = 4 ), intent( in out ) :: i1, i2
  integer( kind = 4 ) :: itmp
  itmp = i1
  i1 = i2
  i2 = itmp
end subroutine swap_int

subroutine swap_real( f1, f2 )
  implicit none
  real( kind = 4 ), intent( in out ) :: f1, f2
  real( kind = 4 ) :: ftmp
  ftmp = f1 
  f1 = f2
  f2 = ftmp
end subroutine swap_real

subroutine swap_double( d1, d2 )
  implicit none
  real( kind = 8 ), intent( in out ) :: d1, d2
  real( kind = 8 ) :: dtmp
  dtmp = d1 
  d1 = d2
  d2 = dtmp
end subroutine swap_double
