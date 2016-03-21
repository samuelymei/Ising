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
