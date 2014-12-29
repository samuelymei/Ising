program Ising
  use precision_m
  use Ising2d_mod
  implicit none
  type ( IsingSys_t ), allocatable :: replica(:)
  integer(kind=4) :: i, j, k
  allocate(replica(1))
  replica(1) = constructor(10,10,1,0.0d0,2.d0)
  call replica(1)%PrintState
  do i = 1, 10000
    call replica(1)%MCmove
    call replica(1)%PrintState
  end do
end program Ising
