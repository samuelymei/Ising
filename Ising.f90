program Ising
  use precision_m
  use Ising2d_mod
  use mpi
  implicit none
  integer ( kind = 4 ) :: ierr
  integer ( kind = 4 ) :: numprocs ! number of MPI processes
  integer ( kind = 4 ) :: mytid ! id of current MPI process
  integer ( kind = 4 ) :: mpiStatus( MPI_Status_size )
  logical :: master

  type ( IsingSys_t ) :: replica
  real( kind = fp_kind ) :: bExt
  real( kind = fp_kind ), allocatable :: temperatures( : )
  real( kind = fp_kind ), allocatable :: energies( : )
  integer( kind = 4 ), allocatable :: nTotalSpins( : )
  real( kind = fp_kind ), parameter :: T_Max = 3, T_Min = 1
  real( kind = fp_kind ) :: currTemperature
  integer( kind = 4 ), allocatable :: temperature_index( : )
  real( kind = fp_kind ), allocatable :: temperature_target( : )
  real( kind = fp_kind ) :: newTemperature

  integer( kind = 4 ) :: numEx = 1000000
  integer( kind = 4 ) :: numRelax = 50
  integer( kind = 4 ) :: indexP, indexEx, indexMove
  logical :: IsOdd
  logical, allocatable :: willBeExchanged(:)

  real( kind = fp_kind ), allocatable :: exchangeRate(:)

  real( kind = fp_kind ) :: MyUniformRand
  
  character( len = : ), allocatable :: outstring

! initial MPI environment
  call MPI_Init ( ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, mytid, ierr )

  master = mytid == 0

! initialize temperatures
  if( master ) then
    allocate( temperatures( numprocs ) )
    allocate( energies( numprocs ) )
    allocate( nTotalSpins( numprocs ) )
    allocate( temperature_index( numprocs ) )
    allocate( exchangeRate( numprocs ) )
    allocate( willBeExchanged( numprocs ) )
    allocate( temperature_target( numprocs ) )

    call AssignTemperature(numprocs, temperatures, T_Max, T_Min)
    do indexP = 1, numprocs - 1 
      call MPI_Send( temperatures( indexP + 1 ), 1, MPI_REAL8, indexP, 0, MPI_COMM_WORLD, ierr)
    end do
    currTemperature = temperatures( 1 )
    forall( indexP = 1 : numprocs )
      temperature_index( indexP ) = indexP
    end forall
  else
    call MPI_Recv( currTemperature, 1, MPI_REAL8, 0, 0, MPI_COMM_WORLD, mpiStatus, ierr)
  end if

  
!  print*,' Process ID ', mytid, ' Temperature:', currTemperature

  bExt = 0.5d0
  replica = constructor( 10, 10, 2, bExt, currTemperature )

! print initial energies
  call MPI_Gather( replica%Epot, 1, MPI_REAL8, energies, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
  call MPI_Gather( replica%nTotalSpin, 1, MPI_INTEGER, nTotalSpins, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
  IndexEx = 0 
  if( master ) then
    call printStates( numprocs, IndexEx, energies )
  end if

  do indexEx = 1, numEx
! Monte Carlo evolution
    do indexMove = 1, numRelax
      call replica%MCmove
    end do

! gather energies for exchange
    call MPI_Gather( replica%Epot, 1, MPI_REAL8, energies, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
    call MPI_Gather( replica%nTotalSpin, 1, MPI_INTEGER, nTotalSpins, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
!    if(indexEx == 10000)then
!      outstring = 'Replica 9999'
!      write(outstring,*)mytid
!      outstring = 'Replica '//trim(adjustl(outstring))
!      call replica%PrintState(outstring)
!      if(master) write(*,'(8I10)') temperature_index
!      if(master) write(*,'(8F10.3)') energies( temperature_index )
!    end if
    if(master)write(88,'(I10,8F10.3)') indexEx, energies(temperature_index)    
    if(master)write(99,'(I10,8I10)') indexEx, nTotalSpins(temperature_index)    

    if( master ) then
      if ( IsOdd( indexEx ) ) then
! check 1:2, 3:4, n-1:n
        willBeExchanged( : ) = .false.
        do indexP = 1, numprocs - 1, 2
          exchangeRate( temperature_index( indexP ) ) = &
                & ( 1.d0 / temperatures( temperature_index( indexP ) ) - & 
                    1.d0 / temperatures( temperature_index( indexP + 1 ) ) ) &
               & *( energies( temperature_index( indexP + 1 ) ) - &
                  & energies( temperature_index( indexP ) ) )
          exchangeRate( temperature_index( indexP ) ) = &
             & exp( - exchangeRate( temperature_index( indexP ) ) )
          if( exchangeRate( temperature_index( indexP ) )  > 1.d0 ) &
              exchangeRate( temperature_index( indexP ) ) = 1.d0
          exchangeRate( temperature_index( indexP + 1 ) ) = &
            & exchangeRate( temperature_index( indexP ) )
          if ( exchangeRate( temperature_index( indexP ) ) > MyUniformRand() ) then
            willBeExchanged( temperature_index( indexP ) ) = .true.   
            willBeExchanged( temperature_index( indexP + 1 ) ) = .true.   
            temperature_target ( temperature_index( indexP ) ) = temperatures( indexP + 1 )       
            temperature_target ( temperature_index( indexP + 1 ) ) = temperatures( indexP )
            call swap_int( temperature_index( indexP ), temperature_index( indexP + 1 ))
          else
            willBeExchanged( temperature_index( indexP ) ) = .false.   
            willBeExchanged( temperature_index( indexP + 1 ) ) = .false.   
            temperature_target ( temperature_index( indexP ) ) = temperatures( indexP )
            temperature_target ( temperature_index( indexP + 1 ) ) = temperatures( indexP + 1 )
          end if
        end do
      else
! check 2:3, 3:4, n-2:n-1
        willBeExchanged( : ) = .false.
        do indexP = 2, numprocs - 2, 2
          exchangeRate( temperature_index( indexP ) ) = &
                & ( 1.d0 / temperatures( temperature_index( indexP ) ) - &  
                    1.d0 / temperatures( temperature_index( indexP + 1 ) ) ) & 
               & *( energies( temperature_index( indexP + 1 ) ) - &
                  & energies( temperature_index( indexP ) ) )
          exchangeRate( temperature_index( indexP ) ) = &
             & exp( - exchangeRate( temperature_index( indexP ) ) )
          if( exchangeRate( temperature_index( indexP ) )  > 1.d0 ) &
              exchangeRate( temperature_index( indexP ) ) = 1.d0
          exchangeRate( temperature_index( indexP + 1 ) ) = &
            & exchangeRate( temperature_index( indexP ) )
          if ( exchangeRate( temperature_index( indexP ) ) > MyUniformRand() ) then
            willBeExchanged( temperature_index( indexP ) ) = .true.
            willBeExchanged( temperature_index( indexP + 1 ) ) = .true.
            temperature_target( temperature_index( indexP ) ) = temperatures( indexP + 1 )
            temperature_target( temperature_index( indexP + 1 ) ) = temperatures( indexP )
            call swap_int( temperature_index( indexP ), temperature_index( indexP + 1 ))
          else
            willBeExchanged( temperature_index( indexP ) ) = .false.   
            willBeExchanged( temperature_index( indexP + 1 ) ) = .false.   
            temperature_target ( temperature_index( indexP ) ) = temperatures( indexP )
            temperature_target ( temperature_index( indexP + 1 ) ) = temperatures( indexP + 1 )
          end if
        end do
        exchangeRate( temperature_index( 1 ) ) = 0.d0
        temperature_target( temperature_index( 1 ) ) = temperatures( 1 )
        exchangeRate( temperature_index( numprocs ) ) = 0.d0
        temperature_target( temperature_index( numprocs ) ) = temperatures( numprocs )
      end if
    end if

    if( master ) then
      do indexP = 1, numprocs - 1 
        call MPI_Send( temperature_target( temperature_index( indexP + 1 ) ), 1, MPI_REAL8, indexP, 0, MPI_COMM_WORLD, ierr)
      end do
      newTemperature = temperature_target( temperature_index( 1 ) )
      call replica%ChangeTemperature ( newTemperature )
!      write(77,*) willBeExchanged(temperature_index)
!      write(77,'(6F10.3)') temperature_target
    else 
      call MPI_Recv( newTemperature, 1, MPI_REAL8, 0, 0, MPI_COMM_WORLD, mpiStatus, ierr)
      call replica%ChangeTemperature ( newTemperature )
    end if
    
  end do

  call replica%Finalize()

  if( master ) then
    deallocate( temperatures )
    deallocate( energies )
    deallocate( nTotalSpins )
    deallocate( temperature_index )
    deallocate( exchangeRate )
    deallocate( willBeExchanged )
    deallocate( temperature_target )
  end if
  call MPI_Finalize ( ierr )
end program Ising

subroutine printStates( n, indexEx, energies )
  use precision_m
  implicit none
  integer( kind = 4 ), intent(in) :: n
  integer( kind = 4 ), intent(in) :: indexEx
  real( kind = fp_kind ) :: energies( n )
  write(*, '(A,I5)')'T = ', IndexEx
  write(*, '(8F10.3)') energies
end subroutine printStates

subroutine AssignTemperature( n, temperatures, t_max, t_min )
  use precision_m
  implicit none
  integer( kind = 4 ), intent( in ) :: n
  real( kind = fp_kind ), intent( out ) :: temperatures( n )
  real( kind = fp_kind ), intent( in ) :: t_max, t_min
  real( kind = fp_kind ) :: c, alpha
  integer( kind = 4 ) :: i, j, k
  alpha = ( n - 1 ) / ( log( t_max ) - log( t_min ) )
  c = t_min / exp( 1.d0 / alpha )
  do i = 1, n
    temperatures( i ) = c*exp( i / alpha )
  end do
end subroutine AssignTemperature
