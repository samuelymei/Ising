!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A 2D Ising model                                                       !
! Written by                                                             !
!                                   Ye Mei                               !
!                        East China Normal University                    !
!                                 12/30/2014                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module Ising2d_mod
  use precision_m
  implicit none
  private
  public :: IsingSys_t, Constructor
  type particle_t
    integer( kind = 4 ) :: iSpin
    integer( kind = 4 ) :: nNeighbor
    integer( kind = 4 ), pointer :: neighbors_p(:) => null()
    real( kind = fp_kind ) :: energy
  end type particle_t

  type IsingSys_t
    integer( kind =4 ) :: nx, ny
    integer( kind =4 ) :: nParticles
    type(particle_t), allocatable :: particles(:)
    real( kind = fp_kind ) :: ePot  ! \sigma(i,j>i) -0.5*ispin(i)*ispin(j) + \sigma(i) B*ispin(i)
    real( kind = fp_kind ) :: bExt
    real( kind = fp_kind ) :: temperature
    integer( kind = 4 )  :: nTotalSpin

  contains
    procedure :: Initialize
    procedure :: FindNeighbor
    procedure :: CalcIndivEner
    procedure :: GetEnergyQuick
    procedure, public :: GetEnergy
    procedure :: GetTotalSpin
    procedure, public :: MCmove
    procedure, public :: Finalize
    procedure, public :: PrintState
    procedure, public :: ChangeTemperature
  end type IsingSys_t

  contains
    function Constructor( nrow, ncol, imean, Bext, temper )
      implicit none
      type ( IsingSys_t ) :: Constructor
      integer( kind = 4 ) , intent(in) :: nrow, ncol
      integer( kind = 4 ) , intent(in) :: imean
      real( kind = fp_kind ), intent(in) :: Bext, temper
      call Constructor%initialize( nrow, ncol, imean, Bext, temper )
    end function Constructor

    subroutine Initialize( this, nrow, ncol, imean, Bext, temper )
      implicit none
      class( IsingSys_t ) :: this
      integer( kind = 4 ) , intent(in) :: nrow, ncol
      integer( kind = 4 ) , intent(in) :: imean
      real( kind = fp_kind ), intent(in) :: Bext, temper
      integer( kind = 4 )  :: indexP
      integer( kind = 4 )  :: indexX, indexY
      real( kind = fp_kind ) :: myUniformRand
      real( kind = fp_kind ) :: rannum
      this%nx = nrow
      this%ny = ncol
      this%nParticles = this%nx * this%ny
      this%bExt = Bext
      this%temperature = temper
      allocate(this%particles(this%nParticles))
      if ( imean == 1 ) then
        this%particles(:)%iSpin = 1
      else
        do indexY = 1, this%ny
          do indexX = 1, this%nx
            indexP = (indexY - 1) * this%nx + indexX
            rannum = myUniformRand()
            if ( rannum >= 0.5d0 ) then
              this%particles(indexP)%iSpin = -1
            else
              this%particles(indexP)%iSpin = 1
            end if
          end do
        end do
      end if
      call FindNeighbor( this )
      do indexP = 1, this%nParticles
        call CalcIndivEner( this, indexP )
      end do
      call GetEnergy( this, this%ePot )
      call GetTotalSpin( this, this%nTotalSpin )
    end subroutine Initialize
 
    subroutine FindNeighbor(this)
      implicit none
      class( IsingSys_t ) :: this
      integer( kind = 4 )  :: indexP
      integer( kind = 4 )  :: indexX, indexY
 
      do indexY = 1, this%ny
        do indexX = 1, this%nx
          indexP = (indexY - 1) * this%nx + indexX
          if ( indexX == 1 .and. indexY == 1 ) then
            this%particles(indexP)%nNeighbor = 2
            allocate(this%particles(indexP)%neighbors_p(this%particles(indexP)%nNeighbor))
            this%particles(indexP)%neighbors_p(1) = indexP + 1
            this%particles(indexP)%neighbors_p(2) = indexP + this%nx
          else if ( indexX == this%nx .and. indexY == this%ny ) then
            this%particles(indexP)%nNeighbor = 2
            allocate(this%particles(indexP)%neighbors_p(this%particles(indexP)%nNeighbor))
            this%particles(indexP)%neighbors_p(1) = indexP - 1
            this%particles(indexP)%neighbors_p(2) = indexP - this%nx
          else if ( indexX == 1 .and. indexY == this%ny ) then
            this%particles(indexP)%nNeighbor = 2
            allocate(this%particles(indexP)%neighbors_p(this%particles(indexP)%nNeighbor))
            this%particles(indexP)%neighbors_p(1) = indexP + 1
            this%particles(indexP)%neighbors_p(2) = indexP - this%nx
          else if ( indexX == this%nx .and. indexY == 1 ) then
            this%particles(indexP)%nNeighbor = 2
            allocate(this%particles(indexP)%neighbors_p(this%particles(indexP)%nNeighbor))
            this%particles(indexP)%neighbors_p(1) = indexP - 1
            this%particles(indexP)%neighbors_p(2) = indexP + this%nx
          else if ( indexX == 1 .and. indexY > 1 .and. indexY < this%ny ) then
            this%particles(indexP)%nNeighbor = 3
            allocate(this%particles(indexP)%neighbors_p(this%particles(indexP)%nNeighbor))
            this%particles(indexP)%neighbors_p(1) = indexP - this%nx
            this%particles(indexP)%neighbors_p(2) = indexP + 1
            this%particles(indexP)%neighbors_p(3) = indexP + this%nx
          else if ( indexX == this%nx .and. indexY > 1 .and. indexY < this%ny ) then
            this%particles(indexP)%nNeighbor = 3
            allocate(this%particles(indexP)%neighbors_p(this%particles(indexP)%nNeighbor))
            this%particles(indexP)%neighbors_p(1) = indexP - this%nx
            this%particles(indexP)%neighbors_p(2) = indexP - 1
            this%particles(indexP)%neighbors_p(3) = indexP + this%nx
          else if ( indexX > 1 .and. indexX < this%nx .and. indexY == 1 ) then
            this%particles(indexP)%nNeighbor = 3
            allocate(this%particles(indexP)%neighbors_p(this%particles(indexP)%nNeighbor))
            this%particles(indexP)%neighbors_p(1) = indexP - 1
            this%particles(indexP)%neighbors_p(2) = indexP + this%nx
            this%particles(indexP)%neighbors_p(3) = indexP + 1
          else if ( indexX > 1 .and. indexX < this%nx .and. indexY == this%ny ) then
            this%particles(indexP)%nNeighbor = 3
            allocate(this%particles(indexP)%neighbors_p(this%particles(indexP)%nNeighbor))
            this%particles(indexP)%neighbors_p(1) = indexP - 1
            this%particles(indexP)%neighbors_p(2) = indexP - this%nx
            this%particles(indexP)%neighbors_p(3) = indexP + 1
          else
            this%particles(indexP)%nNeighbor = 4
            allocate(this%particles(indexP)%neighbors_p(this%particles(indexP)%nNeighbor))
            this%particles(indexP)%neighbors_p(1) = indexP - 1
            this%particles(indexP)%neighbors_p(2) = indexP - this%nx
            this%particles(indexP)%neighbors_p(3) = indexP + 1
            this%particles(indexP)%neighbors_p(4) = indexP + this%nx
          end if
        end do
      end do
    end subroutine FindNeighbor
 
    subroutine CalcIndivEner( this, indexP )
      implicit none
      class( IsingSys_t ) :: this
      integer( kind = 4 ) , intent(in) :: indexP
      real( kind = fp_kind ) :: indivEner
      integer( kind = 4 )  :: indexN
      this%particles(indexP)%energy = this%bExt * this%particles(indexP)%iSpin
      do indexN = 1, this%particles(indexP)%nNeighbor
        this%particles(indexP)%energy = this%particles(indexP)%energy - 0.5* this%particles(indexP)%iSpin * &
                   & this%particles(this%particles(indexP)%neighbors_p(indexN))%iSpin
      end do
    end subroutine CalcIndivEner
 
    subroutine GetEnergyQuick( this, energy )
      implicit none
      class( IsingSys_t ) :: this
      real( kind = fp_kind ), intent(out) :: energy
      integer( kind = 4 )  :: indexP
      energy = 0.d0
      do indexP = 1, this%nParticles
        energy = energy + this%particles(indexP)%energy
        energy = energy - 0.5 * this%bExt * this%particles(indexP)%iSpin
      end do
      energy = energy/2.d0
    end subroutine GetEnergyQuick
 
    subroutine GetEnergy( this, energy )
      implicit none
      class( IsingSys_t ) :: this
      real( kind = fp_kind ), intent(out) :: energy
      integer( kind = 4 )  :: indexX, indexY
      integer( kind = 4 )  :: indexP, indexN
      energy = 0.d0
      do indexX = 1, this%nx
        do indexY = 1, this%ny
          indexP = (indexY - 1) * this%nx + indexX
          do indexN = 1, this%particles(indexP)%nNeighbor
            if( this%particles(indexP)%neighbors_p(indexN) <= indexP ) cycle
            energy = energy - 0.5* this%particles(indexP)%iSpin &
                                & * this%particles(this%particles(indexP)%neighbors_p(indexN))%iSpin
          end do
          energy = energy + this%bExt * this%particles(indexP)%iSpin
        end do
      end do
    end subroutine GetEnergy
 
    subroutine GetTotalSpin( this, nSpin )
      implicit none
      class( IsingSys_t ) :: this
      integer( kind = 4 ) , intent(out) :: nSpin
      nSpin = sum(this%particles(:)%iSpin)
    end subroutine GetTotalSpin

    subroutine MCmove( this )
      implicit none
      class ( IsingSys_t ) :: this
      integer( kind = 4 )  :: indexP ! index for particles
      integer( kind = 4 )  :: indexN ! index for neighbors
      real( kind = fp_kind ) :: deltaE
      real( kind = fp_kind ) :: myUniformRand
      logical :: iAccept
      real( kind = fp_kind ) :: rannum
      indexP = myUniformRand() * this%nParticles + 1
      deltaE = -2 * this%particles(indexP)%energy
      if( deltaE < 0.d0 ) then
        iAccept = .True.
      else
        rannum = myUniformRand()
        if ( exp(-deltaE/this%temperature) > rannum ) then
          iAccept = .True.
        else
          iAccept = .False.
        end if
      end if
      if ( iAccept ) then
        this%nTotalSpin = this%nTotalSpin - (this%particles(indexP)%iSpin)*2
        this%particles(indexP)%iSpin = - this%particles(indexP)%iSpin
        this%particles(indexP)%energy = - this%particles(indexP)%energy
        do indexN = 1, this%particles(indexP)%nNeighbor
          call calcIndivEner(this,this%particles(indexP)%neighbors_p(indexN))
        end do
        this%ePot = this%ePot + deltaE
      end if
    end subroutine MCmove
 
    subroutine Finalize( this )
      implicit none
      class( IsingSys_t ) :: this
      integer( kind = 4 )  :: indexP
      deallocate(this%particles)
    end subroutine Finalize

    subroutine PrintState( this, comments )
      implicit none
      class( IsingSys_t ) :: this
      character( len = * ) :: comments
      write( *, '(A,3F10.3,I10)' ) trim(adjustl(comments)), this%ePot, this%bExt, this%temperature, this%nTotalSpin
    end subroutine PrintState

    subroutine ChangeTemperature( this, tNew )
      implicit none
      class( IsingSys_t ) :: this
      real( kind = fp_kind ) :: tNew
      this%temperature = tNew
    end subroutine ChangeTemperature

end module Ising2d_mod
