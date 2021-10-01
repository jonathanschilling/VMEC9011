module name1

      use stel_kinds, only: dp

      implicit none

      integer, parameter :: nsd     = 10 ! 41       ! number of flux surfaces
      integer, parameter :: mpol    =  7 ! 9        ! number of poloidal modes incl. 0
      integer, parameter :: nmax    =  6 ! 13       ! number of toroidal modes 2*nmax+1 (?)
      integer, parameter :: ntheta  = 32 ! 36       ! number of poloidal grid points
      integer, parameter :: nzeta   = 22 ! 44       ! number of toroidal grid points
      integer, parameter :: nvac    =  0            ! fixed-boundary (0) or free-boundary (1)

      integer, parameter :: nsd1    = nsd+1         ! number of flux surfaces + 1
      integer, parameter :: ntheta1 = 2*(ntheta/2)  ! even version of ntheta
      integer, parameter :: ntheta2 = 1+ntheta1/2   ! stellarator-symmetric reduction of ntheta1
      integer, parameter :: nznt    = nzeta*ntheta2 ! number of surface grid points
      integer, parameter :: nrztd   = nznt*nsd+1    ! number of grid points in volume (+1 ?)
      integer, parameter :: mpol1   = mpol-1        ! highest occuring poloidal mode number
      integer, parameter :: nmax1   = nmax+1        ! number of toroidal coefficients per set of basis functions
      integer, parameter :: mnd     = mpol*nmax1    ! total number of Fourier coefficients per set of basis functions
      integer, parameter :: mnd1    = mnd-1         ! total number of Fourier coefficients per set of basis functions - 1
      integer, parameter :: mnd2    = 2*mnd         ! 2 x total number of Fourier coefficients per set of basis functions
      integer, parameter :: mnmax   = nmax1+mpol1*(1+2*nmax) ! total number of Fourier coefficients per surface
      integer, parameter :: neq     = 6*nsd*mnd     ! total number of free parameters

end module name1
