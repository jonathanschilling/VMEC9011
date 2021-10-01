subroutine evolve(ierflag)

      use stel_kinds, only: dp
      use name0,      only: cp15,  &
                            cp5,   &
                            c1p0
      use scalars,    only: iter1, &
                            iter2, &
                            irst,  &
                            ndamp, &
                            neqs
      use time,       only: otau,  &
                            otav,  &
                            delt
      use fsqu,       only: fsq,   &
                            fsqr1, &
                            fsqz1, &
                            fsql1
      use xstuff,     only: xc,    &
                            xcdot, &
                            gc

      implicit none

      ! TODO: more elegant definition of these BLAS functions
      real(kind=dp) :: dsum

      integer, intent(out) :: ierflag

      integer       :: i, l
      real(kind=dp) :: fsq1, dtau
      real(kind=dp) :: b1, fac
      integer       :: ndamp1 = 10 ! initial value needed here!

      ! COMPUTE MHD FORCES
      call funct3d

      if (iter2.eq.1 .and. irst.eq.2) then
        ! evaluation of funct3d failed; no need for time step
        ierflag = 1
        return
      endif

      ! COMPUTE DAMPING PARAMETER (DTAU) AND EVOLVE
      ! R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
      if (iter2.eq.iter1) then
        ndamp1 = min(ndamp, 15)
        do i = 1, ndamp1
          otau(i) = cp15/delt
        enddo
      endif

      fsq1 = fsqr1 + fsqz1 + fsql1

      if (iter2.gt.iter1) &
        dtau = min(abs(log(fsq/fsq1)), cp15)

      fsq = fsq1

      if (iter2.le.1) &
        return

      ! shift otau by one index into history ...
      do i = 1, ndamp1-1
        otau(i) = otau(i+1)
      enddo

      ! ... and append new value
      if (iter2.gt.iter1) &
        otau(ndamp1) = dtau/delt

      ! average over last ndamp values in otau
      otav = sum(otau(1:ndamp1))/real(ndamp1)

      dtau = delt*otav

      b1  =       c1p0 - cp5*dtau
      fac = c1p0/(c1p0 + cp5*dtau)

      do l = 1, neqs
        xcdot(l) = fac*(xcdot(l)*b1 + delt*gc(l))
        xc(l)    = xc(l) + xcdot(l)*delt
      enddo

      return

end
