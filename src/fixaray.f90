subroutine fixaray(ierflag)

      use stel_kinds, only: dp
      use name0,      only: czero,   &
                            cp5,     &
                            c1p0,    &
                            c2p0,    &
                            c3p0,    &
                            c8p0
      use name1,      only: nzeta,   &
                            ntheta1, &
                            ntheta2, &
                            mpol1,   &
                            nmax,    &
                            nmax1
      use scalars,    only: twopi,   &
                            dnorm,   &
                            isigng
      use inputdat,   only: nfp, rb, zb
      use mnarray,    only: mscale,  &
                            nscale,  &
                            xmpq,    &
                            xrz3,    &
                            xrz4
      use spectra,    only: faccon
      use trignew,    only: cosmu,   &
                            sinmu,   &
                            cosmui,  &
                            cosmum,  &
                            sinmum,  &
                            cosmumi, &
                            cosnv,   &
                            sinnv,   &
                            cosnvn,  &
                            sinnvn

      implicit none

      ! TODO: more elegant definition of these BLAS functions
      real(kind=dp) :: dsum

      integer, intent(out) :: ierflag

      integer       :: i, j, m, n
      real(kind=dp) :: arg, argi, argj
      real(kind=dp) :: power
      real(kind=dp) :: t1
      real(kind=dp) :: rtest, ztest

      !         INDEX OF LOCAL VARIABLES
      !
      ! mscale  array for norming theta-trig functions (internal use only)
      ! nscale  array for norming zeta -trig functions (internal use only)

      ! COMPUTE TRIGONOMETRIC FUNCTION ARRAYS

      twopi = c8p0*atan(c1p0)
      dnorm = c1p0/real(nzeta*(ntheta2-1))

      do i = 1, ntheta2
        argi = twopi * real(i-1)/real(ntheta1)
        do m = 0, mpol1
          arg = argi * real(m)

          cosmu(i,m) = cos(arg)*mscale(m)
          sinmu(i,m) = sin(arg)*mscale(m)

          cosmui(i,m) = cosmu(i,m)
          if (i.eq.1 .or. i.eq.ntheta2) &
            cosmui(i,m) = cp5*cosmui(i,m)

          cosmum (i,m) =  cosmu (i,m)*real(m)
          sinmum (i,m) = -sinmu (i,m)*real(m)
          cosmumi(i,m) =  cosmui(i,m)*real(m)
        enddo
      enddo

      do j = 1, nzeta
        argj = twopi * real(j-1)/real(nzeta)
        do n = 0, nmax
          arg = argj * real(n)

          cosnv (j,n) = cos(arg)*nscale(n)
          sinnv (j,n) = sin(arg)*nscale(n)
          cosnvn(j,n) = cosnv(j,n)*real(n*nfp)
          sinnvn(j,n) =-sinnv(j,n)*real(n*nfp)
        enddo
      enddo

      do m = 0, mpol1
        power = -cp5*real(m) ! -m/2

        xmpq(m,1) = real(m*(m-1))
        xmpq(m,2) = real(m**4)
        xmpq(m,3) = real(m**5)

        ! later: 1/xpmq(m,1)**2 = 1/[m(m-1)]^2 = 1/{(m^2-m)^2} = 1/{m^4 - 2*m^3 + m^2}
        t1 = -cp5*real(isigng)*dnorm/real((1+m)**4) ! TODO: what is (1+m)^4 ???
        do n = 0, nmax

          faccon(n,m) = t1/(mscale(m)*nscale(n))**2
          if (m.eq.0 .or. n.eq.0) &
            faccon(n,m) = cp5*faccon(n,m)

          xrz3(n,m) = c2p0**(power + c1p0)
          xrz4(n,m) =-c3p0**power
        enddo
      enddo

      ! CHECK SIGN OF JACOBIAN (SHOULD BE SAME AS ISIGNG)
      rtest = sum(rb(0:nmax,1,1))
      ztest = sum(zb(0:nmax,1,2))
      print *, "rtest = ", rtest
      print *, "ztest = ", ztest
      if ((rtest*ztest*real(isigng)) .ge. czero) &
        ierflag = 5

      return

end
