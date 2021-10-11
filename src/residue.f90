subroutine residue(gcr, gcz, gcl, bz0, work)

      use stel_kinds, only: dp
      use name0, only: czero, c1p0, c2p0
      use name1, only: mpol1, nmax, mnd2
      use scalars, only: ns, mns, iter1, hs, meven, modd
      use fsqu, only: fsqr, fsqz, fsql, fsqr1, fsqz1, fsql1, &
                      fnorm, fnorm1, fedge
      use precond, only: arm, brm, ard, brd, cr, &
                          azm, bzm, azd, bzd
      use xstuff, only: gc
      use scalefac, only: faclam

      implicit none

      ! TODO: more elegant definition of these BLAS functions
      real(kind=dp) :: ddot
      external dscal

      real(kind=dp) :: gcr(ns,0:nmax,0:mpol1,2)
      real(kind=dp) :: gcz(ns,0:nmax,0:mpol1,2)
      real(kind=dp) :: gcl(ns,0:nmax,0:mpol1,2)

      real(kind=dp) :: work(mns,6) ! TODO: why 6 --> only 5 used.. ???
      real(kind=dp) :: bz0

      integer       :: js, n, l
      real(kind=dp) :: fac

      ! IMPOSE M=1 MODE CONSTRAINT (ON Z CS, NTYPE = 1)
      if (fsqz.lt.1.e-8 .and. iter1.ne.1) then
        do n = 1, nmax
          do js = 2, ns
            gcz(js,n,1,1) = czero
          end do
        end do
      end if

      ! CONSTRUCT INVARIANT RESIDUALS
      bz0  = hs * c2p0/bz0**2
      fsql = bz0 * ddot(2*mns, gcl, 1, gcl, 1)
      call getfsq(gcr, gcz, fsqr, fsqz, fnorm, meven)

      fedge = fnorm * (  ddot(mnd2, gcr(ns,0,0,1), ns, gcr(ns,0,0,1), ns) &
                       + ddot(mnd2, gcz(ns,0,0,1), ns, gcz(ns,0,0,1), ns) )

      ! PERFORM PRECONDITIONING AND COMPUTE RESIDUES
      call scalfor(gcr, arm, brm, ard, brd, cr, work, work(1,2), work(1,3), work(1,4), work(1,5))
      call scalfor(gcz, azm, bzm, azd, bzd, cr, work, work(1,2), work(1,3), work(1,4), work(1,5))

      ! REDUCE R-Z FORCES IF INITIALLY VERY LARGE
      ! if in first iteration: scale  after computing fsqX1
      !      other iterations: scale before computing fsqX1
      fac = c1p0/(c1p0 + (fsqr+fsqz))
      if (iter1.gt.1) then
        ! any further iteration (but not first)
        call dscal(4*mns, fac, gc, 1)
      end if

      call getfsq(gcr, gcz, fsqr1, fsqz1, fnorm1, modd)

      if (iter1.eq.1) then
        ! only in first iteration
        call dscal(4*mns, fac, gc, 1)
      end if

      ! CONSTRUCT PRECONDITIONED (SCALED) LAMBDA FORCES
      do l = 1, 2*mns
        gc(l+4*mns) = faclam(l)*gc(l+4*mns)
      end do
      fsql1 = hs * ddot(2*mns, gcl, 1, gcl, 1)

      ! output:
      ! 1.  "invariant" force residuals: fsqr,  fsqz  (scaled by fnorm) and fsql  (scaled by hs/bz0^2)
      ! 2. precondioned force residuals: fsqr1, fsqz1 (scaled by fnorm) and fsql1 (scaled by hs)

      return
end
