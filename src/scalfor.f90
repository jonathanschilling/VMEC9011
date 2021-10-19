subroutine scalfor(gcx, &
                   axm, bxm, axd, bxd, cx, &
                   bx, dx, ax, gm, alf)

      use stel_kinds, only: dp
      use name0, only: cp9, c1p1
      use name1, only: mpol1, nmax, nsd, nsd1
      use scalars, only: ns
      use extfld, only: ivac
      use mnarray, only: jmin2
      use inputdat, only: nfp

      implicit none

      real(kind=dp), intent(inout) :: gcx(ns,0:nmax,0:mpol1,2) ! force Fourier coefficients to precondition

      ! preconditioner matrix elements from precondn()
      real(kind=dp), intent(in)    :: axm(nsd1,2)
      real(kind=dp), intent(in)    :: bxm(nsd1,2)
      real(kind=dp), intent(in)    :: axd(nsd1,2)
      real(kind=dp), intent(in)    :: bxd(nsd1,2)
      real(kind=dp), intent(in)    :: cx (nsd   )

      real(kind=dp), intent(out) :: bx(ns,0:nmax,0:mpol1) ! temporary arrays
      real(kind=dp), intent(out) :: dx(ns,0:nmax,0:mpol1)
      real(kind=dp), intent(out) :: ax(ns,0:nmax,0:mpol1)
      real(kind=dp), intent(out) :: gm(*)
      real(kind=dp), intent(out) :: alf(*)

      integer :: jmax, js, m, mp, n, ntype

      jmax = ns
      if (ivac.lt.1) &
        jmax = ns-1

      do m = 0, mpol1
        mp = mod(m, 2) + 1 ! m-parity: even m --> 1, odd-m --> 2
        do n = 0, nmax

          ! assemble preconditioning matrix elements
          do js = jmin2(m), jmax
            bx(js,n,m) = axm(js  ,mp) + bxm(js  ,mp)*m**2
            dx(js,n,m) = axd(js  ,mp) + bxd(js  ,mp)*m**2 + cx(js)*(n*nfp)**2
            ax(js,n,m) = axm(js+1,mp) + bxm(js+1,mp)*m**2
          end do

          if (m.eq.1) then
            ! TODO: what is this?
            ! m=1-component of ns=2 for all n
            ! maybe has to do with m=1-constraint for innermost flux surface at js=2?
            dx(2,n,1) = dx(2,n,1) + bx(2,n,1)
          end if

          ! TODO: What happens here ?
          ! corrections at LCFS ???
          ! only relevant in free-boundary mode?
          dx(ns,n,m) = c1p1*dx(ns,n,m)
          bx(ns,n,m) =  cp9*bx(ns,n,m)
        end do
      end do

      ! solve tri-diagonal system for cc/ss and cs/sc coefficients separately
      ! in order to precondition force matrices
      do ntype = 1, 2
        ! ax: super-diagonal
        ! dx:       diagonal
        ! bx: sub  -diagonal elements
        ! gcx: right-hand side and solution
        call trid(ax, dx, bx, gcx(1,0,0,ntype), gm, alf, jmin2, jmax)
      end do

      return
end
