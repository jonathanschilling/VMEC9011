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

      real(kind=dp) :: gcx(ns,0:nmax,0:mpol1,2)
      real(kind=dp) :: axm(nsd1,2)
      real(kind=dp) :: bxm(nsd1,2)
      real(kind=dp) :: axd(nsd1,2)
      real(kind=dp) :: bxd(nsd1,2)
      real(kind=dp) :: cx(nsd)
      real(kind=dp) :: bx(ns,0:nmax,0:mpol1)
      real(kind=dp) :: dx(ns,0:nmax,0:mpol1)
      real(kind=dp) :: ax(ns,0:nmax,0:mpol1)
      real(kind=dp) :: gm(*)
      real(kind=dp) :: alf(*)

      integer :: jmax, js, m, mp, n, ntype

      jmax = ns
      if (ivac.lt.1) &
        jmax = ns-1

      do m = 0, mpol1
        mp = mod(m, 2) + 1
        do n = 0, nmax

          do js = jmin2(m), jmax
            bx(js,n,m) = axm(js  ,mp) + bxm(js  ,mp)*m**2
            dx(js,n,m) = axd(js  ,mp) + bxd(js  ,mp)*m**2 + cx(js)*(n*nfp)**2
            ax(js,n,m) = axm(js+1,mp) + bxm(js+1,mp)*m**2
          end do

          if (m.eq.1) &
            dx(2,n,1) = dx(2,n,1) + bx(2,n,1)

          ! TODO: What happens here ?
          dx(ns,n,m) = c1p1*dx(ns,n,m)
          bx(ns,n,m) =  cp9*bx(ns,n,m)
        end do
      end do

      do ntype = 1, 2
        call trid(ax, dx, bx, gcx(1,0,0,ntype), gm, alf, jmin2, jmax)
      end do

      return
end
