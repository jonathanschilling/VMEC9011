subroutine pressure(gsqrt, bsq, wint)

      use stel_kinds, only: dp
      use name1, only: nznt
      use scalars, only: ns, hs, dnorm
      use profs, only: vp, mass, pres
      use inputdat, only: gam
      use fsqu, only: wp

      implicit none

      ! TODO: more elegant definition of these BLAS functions
      real(kind=dp) :: ddot

      real(kind=dp), intent(in)  :: gsqrt(ns,nznt)
      real(kind=dp), intent(out) :: bsq  (ns,nznt)
      real(kind=dp), intent(in)  :: wint (ns,nznt)

      integer :: js, lk

      do js = 2, ns
        vp(js) = ddot(nznt, gsqrt(js,1), ns, wint(js,1), ns)
        vp(js) = dnorm*abs(vp(js))
      end do

      do js = 2, ns
        pres(js) = mass(js)/vp(js)**gam
      end do

      wp = hs*ddot(ns, vp, 1, pres, 1)

      ! bsq gets a scalar offset of the kinetic pressure on each flux surface
      do js = 2, ns
        do lk = 1, nznt
          bsq(js, lk) = pres(js)
        end do
      end do

      return
end
