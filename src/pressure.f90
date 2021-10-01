subroutine pressure(gsqrt, bsq, wint)

      use stel_kinds, only: dp
      use name1, only: nznt
      use scalars, only: ns, dnorm
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
      end do

      do js = 2, ns
        vp(js)   = dnorm*abs(vp(js))
        pres(js) = mass(js)/vp(js)**gam
      end do

      wp = ddot(ns, vp, 1, pres, 1)/real(ns-1)

      do js = 2, ns
        do lk = 1, nznt
          bsq(js, lk) = pres(js)
        end do
      end do

      return
end
