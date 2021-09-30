subroutine getfsq(gcr, gcz, gnormr, gnormz, gnorm, mprecon)

      use stel_kinds, only: dp
      use name0, only: czero
      use name1, only: mnd1
      use scalars, only: ns

      implicit none

      real(kind=dp), intent(in)  :: gcr(ns,0:mnd1,2)
      real(kind=dp), intent(in)  :: gcz(ns,0:mnd1,2)
      real(kind=dp), intent(out) :: gnormr
      real(kind=dp), intent(out) :: gnormz
      real(kind=dp), intent(in)  :: gnorm
      integer,       intent(in)  :: mprecon

      integer :: jsmax, mn, js

      gnormr = czero
      gnormz = czero

      jsmax = (ns-1) + mprecon
      do mn = 0, mnd1
        do js = 1, jsmax
          gnormr = gnormr + gnorm * (gcr(js,mn,1)**2 + gcr(js,mn,2)**2)
          gnormz = gnormz + gnorm * (gcz(js,mn,1)**2 + gcz(js,mn,2)**2)
        end do
      end do

      return
end
