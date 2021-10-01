subroutine spectrum(rcc, rss, zcs, zsc)

      use stel_kinds, only: dp
      use name0, only: czero
      use name1, only: mpol1, nmax, nsd
      use scalars, only: ns
      use mnarray, only: mscale, nscale, xmpq
      use spectra, only: specw

      implicit none

      real(kind=dp), intent(in) :: rcc(ns,0:nmax,0:mpol1)
      real(kind=dp), intent(in) :: rss(ns,0:nmax,0:mpol1)
      real(kind=dp), intent(in) :: zsc(ns,0:nmax,0:mpol1)
      real(kind=dp), intent(in) :: zcs(ns,0:nmax,0:mpol1)

      real(kind=dp) :: t1(nsd)
      real(kind=dp) :: dnumer(nsd)
      real(kind=dp) :: denom (nsd)

      integer       :: js, n, m
      real(kind=dp) :: scalefac

      do js = 2, ns
        dnumer(js) = czero
        denom (js) = czero
      end do

      do n = 0, nmax
        do m = 1, mpol1

          scalefac = (mscale(m)*nscale(n))**2

          do js = 2, ns
            t1(js) = (  rcc(js,n,m)**2 + rss(js,n,m)**2 &
                      + zcs(js,n,m)**2 + zsc(js,n,m)**2  ) * scalefac
          end do

          do js = 2, ns
            dnumer(js) = dnumer(js) + t1(js)*xmpq(m,3)
            denom (js) = denom (js) + t1(js)*xmpq(m,2)
          end do

        end do
      end do

      do js = 2, ns
        specw(js) = dnumer(js)/denom(js)
      end do

      return
end
