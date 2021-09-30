subroutine interp(xnew, xold, scalefac, nsin)

      use stel_kinds, only: dp
      use name0, only: czero, c1p0, c2p0, c1pm10
      use name1, only: mpol1, nmax, mnd1
      use scalars, only: ns, hs

      implicit none

      real(kind=dp), intent(out)   ::  xnew(ns,  0:nmax,0:mpol1,6)
      real(kind=dp), intent(inout) ::  xold(nsin,0:nmax,0:mpol1,6)
      real(kind=dp), intent(in)    :: scalefac(ns,  0:nmax,0:mpol1  )
      integer,       intent(in)    :: nsin

      integer       :: ntype, js, js1, js2, n, m
      real(kind=dp) :: hsold, sj, s1, xint

      hsold = c1p0/real(nsin-1)

      do ntype = 1, 6

        ! extrapolate to axis before interpolation
        do n = 0, nmax
          xold(1,n,1,ntype) = c2p0*xold(2,n,1,ntype) - xold(3,n,1,ntype)
        end do

        ! INTERPOLATE R,Z AND LAMBDA ON FULL GRID
        do js = 1, ns

          sj = (js-1)*hs
          js1 = 1 + int(sj/hsold + c1pm10) ! TODO: why c1pm10 ?
          js2 = min(js1+1, nsin)
          s1 = (js1-1)*hsold
          xint = (sj-s1)/hsold

          do n = 0, nmax
            do m = 0, mpol1
              xnew(js,n,m,ntype) = (  (c1p0-xint)*xold(js1,n,m,ntype) &
                                    +       xint *xold(js2,n,m,ntype) ) / scalefac(js,n,m)
            end do ! m
          end do ! n
        end do ! js

        ! set interpolated axis values to zero
        do n = 0, nmax
          xnew(1,n,1,ntype) = czero
        end do ! n

      end do ! ntype

      return
end

