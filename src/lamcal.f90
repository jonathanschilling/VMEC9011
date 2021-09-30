subroutine lamcal(phipog, guu, guv, gvv)

      use stel_kinds, only: dp
      use name0, only: cp5, c1p0, c2p0
      use name1, only: nsd1, nznt, mpol1, nmax, nmax1
      use scalars, only: nrzt, ns
      use inputdat, only: nfp
      use mnarray, only: jlam
      use scalefac, only: faclam
      use fsqu, only: mns

      implicit none

      ! TODO: more elegant definition of these BLAS functions
      real(kind=dp) :: dsum, ddot

      real(kind=dp), intent(in) :: phipog(nrzt)
      real(kind=dp), intent(in) :: guu(nrzt)
      real(kind=dp), intent(in) :: guv(nrzt)
      real(kind=dp), intent(in) :: gvv(nrzt)

      real(kind=dp) :: blam(nsd1)
      real(kind=dp) :: clam(nsd1)
      real(kind=dp) :: dlam(nsd1)
      data blam/nsd1*0.0_dp/, &
           clam/nsd1*0.0_dp/, &
           dlam/nsd1*0.0_dp/

      integer       :: js, l, m, n, lmn
      real(kind=dp) :: tnn, tmn, tmm

      do js = 2, ns
        blam(js) = ddot(nznt, guu(js), ns, phipog(js), ns)
        dlam(js) = ddot(nznt, guv(js), ns, phipog(js), ns)*c2p0*real(nfp)
        clam(js) = ddot(nznt, gvv(js), ns, phipog(js), ns)
      end do

      do m = 0, mpol1
        do n = 0, nmax
          lmn = ns*(n + nmax1*m)

          tnn = real( (n*nfp)**2 )
          if (m.eq.0 .and. n.eq.0) &
            tnn = -c1p0
          tmn = real(m*n)
          tmm = real(m*m)

          do js = jlam(m), ns
            faclam(js+lmn) = -c2p0*c2p0/(       (blam(js) + blam(js+1))            * tnn   &
                                         + sign((dlam(js) + dlam(js+1)), blam(js)) * tmn   &
                                         +      (clam(js) + clam(js+1))            * tmm )
          end do
          faclam(ns+lmn) = cp5*faclam(ns+lmn)
        end do
      end do

      do l = 1, mns
        faclam(l+mns) = faclam(l)
      end do

      return
end
