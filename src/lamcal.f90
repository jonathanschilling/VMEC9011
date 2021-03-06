subroutine lamcal(phipog, guu, guv, gvv)

      use stel_kinds, only: dp
      use name0, only: cp5, c1p0, c2p0
      use name1, only: nsd1, nznt, mpol1, nmax, nmax1
      use scalars, only: nrzt, ns, mns
      use inputdat, only: nfp
      use mnarray, only: jlam
      use scalefac, only: faclam

      implicit none

      ! TODO: more elegant definition of these BLAS functions
      real(kind=dp) :: dsum, ddot

      real(kind=dp), intent(in) :: phipog(nrzt) ! phi'/sqrt(g) * wint --> ready for integration; on half-grid
      real(kind=dp), intent(in) :: guu(nrzt)    !
      real(kind=dp), intent(in) :: guv(nrzt)    ! metric elements on half--grid
      real(kind=dp), intent(in) :: gvv(nrzt)    !

      real(kind=dp) :: blam(nsd1)
      real(kind=dp) :: clam(nsd1)
      real(kind=dp) :: dlam(nsd1)
      data blam/nsd1*0.0_dp/, & ! initial data for startup ??
           clam/nsd1*0.0_dp/, &
           dlam/nsd1*0.0_dp/

      integer       :: js, l, m, n, lmn
      real(kind=dp) :: tnn, tmn, tmm

      ! A(s), B(s), C(s) in Betancourt (1988) "BETAS", p. 559
      do js = 2, ns
        blam(js) = ddot(nznt, guu(js), ns, phipog(js), ns)
        dlam(js) = ddot(nznt, guv(js), ns, phipog(js), ns)*c2p0*real(nfp)
        clam(js) = ddot(nznt, gvv(js), ns, phipog(js), ns)
      end do

      do m = 0, mpol1
        do n = 0, nmax
          lmn = ns*(n + nmax1*m)

          tnn = real( (n*nfp)**2 )
          tmn = real(m*n)
          tmm = real(m*m)

          if (m.eq.0 .and. n.eq.0) then
            ! only a single term (by t_nn) contributes to faclam for m=0, n=0
            ! since t_mn and t_mm are zero
            ! --> -blam
            tnn = -c1p0
          end if


          ! note: sign(a,b) = |a| * sgn(b)
          ! note: 2/x = 1/(0.5*x)

          do js = jlam(m), ns
            faclam(js+lmn) = -c2p0*c2p0/(       (blam(js) + blam(js+1))            * tnn   &
                                         + sign((dlam(js) + dlam(js+1)), blam(js)) * tmn   &
                                         +      (clam(js) + clam(js+1))            * tmm )
          end do

          ! edge inter-/extrapolation correction at LCFS
          faclam(ns+lmn) = cp5*faclam(ns+lmn)
        end do
      end do

      ! same "lambda factors" for flsc and flcs
      do l = 1, mns
        faclam(l+mns) = faclam(l)
      end do

      return
end
