subroutine tomnsp(frcc, frss, fzcs, fzsc, flcs, flsc, &
                  armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn, arcon, azcon, &
                  work1, work2, work3)

      use stel_kinds, only: dp
      use name0, only: czero
      use name1, only: mpol1, nmax, nzeta, ntheta2
      use scalars, only: ns
      use extfld, only: ivac
      use mnarray, only: xmpq, jmin2, jlam
      use trignew, only: cosmui, sinmu, cosmumi, sinmum, &
                         cosnv,  sinnv, cosnvn,  sinnvn

      implicit none

      ! inout needed since the force vector gc is initialized to 0 outside (funct3d)
      real(kind=dp), intent(inout) :: frcc(ns,0:nmax,0:mpol1)
      real(kind=dp), intent(inout) :: frss(ns,0:nmax,0:mpol1)
      real(kind=dp), intent(inout) :: fzcs(ns,0:nmax,0:mpol1)
      real(kind=dp), intent(inout) :: fzsc(ns,0:nmax,0:mpol1)
      real(kind=dp), intent(inout) :: flcs(ns,0:nmax,0:mpol1)
      real(kind=dp), intent(inout) :: flsc(ns,0:nmax,0:mpol1)

      real(kind=dp), intent(in)    :: armn (ns*nzeta,ntheta2,0:1)
      real(kind=dp), intent(in)    :: brmn (ns*nzeta,ntheta2,0:1)
      real(kind=dp), intent(in)    :: crmn (ns*nzeta,ntheta2,0:1)
      real(kind=dp), intent(in)    :: azmn (ns*nzeta,ntheta2,0:1)
      real(kind=dp), intent(in)    :: bzmn (ns*nzeta,ntheta2,0:1)
      real(kind=dp), intent(in)    :: czmn (ns*nzeta,ntheta2,0:1)
      real(kind=dp), intent(in)    :: blmn (ns*nzeta,ntheta2,0:1)
      real(kind=dp), intent(in)    :: clmn (ns*nzeta,ntheta2,0:1)
      real(kind=dp), intent(in)    :: arcon(ns*nzeta,ntheta2,0:1)
      real(kind=dp), intent(in)    :: azcon(ns*nzeta,ntheta2,0:1)

      ! all these refer to the global "workb" array
      real(kind=dp), intent(inout) :: work1(ns*nzeta*12)
      real(kind=dp), intent(inout) :: work2(ns*nzeta,12)
      real(kind=dp), intent(inout) :: work3(ns,nzeta,12)

      integer       :: jmax, i, jk, k, l, m, n, js, mparity
      real(kind=dp) :: temp1, temp3

      ! include contribution from LCFS only in free-boundary mode
      jmax = ns
      if (ivac.lt.1) &
        jmax = ns-1

      ! BEGIN INVERSE FOURIER TRANSFORM
      ! DO THETA (U) INTEGRATION FIRST
      do m = 0, mpol1
        mparity = mod(m, 2)

        ! clean up workspace first
        do l = 1, 12*ns*nzeta
          work1(l) = czero
        end do

        do i = 1, ntheta2
          do jk = 1, ns*nzeta

            ! effective A forces on R and Z: MHD + Spectral Condensation
            temp1        =                armn(jk,i,mparity) + xmpq(m,1)*arcon(jk,i,mparity)
            temp3        =                azmn(jk,i,mparity) + xmpq(m,1)*azcon(jk,i,mparity)

            work2(jk, 1) = work2(jk, 1) + temp1             *cosmui (i,m) + brmn(jk,i,mparity)*sinmum (i,m)
            work2(jk, 2) = work2(jk, 2) - crmn(jk,i,mparity)*cosmui (i,m)
            work2(jk, 3) = work2(jk, 3) + temp1             *sinmu  (i,m) + brmn(jk,i,mparity)*cosmumi(i,m)
            work2(jk, 4) = work2(jk, 4) - crmn(jk,i,mparity)*sinmu  (i,m)
            work2(jk, 5) = work2(jk, 5) + temp3             *cosmui (i,m) + bzmn(jk,i,mparity)*sinmum (i,m)
            work2(jk, 6) = work2(jk, 6) - czmn(jk,i,mparity)*cosmui (i,m)
            work2(jk, 7) = work2(jk, 7) + temp3             *sinmu  (i,m) + bzmn(jk,i,mparity)*cosmumi(i,m)
            work2(jk, 8) = work2(jk, 8) - czmn(jk,i,mparity)*sinmu  (i,m)
            work2(jk, 9) = work2(jk, 9) + blmn(jk,i,mparity)*sinmum (i,m)
            work2(jk,10) = work2(jk,10) - clmn(jk,i,mparity)*cosmui (i,m)
            work2(jk,11) = work2(jk,11) + blmn(jk,i,mparity)*cosmumi(i,m)
            work2(jk,12) = work2(jk,12) - clmn(jk,i,mparity)*sinmu  (i,m)
          end do ! jk
        end do ! i

        ! NEXT, DO ZETA (V) INTEGRATION
        do n = 0, nmax

          ! ns*nzeta from first half of transform is split up here since
          ! R,Z forces are computed from jmin2(m) on and lambda forces from jlam(m) on
          ! --> for axis, inner most surface: forces only from restricted poloidal modes
          do k = 1, nzeta
            do js = jmin2(m), jmax
              frcc(js,n,m) = frcc(js,n,m) + work3(js,k, 1)*cosnv (k,n) &
                                          + work3(js,k, 2)*sinnvn(k,n)
              frss(js,n,m) = frss(js,n,m) + work3(js,k, 3)*sinnv (k,n) &
                                          + work3(js,k, 4)*cosnvn(k,n)
              fzcs(js,n,m) = fzcs(js,n,m) + work3(js,k, 5)*sinnv (k,n) &
                                          + work3(js,k, 6)*cosnvn(k,n)
              fzsc(js,n,m) = fzsc(js,n,m) + work3(js,k, 7)*cosnv (k,n) &
                                          + work3(js,k, 8)*sinnvn(k,n)
            end do ! js

            do js = jlam(m), ns
              flcs(js,n,m) = flcs(js,n,m) + work3(js,k, 9)*sinnv (k,n) &
                                          + work3(js,k,10)*cosnvn(k,n)
              flsc(js,n,m) = flsc(js,n,m) + work3(js,k,11)*cosnv (k,n) &
                                          + work3(js,k,12)*sinnvn(k,n)
            end do ! js
          end do ! k
        end do ! n

      end do ! m

      return
end
