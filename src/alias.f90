subroutine alias(gcon, zcon, work1, work2, work3, gcs, gsc)

      use stel_kinds, only: dp
      use name0, only: czero
      use name1, only: nzeta, ntheta2, nmax, mpol1
      use scalars, only: ns
      use trignew, only: cosmu, sinmu, cosmui, cosnv, sinnv
      use spectra, only: faccon, tcon

      implicit none

      real(kind=dp), intent(inout) :: gcon(ns*nzeta,ntheta2) ! in: zeroed; out: de-aliased gcon
      real(kind=dp), intent(in)    :: zcon(ns*nzeta,ntheta2) ! rcon*ru + zcon*zu
      real(kind=dp), intent(out)   :: work1(4*ns*nzeta)
      real(kind=dp), intent(out)   :: work2(ns*nzeta,4)
      real(kind=dp), intent(out)   :: work3(ns,nzeta,4)
      real(kind=dp), intent(inout) :: gcs(ns,0:nmax,0:mpol1) ! temp
      real(kind=dp), intent(inout) :: gsc(ns,0:nmax,0:mpol1) ! temp

      integer       :: js, jk, k, l, m, n
      real(kind=dp) :: fm

      ! BEGIN DE-ALIASING (TRUNCATION OF GCON IN FOURIER-SPACE
      do m = 1, mpol1-1 ! --> exclude m=0, mpol-2 (lowest and highest)

        ! clean workspace
        do l = 1,4*ns*nzeta
          work1(l) = czero
        end do

        do i = 1,ntheta2
          do jk = 1,ns*nzeta
            work2(jk,1) = work2(jk,1) + zcon(jk,i)*cosmui(i,m)
            work2(jk,2) = work2(jk,2) + zcon(jk,i)*sinmu (i,m)
          end do
        end do

        do n = 0,nmax
          fm = faccon(n,m)
          do k = 1,nzeta
            do js= 2,ns
              gcs(js,n,m) = gcs(js,n,m) + fm*tcon(js) * work3(js,k,1)*sinnv(k,n)
              gsc(js,n,m) = gsc(js,n,m) + fm*tcon(js) * work3(js,k,2)*cosnv(k,n)
            end do
          end do
        end do

        !  RECONSTRUCT DE-ALIASED GCON
        do n = 0,nmax
          do k = 1,nzeta
            do js= 2,ns
              work3(js,k,3) = work3(js,k,3) + gcs(js,n,m)*sinnv(k,n)
              work3(js,k,4) = work3(js,k,4) + gsc(js,n,m)*cosnv(k,n)
            end do
          end do
        end do
        do i = 1,ntheta2
          do jk= 1,ns*nzeta
            gcon(jk,i) = gcon(jk,i) + work2(jk,3)*cosmu(i,m) &
                                    + work2(jk,4)*sinmu(i,m)
          end do
        end do

      end do ! loop over m = 1, mpol1-1

      return

end
