subroutine alias(gcon, zcon, work1, work2, work3, gcs, gsc)

      use stel_kinds, only: dp

      implicit none

      real(kind=dp) :: gcon(ns*nzeta,ntheta2)
      real(kind=dp) :: zcon(ns*nzeta,ntheta2)
      real(kind=dp) :: work1(4*ns*nzeta)
      real(kind=dp) :: work2(ns*nzeta,4)
      real(kind=dp) :: work3(ns,nzeta,4)
      real(kind=dp) :: gcs(ns,0:nmax,0:mpol1)
      real(kind=dp) :: gsc(ns,0:nmax,0:mpol1)

      ! BEGIN DE-ALIASING (TRUNCATION OF GCON IN FOURIER-SPACE
      do m = 1, mpol1-1 ! --> exclude m=0, mpol-2 (lowest and highest)

        ! clean workspace
        do l = 1,4*ns*nzeta
          work1(l) = czero
        end do

        do i = 1,ntheta2
          do jk = 1,ns*nzeta
            work2(jk,01) = work2(jk,01) + zcon(jk,i)*cosmui(i,m)
            work2(jk,02) = work2(jk,02) + zcon(jk,i)*sinmu (i,m)
          end do
        end do

        do n = 0,nmax
          fm = faccon(n,m)
          do k = 1,nzeta
            do js= 2,ns
              gcs(js,n,m) =gcs(js,n,m) +fm*tcon(js)*work3(js,k,01)*sinnv(k,n)
              gsc(js,n,m) =gsc(js,n,m) +fm*tcon(js)*work3(js,k,02)*cosnv(k,n)
            end do
          end do
        end do

        !  RECONSTRUCT DE-ALIASED GCON
        do n = 0,nmax
          do k = 1,nzeta
            do js= 2,ns
              work3(js,k,03) = work3(js,k,03) + gcs(js,n,m)*sinnv(k,n)
              work3(js,k,04) = work3(js,k,04) + gsc(js,n,m)*cosnv(k,n)
            end do
          end do
        end do
        do i = 1,ntheta2
          do jk= 1,ns*nzeta
            gcon(jk,i) = gcon(jk,i) + work2(jk,03)*cosmu(i,m) &
                                    + work2(jk,04)*sinmu(i,m)
          end do
        end do

      end do ! loop over m = 1, mpol1-1

      return

end
