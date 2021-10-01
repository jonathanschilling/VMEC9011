subroutine trid(a, d, b, c, gam, alf, jmin, nn)

      use stel_kinds, only: dp
      use name1, only: mnd1
      use bounds, only: mlower, mupper
      use scalars, only: ns

      implicit none

      real(kind=dp), intent(in)    :: a  (ns,0:mnd1)
      real(kind=dp), intent(in)    :: b  (ns,0:mnd1)
      real(kind=dp), intent(inout) :: c  (ns,0:mnd1)
      real(kind=dp), intent(in)    :: d  (ns,0:mnd1)
      real(kind=dp), intent(out)   :: gam(0:mnd1,ns) ! only internal temporary storage
      real(kind=dp), intent(out)   :: alf(0:mnd1,ns) ! only internal temporary storage
      integer,       intent(in)    :: jmin(0:3)
      integer,       intent(in)    :: nn

      integer :: imodes, in, in1, mn, i, i1, n2

      ! SOLVES B(I)*X(I-1)+D(I)*X(I)+A(I)*X(I+1)=C(I), I=IN,NN
      ! AND RETURNS ANSWER IN C(I)

      ! SEPARATE M=0 (IMODES=1), M=1 (IMODES=2), M>1 (IMODES=3)
      do imodes = 1, 3
        in  = jmin(imodes-1)
        in1 = in + 1

        do mn = mlower(imodes), mupper(imodes)
          gam(mn,in) = d(in,mn)
        end do

        do i = in1, nn
          do mn = mlower(imodes), mupper(imodes)
            alf(mn,i-1) = a(i-1,mn)/gam(mn,i-1)
            gam(mn,i  ) = d(i  ,mn) - b(i,mn)*alf(mn,i-1)
          end do
        end do

        do mn = mlower(imodes), mupper(imodes)
          c(in,mn) = c(in,mn)/gam(mn,in)
        end do

        do i = in1, nn
          do mn = mlower(imodes), mupper(imodes)
            c(i,mn) = (c(i,mn) - b(i,mn)*c(i-1,mn))/gam(mn,i)
          end do
        end do

        n2 = nn + in
        do i = in1, nn
          i1 = n2-i

          do mn = mlower(imodes), mupper(imodes)
            c(i1,mn) = c(i1,mn) - alf(mn,i1)*c(i1+1,mn)
          end do
        end do

      end do

      return
end
