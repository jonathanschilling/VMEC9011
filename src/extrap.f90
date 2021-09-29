subroutine extrap(rmncc, rmnss, zmncs, zmnsc, lmncs, lmnsc, x3, x4, ns)

      use stel_kinds, only: dp
      use name1, only: mnd1, nmax, nmax1

      implicit none

      real(kind=dp), intent(inout) :: rmncc(ns,0:mnd1)
      real(kind=dp), intent(inout) :: rmnss(ns,0:mnd1)
      real(kind=dp), intent(in)    ::       x3(0:mnd1)
      real(kind=dp), intent(inout) :: zmncs(ns,0:mnd1)
      real(kind=dp), intent(inout) :: zmnsc(ns,0:mnd1)
      real(kind=dp), intent(in)    ::       x4(0:mnd1)
      real(kind=dp), intent(inout) :: lmncs(ns,0:mnd1)
      real(kind=dp), intent(inout) :: lmnsc(ns,0:mnd1)
      integer,       intent(in)    ::       ns

      integer :: mn, n

      do mn = 2*nmax1, mnd1
        rmncc(2,mn) = x3(mn)*rmncc(3,mn) + x4(mn)*rmncc(4,mn)
        rmnss(2,mn) = x3(mn)*rmnss(3,mn) + x4(mn)*rmnss(4,mn)
        zmncs(2,mn) = x3(mn)*zmncs(3,mn) + x4(mn)*zmncs(4,mn)
        zmnsc(2,mn) = x3(mn)*zmnsc(3,mn) + x4(mn)*zmnsc(4,mn)
      enddo

      do mn = nmax1, mnd1
        lmncs(2,mn) = x3(mn)*lmncs(3,mn) + x4(mn)*lmncs(4,mn)
        lmnsc(2,mn) = x3(mn)*lmnsc(3,mn) + x4(mn)*lmnsc(4,mn)
      enddo

      do n = 1,nmax
        lmncs(1,n) = x3(0)*lmncs(2,n) - lmncs(3,n)
      enddo

      return
end
