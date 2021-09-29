subroutine getiota(phipog, guu, guv, lu, lv)

      use stel_kinds, only: dp
      use name0, only: czero
      use name1, only: nznt
      use scalars, only: ns
      use inputdat, only: ncurr
      use current, only: jv
      use scalefac, only: wint
      use profs, only: iotas

      implicit none

      real(kind=dp), intent(in)    :: phipog(ns, nznt)
      real(kind=dp), intent(in)    :: guu   (ns, nznt)
      real(kind=dp), intent(in)    :: guv   (ns, nznt)
      real(kind=dp), intent(in)    :: lu    (ns, nznt)
      real(kind=dp), intent(inout) :: lv    (ns, nznt)

      integer       :: js, lk
      real(kind=dp) :: top, bot

      if (ncurr.ne.0) then
        do js = 2, ns

          top = jv(js)
          bot = czero

          do lk = 1, nznt
            top = top - wint(js,lk)*(  guu(js,lk)*lv(js,lk)  &
                                     + guv(js,lk)*lu(js,lk) )
            bot = bot + wint(js,lk)*   guu(js,lk)*phipog(js,lk)
          end do

          iotas(js) = top/bot

        end do
      end if

      ! ADD IOTA TO LAMBDA (ON HALF MESH NOW)
      do js = 2, ns
        do lk = 1, nznt
          lv(js,lk) = lv(js,lk) + phipog(js,lk)*iotas(js)
        end do
      end do

      return
end
