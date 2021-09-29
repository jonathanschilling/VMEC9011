module realsp

      use stel_kinds, only: dp
      use name1,      only: nrztd

      implicit none

      ! 2*... --> split into contributions
      ! from even-m and odd-m poloidal modes
      real(kind=dp), target ::   r1(2*nrztd)
      real(kind=dp)         ::   ru(2*nrztd)
      real(kind=dp)         ::   rv(2*nrztd)
      real(kind=dp), target ::   z1(2*nrztd)
      real(kind=dp)         ::   zu(2*nrztd)
      real(kind=dp)         ::   zv(2*nrztd)
      real(kind=dp), target :: rcon(2*nrztd)
      real(kind=dp), target :: zcon(2*nrztd)

      real(kind=dp)         :: ru0(nrztd)
      real(kind=dp)         :: zu0(nrztd)

      real(kind=dp)         :: rcon0(nrztd)
      real(kind=dp)         :: zcon0(nrztd)

      real(kind=dp)         :: gcon(nrztd)

      contains

subroutine clear_realsp

      use name0, only: czero
      use name1, only: nrztd

      implicit none

      integer :: l

      do l = 1, 2*nrztd
        r1(l)   = czero
        ru(l)   = czero
        rv(l)   = czero
        z1(l)   = czero
        zu(l)   = czero
        zv(l)   = czero
        rcon(l) = czero
        zcon(l) = czero
      end do

end subroutine

end module realsp
