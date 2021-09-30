module realsp

      use stel_kinds, only: dp
      use name1,      only: nrztd

      implicit none

      real(kind=dp), target  :: workb(16*nrzd)

      real(kind=dp), pointer ::   r1(2*nrztd) => workb(1+ 0*nzrd)
      real(kind=dp), pointer ::   ru(2*nrztd) => workb(1+ 2*nzrd)
      real(kind=dp), pointer ::   rv(2*nrztd) => workb(1+ 4*nzrd)
      real(kind=dp), pointer ::   z1(2*nrztd) => workb(1+ 6*nzrd)
      real(kind=dp), pointer ::   zu(2*nrztd) => workb(1+ 8*nzrd)
      real(kind=dp), pointer ::   zv(2*nrztd) => workb(1+10*nzrd)
      real(kind=dp), pointer :: rcon(2*nrztd) => workb(1+12*nzrd)
      real(kind=dp), pointer :: zcon(2*nrztd) => workb(1+14*nzrd)

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
        zv(l)   = czero ! 12*nrzd up to here
        rcon(l) = czero
        zcon(l) = czero ! 16*nrzd up to here
      end do

end subroutine

end module realsp
