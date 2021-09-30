module realsp

      use stel_kinds, only: dp
      use name1,      only: nrztd

      implicit none

      real(kind=dp), target  :: workb(16*nrztd)

      real(kind=dp), pointer ::   r1(:) => workb(1+ 0*nrztd: 2*nrztd)
      real(kind=dp), pointer ::   ru(:) => workb(1+ 2*nrztd: 4*nrztd)
      real(kind=dp), pointer ::   rv(:) => workb(1+ 4*nrztd: 6*nrztd)
      real(kind=dp), pointer ::   z1(:) => workb(1+ 6*nrztd: 8*nrztd)
      real(kind=dp), pointer ::   zu(:) => workb(1+ 8*nrztd:10*nrztd)
      real(kind=dp), pointer ::   zv(:) => workb(1+10*nrztd:12*nrztd)
      real(kind=dp), pointer :: rcon(:) => workb(1+12*nrztd:14*nrztd)
      real(kind=dp), pointer :: zcon(:) => workb(1+14*nrztd:16*nrztd)

      real(kind=dp) :: ru0(nrztd)
      real(kind=dp) :: zu0(nrztd)

      real(kind=dp) :: rcon0(nrztd)
      real(kind=dp) :: zcon0(nrztd)

      real(kind=dp) :: gcon(nrztd)

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
        zv(l)   = czero ! 12*nrztd up to here
        rcon(l) = czero
        zcon(l) = czero ! 16*nrztd up to here
      end do

end subroutine ! clear_realsp

end module ! realsp
