module realsp

      use stel_kinds, only: dp
      use name1,      only: nrztd

      implicit none

      real(kind=dp) :: workb(16*nrztd)

      ! 2*nrztd per array
      real(kind=dp) ::   r1(2*nrztd)
      real(kind=dp) ::   ru(2*nrztd)
      real(kind=dp) ::   rv(2*nrztd)
      real(kind=dp) ::   z1(2*nrztd)
      real(kind=dp) ::   zu(2*nrztd)
      real(kind=dp) ::   zv(2*nrztd)
      real(kind=dp) :: rcon(2*nrztd)
      real(kind=dp) :: zcon(2*nrztd)
      equivalence (  r1, workb(1+ 0*nrztd)), &
                  (  ru, workb(1+ 2*nrztd)), &
                  (  rv, workb(1+ 4*nrztd)), &
                  (  z1, workb(1+ 6*nrztd)), &
                  (  zu, workb(1+ 8*nrztd)), &
                  (  zv, workb(1+10*nrztd)), &
                  (rcon, workb(1+12*nrztd)), &
                  (zcon, workb(1+14*nrztd))

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
