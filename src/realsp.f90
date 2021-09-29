module realsp

      use stel_kinds, only: dp
      use name1,      only: nrztd

      implicit none

      real(kind=dp), target :: r1(2*nrztd)
      real(kind=dp)         :: ru(2*nrztd)
      real(kind=dp)         :: rv(2*nrztd)

      real(kind=dp), target :: z1(2*nrztd)
      real(kind=dp)         :: zu(2*nrztd)
      real(kind=dp)         :: zv(2*nrztd)

      real(kind=dp)         :: ru0(nrztd)
      real(kind=dp)         :: zu0(nrztd)

      real(kind=dp), target :: rcon(2*nrztd)
      real(kind=dp), target :: zcon(2*nrztd)
      real(kind=dp)         :: rcon0(nrztd)
      real(kind=dp)         :: zcon0(nrztd)

      real(kind=dp)         :: gcon(nrztd)
end module realsp
