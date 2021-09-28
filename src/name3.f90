module realsp
      use stel_kinds, only: dp
      implicit none
      real(kind=dp) :: r1(2*nrztd)
      real(kind=dp) :: ru(2*nrztd)
      real(kind=dp) :: rv(2*nrztd)
      real(kind=dp) :: z1(2*nrztd)
      real(kind=dp) :: zu(2*nrztd)
      real(kind=dp) :: zv(2*nrztd)
      real(kind=dp) :: rcon(2*nrztd)
      real(kind=dp) :: zcon(2*nrztd)
      real(kind=dp) :: gcon(nrztd)
      real(kind=dp) :: ru0(nrztd)
      real(kind=dp) :: zu0(nrztd)
      real(kind=dp) :: rcon0(nrztd)
      real(kind=dp) :: zcon0(nrztd)
end module realsp
