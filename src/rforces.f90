module rforces

      use stel_kinds, only: dp
      use name1,      only: nrztd

      implicit none

      ! worka temporary array: 12*nrztd
      real(kind=dp), target :: armn(2*nrztd)
      real(kind=dp)         :: brmn(2*nrztd)
      real(kind=dp), target :: crmn(2*nrztd)
      real(kind=dp)         :: azmn(2*nrztd)
      real(kind=dp)         :: bzmn(2*nrztd)
      real(kind=dp), target :: czmn(2*nrztd)

      real(kind=dp)         :: blmn(2*nrztd)
      real(kind=dp)         :: clmn(2*nrztd)

end module rforces
