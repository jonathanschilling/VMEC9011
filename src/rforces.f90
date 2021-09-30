module rforces

      use stel_kinds, only: dp
      use name1,      only: nrztd

      implicit none

      real(kind=dp), target  :: worka(12*nrztd)

      real(kind=dp), pointer :: armn(2*nrztd) => worka(1+ 0*nrzd)
      real(kind=dp), pointer :: brmn(2*nrztd) => worka(1+ 2*nrzd)
      real(kind=dp), pointer :: crmn(2*nrztd) => worka(1+ 4*nrzd)
      real(kind=dp), pointer :: azmn(2*nrztd) => worka(1+ 6*nrzd)
      real(kind=dp), pointer :: bzmn(2*nrztd) => worka(1+ 8*nrzd)
      real(kind=dp), pointer :: czmn(2*nrztd) => worka(1+10*nrzd)

      real(kind=dp)          :: blmn(2*nrztd)
      real(kind=dp)          :: clmn(2*nrztd)

end module rforces
