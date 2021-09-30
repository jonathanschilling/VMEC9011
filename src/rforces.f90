module rforces

      use stel_kinds, only: dp
      use name1,      only: nrztd

      implicit none

      real(kind=dp), target  :: worka(12*nrztd)

      real(kind=dp), pointer :: armn(:) => worka(1+ 0*nrztd: 2*nrztd)
      real(kind=dp), pointer :: brmn(:) => worka(1+ 2*nrztd: 4*nrztd)
      real(kind=dp), pointer :: crmn(:) => worka(1+ 4*nrztd: 6*nrztd)
      real(kind=dp), pointer :: azmn(:) => worka(1+ 6*nrztd: 8*nrztd)
      real(kind=dp), pointer :: bzmn(:) => worka(1+ 8*nrztd:10*nrztd)
      real(kind=dp), pointer :: czmn(:) => worka(1+10*nrztd:12*nrztd)

      real(kind=dp)          :: blmn(2*nrztd)
      real(kind=dp)          :: clmn(2*nrztd)

end module rforces
