module rforces

      use stel_kinds, only: dp
      use name1,      only: nrztd

      implicit none

      real(kind=dp), target  :: worka(12*nrztd)

      real(kind=dp), pointer :: armn(:)
      real(kind=dp), pointer :: brmn(:)
      real(kind=dp), pointer :: crmn(:)
      real(kind=dp), pointer :: azmn(:)
      real(kind=dp), pointer :: bzmn(:)
      real(kind=dp), pointer :: czmn(:)

      real(kind=dp) :: blmn(2*nrztd)
      real(kind=dp) :: clmn(2*nrztd)

      contains

subroutine setup_rforces
      implicit none
      armn => worka(1+ 0*nrztd: 2*nrztd)
      brmn => worka(1+ 2*nrztd: 4*nrztd)
      crmn => worka(1+ 4*nrztd: 6*nrztd)
      azmn => worka(1+ 6*nrztd: 8*nrztd)
      bzmn => worka(1+ 8*nrztd:10*nrztd)
      czmn => worka(1+10*nrztd:12*nrztd)
end subroutine

end module ! rforces
