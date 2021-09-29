function pmass(x)
      use stel_kinds, only: dp
      use inputdat, only: am
      implicit none

      real(kind=dp) :: pmass
      real(kind=dp) :: x

      pmass = am(0)+x*(am(1)+x*(am(2)+x*(am(3)+x*(am(4)+x*am(5)))))

end function

function piota(x)
      use stel_kinds, only: dp
      use inputdat, only: ai
      implicit none

      real(kind=dp) :: piota
      real(kind=dp) :: x

      piota = ai(0)+x*(ai(1)+x*(ai(2)+x*(ai(3)+x*(ai(4)+x*ai(5)))))

end function

function pcurr(x)
      use stel_kinds, only: dp
      use inputdat, only: ac
      implicit none

      real(kind=dp) :: pcurr
      real(kind=dp) :: x

      pcurr = x*(ac(1)+x*(ac(2)+x*(ac(3)+x*(ac(4)+x*ac(5)))))

end function
