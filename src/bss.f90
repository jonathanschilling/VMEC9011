subroutine bss(r12, rs, zs, ru12, zu12, bsubs, &
               lu, lv, br, bphi, bz)

      use stel_kinds, only: dp
      use name0, only: cp25, cp5
      use scalars, only: nrzt
      use realsp, only: r1, ru, rv, z1, zu, zv
      use scalefac, only: shalf

      implicit none

      real(kind=dp), intent(in)  :: r12  (nrzt)
      real(kind=dp), intent(in)  :: rs   (nrzt)
      real(kind=dp), intent(in)  :: zs   (nrzt)
      real(kind=dp), intent(in)  :: ru12 (nrzt)
      real(kind=dp), intent(in)  :: zu12 (nrzt)
      real(kind=dp), intent(in)  :: lu   (nrzt)
      real(kind=dp), intent(in)  :: lv   (nrzt)
      real(kind=dp), intent(out) :: bsubs(nrzt)
      real(kind=dp), intent(out) :: br   (nrzt)
      real(kind=dp), intent(out) :: bphi (nrzt)
      real(kind=dp), intent(out) :: bz   (nrzt)

      integer       :: lodd
      real(kind=dp) :: rv12, zv12
      real(kind=dp) :: gsu, gsv

      do l = 2, nrzt
        lodd = l + nrzt

        rv12 =     cp5 * ( rv(l) + rv(l-1) + shalf(l) * (rv(lodd)+rv(lodd-1)) )
        zv12 =     cp5 * ( zv(l) + zv(l-1) + shalf(l) * (zv(lodd)+zv(lodd-1)) )

        ! metric elements g_su, g_sv on half-grid
        gsu  =   rs(l) * ru12(l)                                       &
               + zs(l) * zu12(l)                                       &
               +  cp25 * (     r1(lodd)*ru(lodd)+r1(lodd-1)*ru(lodd-1) &
                          +    z1(lodd)*zu(lodd)+z1(lodd-1)*zu(lodd-1) &
                          + (  r1(lodd)*ru(l)   +r1(lodd-1)*ru(l-1)    &
                             + z1(lodd)*zu(l)   +z1(lodd-1)*zu(l-1) )/shalf(l) )

        gsv  =   rs(l) * rv12                                          &
               + zs(l) * zv12                                          &
               +  cp25 * (     r1(lodd)*rv(lodd)+r1(lodd-1)*rv(lodd-1) &
                          +    z1(lodd)*zv(lodd)+z1(lodd-1)*zv(lodd-1) &
                          + (  r1(lodd)*rv(l)   +r1(lodd-1)*rv(l-1)    &
                             + z1(lodd)*zv(l)   +z1(lodd-1)*zv(l-1) )/shalf(l) )

        bsubs(l)= lv(l)*gsu     + lu(l)*gsv

        br(l)   = lv(l)*ru12(l) + lu(l)*rv12
        bphi(l) =                 lu(l)*r12(l)
        bz(l)   = lv(l)*zu12(l) + lu(l)*zv12
      end do

      return
end
