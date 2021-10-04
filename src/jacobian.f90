subroutine jacobian(r1, ru, z1, zu, &
                    zu12, ru12, zs, rs, gsqrt, r12, tau)

      use stel_kinds, only: dp
      use name0, only: czero, cp25, cp5
      use scalars, only: nrzt, meven, modd, ohs, irst
      use scalefac, only: shalf, wint

      implicit none

      real(kind=dp), intent(in)  :: r1   (nrzt,0:1)
      real(kind=dp), intent(in)  :: ru   (nrzt,0:1)
      real(kind=dp), intent(in)  :: z1   (nrzt,0:1)
      real(kind=dp), intent(in)  :: zu   (nrzt,0:1)
      real(kind=dp), intent(out) :: ru12 (nrzt)
      real(kind=dp), intent(out) :: zu12 (nrzt)
      real(kind=dp), intent(out) :: rs   (nrzt)
      real(kind=dp), intent(out) :: zs   (nrzt)
      real(kind=dp), intent(out) :: r12  (nrzt)
      real(kind=dp), intent(out) :: gsqrt(nrzt)
      real(kind=dp), intent(out) :: tau  (nrzt)

      integer       :: l, lm
      real(kind=dp) :: taumin, taumax

      ! (RS, ZS)=(R, Z) SUB S, (RU12, ZU12)=(R, Z) SUB THETA(=U) AND GSQRT=SQRT(G)
      ! ARE DIFFERENCED ON HALF MESH
      do l = 2, nrzt
        lm = l-1

        ! dr/du, dz/du on half grid
        ru12(l) = cp5*(ru(l, meven) + ru(lm, meven) + shalf(l)*(ru(l, modd) + ru(lm, modd)))
        zu12(l) = cp5*(zu(l, meven) + zu(lm, meven) + shalf(l)*(zu(l, modd) + zu(lm, modd)))

        ! dr/ds, dz/ds on half grid
        rs(l)   = ohs*(r1(l, meven) - r1(lm, meven) + shalf(l)*(r1(l, modd) - r1(lm, modd)))
        zs(l)   = ohs*(z1(l, meven) - z1(lm, meven) + shalf(l)*(z1(l, modd) - z1(lm, modd)))

        ! r1 (=R) on half grid
        r12(l)  = cp5*(r1(l, meven) + r1(lm, meven) + shalf(l)*(r1(l, modd) + r1(lm, modd)))

        ! Jacobian
        gsqrt(l)= r12(l) * (  ru12(l)*zs(l) - rs(l)*zu12(l)                                                &
                            + cp25 * (     ru(l,modd )*z1(l,modd) + ru(lm,modd )*z1(lm,modd)               &
                                         - zu(l,modd )*r1(l,modd) - zu(lm,modd )*r1(lm,modd)               &
                                      + (  ru(l,meven)*z1(l,modd) + ru(lm,meven)*z1(lm,modd)               &
                                         - zu(l,meven)*r1(l,modd) - zu(lm,meven)*r1(lm,modd))/shalf(l) ) )

        ! TODO; in newer VMEC: gsqrt = R*tau
        tau(l) = wint(l)*gsqrt(l)
      end do

      ! TEST FOR SIGN CHANGE IN JACOBIAN
      taumin = czero
      taumax = czero
      do l = 2, nrzt
        taumin = min(tau(l), taumin)
        taumax = max(tau(l), taumax)
      end do
      if (taumin*taumax.lt.czero) then
        ! if product of min and max is less than zero,
        ! taumin and taumax have different signs
        irst=2
      end if

      return
end
