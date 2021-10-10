subroutine forces(guu, guv, gvv)

      use stel_kinds, only: dp
      use name0, only: czero, cp25, cp5
      use name1, only: nznt, nrztd
      use extfld, only: ivac
      use magfield, only: rbsq
      use realsp, only: r1, ru, rv, z1, zu, zv, ru0, zu0, &
                        rcon, zcon, rcon0, zcon0, gcon
      use rforces, only: worka, &
                         armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn
      use scalars, only: ns, ohs, nrzt
      use scalefac, only: sqrts, shalf

      implicit none

      ! on entry, these contain (lu*lv combinations at end of bcovar):
      ! guu : sqrt(g)*B^u*B^u
      ! guv : sqrt(g)*B^u*B^v
      ! gvv : sqrt(g)*B^v*B^v
      ! on the half-grid.
      ! on exit, they contain: ???
      real(kind=dp), intent(inout) :: guu(nrztd)
      real(kind=dp), intent(inout) :: guv(nrztd)
      real(kind=dp), intent(inout) :: gvv(nrztd)

      real(kind=dp), pointer :: guus(:)
      real(kind=dp), pointer :: guvs(:)
      real(kind=dp)          :: gvvs(nrztd) ! gvvs is stored here ?
      real(kind=dp)          :: bsqr(nrztd)

      ! IN LOOPS, L (L1) INDEX REPRESENTS EVEN (ODD) COMPONENT.
      integer       :: l, l1, lk, m
      real(kind=dp) :: rcon1, zcon1
      real(kind=dp) :: s2, guus2, guvs2, gvvs2

      guus => worka(1+ 5*nrztd: 6*nrztd) ! crmn(lodd)
      guvs => worka(1+11*nrztd:12*nrztd) ! czmn(lodd)

      ! ON ENTRY, ARMN=ZU, BRMN=ZS, AZMN=RU, BZMN=RS, CZMN=R*BSQ.
      ! ^^ on half-grid: zu12, ru12 !!! from jacobian()

      ! IT IS ESSENTIAL THAT CRMN, CZMN AT j=1 ARE ZERO INITIAL.
      ! crmn, czmn contain lv, lu on entry --> need to clear these
      do l = 1, nrzt+1, ns
        crmn(l) = czero
        czmn(l) = czero
      end do

      do l = 1, nrzt+1
        ! 1
        guus(l)  = guu(l)*shalf(l)
        guvs(l)  = guv(l)*shalf(l)
        gvvs(l)  = gvv(l)*shalf(l)

        ! 2
        armn(l)  = ohs*armn(l)*czmn(l)
        azmn(l)  =-ohs*azmn(l)*czmn(l)

        ! 3
        brmn(l)  = brmn(l)*czmn(l)
        bzmn(l)  =-bzmn(l)*czmn(l)

        ! 4
        bsqr(l)  = czmn(l)/shalf(l)
      end do

      ! 5
      do l = 2, nrzt+1
        l1 = l+nrzt
        armn(l1) = armn(l)*shalf(l)
        azmn(l1) = azmn(l)*shalf(l)
        brmn(l1) = brmn(l)*shalf(l)
        bzmn(l1) = bzmn(l)*shalf(l)
      end do

      ! CONSTRUCT CYLINDRICAL (R,Z) FORCE KERNELS (L=EVEN, L1=ODD COMPONENTS)
      ! even-m force contributions
      do l = 1, nrzt
        l1 = l + nrzt

        ! average half-grid quantities onto full grid
        ! 6
        bsqr(l) = cp25*(bsqr(l) + bsqr(l+1)) ! additional factor of 1/2 ???
        czmn(l) = cp25*(czmn(l) + czmn(l+1)) ! additional factor of 1/2 ???
        guu (l) =  cp5*(guu (l) + guu (l+1))
        guus(l) =  cp5*(guus(l) + guus(l+1))
        guv (l) =  cp5*(guv (l) + guv (l+1))
        guvs(l) =  cp5*(guvs(l) + guvs(l+1))
        gvv (l) =  cp5*(gvv (l) + gvv (l+1))
        gvvs(l) =  cp5*(gvvs(l) + gvvs(l+1))

        ! actually do some MHD
        ! These are the "final" forms of even-m aXmn and bXmn where X=r,z
        armn(l) = armn(l+1) - armn(l) + cp5*(crmn(l) + crmn(l+1)) - (gvv(l)*r1(l) + gvvs(l)*r1(l1)) ! 7
        azmn(l) = azmn(l+1) - azmn(l)                                                               ! 8

        brmn(l) =   z1(l1)*bsqr(l)    + cp5*(brmn(l) + brmn(l+1)) - (guu(l)*ru(l) + guus(l)*ru(l1) + guv(l)*rv(l) + guvs(l)*rv(l1)) ! 9
        bzmn(l) =  -r1(l1)*bsqr(l)    + cp5*(bzmn(l) + bzmn(l+1)) - (guu(l)*zu(l) + guus(l)*zu(l1) + guv(l)*zv(l) + guvs(l)*zv(l1)) ! 10

        ! crmn contains some intermediate quantity required below
        crmn(l) = cp5*(crmn(l)*shalf(l) + crmn(l+1)*shalf(l+1)) ! 11
      end do

      ! odd-m force contributions
      do l = 1, nrzt
        l1 = l + nrzt

        ! s = sqrt(s)*sqrt(s) on the full grid
        s2 = sqrts(l)**2

        ! 12
        guus2 = guu(l)*s2
        guvs2 = guv(l)*s2
        gvvs2 = gvv(l)*s2

        ! final "MHD" form of odd-m parts of aXmn, bXmn
        ! 13
        armn(l1) = armn(l1+1) - armn(l1) - zu(l1)*czmn(l) - zu(l)*bsqr(l) + crmn(l) - (gvvs(l)*r1(l) + gvvs2*r1(l1))
        azmn(l1) = azmn(l1+1) - azmn(l1) + ru(l1)*czmn(l) + ru(l)*bsqr(l)

        ! 14
        brmn(l1) =     z1(l1)*czmn(l)    + cp5*(brmn(l1+1) + brmn(l1)) - (guus(l)*ru(l) + guus2*ru(l1) + guvs(l)*rv(l) + guvs2*rv(l1))
        bzmn(l1) =    -r1(l1)*czmn(l)    + cp5*(bzmn(l1+1) + bzmn(l1)) - (guus(l)*zu(l) + guus2*zu(l1) + guvs(l)*zv(l) + guvs2*zv(l1))

        ! actually build final form of cXmn now that temporary stuff in crmn, czmn is no longer needed
        ! 15
        crmn(l)  = guv (l)*ru(l) + guvs(l)*ru(l1) + gvv (l)*rv(l) + gvvs(l)*rv(l1)
        crmn(l1) = guvs(l)*ru(l) + guvs2  *ru(l1) + gvvs(l)*rv(l) + gvvs2  *rv(l1)

        ! 16
        czmn(l)  = guv (l)*zu(l) + guvs(l)*zu(l1) + gvv (l)*zv(l) + gvvs(l)*zv(l1)
        czmn(l1) = guvs(l)*zu(l) + guvs2  *zu(l1) + gvvs(l)*zv(l) + gvvs2  *zv(l1)
      end do

      ! ASSIGN EDGE FORCES (JS = NS) FOR FREE BOUNDARY CALCULATION
      if (ivac.ge.1) then
        do m = 0, 1 ! parity of m
          do lk = 1,nznt
            l = ns*lk
            l1 = l + m*nrzt

            armn(l1) = armn(l1) + zu0(l)*rbsq(lk)
            azmn(l1) = azmn(l1) - ru0(l)*rbsq(lk)
          end do
        end do
      endif

      ! COMPUTE CONSTRAINT FORCE KERNELS
      do l = 1, nrzt
        l1 = l+nrzt
        rcon1    = (rcon(l) - rcon0(l)) * gcon(l)
        zcon1    = (zcon(l) - zcon0(l)) * gcon(l)

        ! take constraint force into account also in bXmn
        brmn(l)  = brmn(l)  + rcon1
        bzmn(l)  = bzmn(l)  + zcon1
        brmn(l1) = brmn(l1) + rcon1*sqrts(l)
        bzmn(l1) = bzmn(l1) + zcon1*sqrts(l)

        ! these become arcon, azcon in tomnsp()
        rcon(l)  =  ru0(l) * gcon(l)
        zcon(l)  =  zu0(l) * gcon(l)
        rcon(l1) = rcon(l) * sqrts(l)
        zcon(l1) = zcon(l) * sqrts(l)
      end do

      return
end
