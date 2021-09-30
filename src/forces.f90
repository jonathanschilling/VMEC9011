subroutine forces(guu, guv, gvv)

      use stel_kinds, only: dp
      use extfld, only: ivac
      use magfield, only: rbsq

      implicit none

      real(kind=dp) :: guu(nrztd)
      real(kind=dp) :: guv(nrztd)
      real(kind=dp) :: gvv(nrztd)

      real(kind=dp), pointer :: guus(nrztd)
      real(kind=dp), pointer :: guvs(nrztd)

      real(kind=dp) :: bsqr(nrztd)
      real(kind=dp) :: gvvs(nrztd)

      guus => crmn(1+nrztd)
      guvs => czmn(1+nrztd

      ! ON ENTRY, ARMN=ZU, BRMN=ZS, AZMN=RU, BZMN=RS, CZMN=R*BSQ.
      ! IN LOOPS, L (L1) INDEX REPRESENTS EVEN (ODD) COMPONENT.

      ! IT IS ESSENTIAL THAT CRMN, CZMN AT j=1 ARE ZERO INITIAL.
      ! crmn, czmn contain lv, lu on entry --> need to clear these
      do l = 1, nrzt+1, ns
        crmn(l) = czero
        czmn(l) = czero
      end do

      do l = 1, nrzt+1
        guus(l)  = guu(l)*shalf(l)
        guvs(l)  = guv(l)*shalf(l)
        gvvs(l)  = gvv(l)*shalf(l)
        armn(l)  = ohs*armn(l)*czmn(l)
        azmn(l)  =-ohs*azmn(l)*czmn(l)
        brmn(l)  = brmn(l)*czmn(l)
        bzmn(l)  =-bzmn(l)*czmn(l)
        bsqr(l)  = czmn(l)/shalf(l)
      end do

      do l = 2, nrzt+1
        l1 = l+nrzt
        armn(l1) = armn(l)*shalf(l)
        azmn(l1) = azmn(l)*shalf(l)
        brmn(l1) = brmn(l)*shalf(l)
        bzmn(l1) = bzmn(l)*shalf(l)
      end do

      ! CONSTRUCT CYLINDRICAL FORCE KERNELS (L=EVEN, L1=ODD COMPONENTS)
      do l = 1, nrzt
        l1 = l+nrzt
        bsqr(l) = cp25*(bsqr(l) + bsqr(l+1))
        czmn(l) = cp25*(czmn(l) + czmn(l+1))
        guu(l)  =  cp5*(guu(l)  +  guu(l+1))
        guus(l) =  cp5*(guus(l) + guus(l+1))
        guv(l)  =  cp5*(guv(l)  +  guv(l+1))
        guvs(l) =  cp5*(guvs(l) + guvs(l+1))
        gvv(l)  =  cp5*(gvv(l)  +  gvv(l+1))
        gvvs(l) =  cp5*(gvvs(l) + gvvs(l+1))
        armn(l) =  armn(l+1) - armn(l)
                 + cp5*(crmn(l) + crmn(l+1))
                 -     (  gvv (l)*r1(l )
                        + gvvs(l)*r1(l1) )
        crmn(l) =  cp5*(crmn(l)*shalf(l) + crmn(l+1)*shalf(l+1))
        azmn(l) = azmn(l+1) - azmn(l)
        brmn(l) =    cp5*(brmn(l) + brmn(l+1))
                  +  z1(l1)*bsqr(l)
                  - (  guu (l) * ru(l )
                     + guus(l) * ru(l1)
                     + guv (l) * rv(l )
                     + guvs(l) * rv(l1) )
        bzmn(l) = cp5*(bzmn(l) + bzmn(l+1))
                  - r1(l1)*bsqr(l)
                  - (  guu (l)*zu(l )
                     + guus(l)*zu(l1)
                     + guv (l)*zv(l )
                     + guvs(l)*zv(l1) )
      end do

      do l=1,nrzt
        l1 = l+nrzt

        s2 = sqrts(l)**2

        guus2 = guu(l)*s2
        guvs2 = guv(l)*s2
        gvvs2 = gvv(l)*s2

        armn(l1) = armn(l1+1) - armn(l1)
                   - (zu  (l1)*czmn(l) + zu(l)*bsqr(l )) + crmn(l)
                   - (gvvs(l )*r1  (l) + gvvs2*r1  (l1))
        azmn(l1) = azmn(l1+1) - azmn(l1)
                   + ru(l1)*czmn(l)
                   + ru(l )*bsqr(l)

        brmn(l1) = cp5*(brmn(l1) + brmn(l1+1))
                   + z1(l1)*czmn(l)
                   - (  guus(l)*ru(l )
                      + guus2  *ru(l1)
                      + guvs(l)*rv(l )
                      + guvs2  *rv(l1) )
        bzmn(l1) = cp5*(bzmn(l1) + bzmn(l1+1))
                   - r1(l1)*czmn(l)
                   - (  guus(l)*zu(l )
                      + guus2  *zu(l1)
                      + guvs(l)*zv(l )
                      + guvs2  *zv(l1) )

        crmn(l)  =   guv (l) * ru(l )
                   + guvs(l) * ru(l1)
                   + gvv (l) * rv(l )
                   + gvvs(l) * rv(l1)
        crmn(l1) =   guvs(l) * ru(l )
                   + guvs2   * ru(l1)
                   + gvvs(l) * rv(l )
                   + gvvs2   * rv(l1)

        czmn(l)  =   guv (l) * zu(l )
                   + guvs(l) * zu(l1)
                   + gvv (l) * zv(l )
                   + gvvs(l) * zv(l1)
        czmn(l1) =   guvs(l) * zu(l )
                   + guvs2   * zu(l1)
                   + gvvs(l) * zv(l )
                   + gvvs2   * zv(l1)
      end do

      ! ASSIGN EDGE FORCES (JS = NS) FOR FREE BOUNDARY CALCULATION
      if( ivac.ge.1 )then
        do m = 0,1
          do lk = 1,nznt
            l = ns*lk
            l1 = l + m*nrzt

            armn(l1) = armn(l1) + zu0(l)*rbsq(lk)
            azmn(l1) = azmn(l1) - ru0(l)*rbsq(lk)
          end do
        end do
      endif

      ! COMPUTE CONSTRAINT FORCE KERNELS
      do l = 1,nrzt
        l1 = l+nrzt
        rcon1   = (rcon(l) - rcon0(l)) * gcon(l)
        zcon1   = (zcon(l) - zcon0(l)) * gcon(l)

        brmn(l) = brmn(l) + rcon1
        bzmn(l) = bzmn(l) + zcon1
        brmn(l1)= brmn(l1)+ rcon1*sqrts(l)
        bzmn(l1)= bzmn(l1)+ zcon1*sqrts(l)

        rcon(l) =  ru0(l) * gcon(l)
        zcon(l) =  zu0(l) * gcon(l)
        rcon(l1)= rcon(l) * sqrts(l)
        zcon(l1)= zcon(l) * sqrts(l)
      end do

      return
end
