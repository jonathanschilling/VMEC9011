subroutine bcovar(bsubu, bsubv, gsqrt, bsq, r12, rs, zs, &
                  ru12, zu12, guu, guv, gvv, phipog, lu, lv)

      use stel_kinds, only: dp

      implicit none

      common/precond/ard(nsd1,2),arm(nsd1,2),brd(nsd1,2),brm(nsd1,2),cr(nsd1),
                     azd(nsd1,2),azm(nsd1,2),bzd(nsd1,2),bzm(nsd1,2)

      real(kind=dp), intent(out) :: bsubu(nrzt,0:1)
      real(kind=dp), intent(out) :: bsubv(nrzt,0:1)
      real(kind=dp)              :: gsqrt(*)
      real(kind=dp)              :: bsq(*)
      real(kind=dp)              :: r12(*)
      real(kind=dp)              :: rs(*)
      real(kind=dp)              :: zs(*)
      real(kind=dp)              :: ru12(*)
      real(kind=dp)              :: zu12(*)
      real(kind=dp), intent(out) :: guu(*)
      real(kind=dp), intent(out) :: guv(*)
      real(kind=dp), intent(out) :: gvv(*)
      real(kind=dp)              :: phipog(*)
      real(kind=dp)              :: lu(nrzt,0:1)
      real(kind=dp)              :: lv(nrzt,0:1)

      real(kind=dp) :: ar(nsd)
      real(kind=dp) :: az(nsd)

      ! INITIALIZATION BLOCK
      do l = 1, nrzt+1
         guu(l) = czero
         guv(l) = czero
         gvv(l) = czero
      end do

      ! COMPUTE METRIC ELEMENTS GIJ ON HALF MESH
      ! loop over is leads to radial averaging over l-1 and l
      do is = -1,0
        do l=2,nrzt
          lme = l + is
          lmo = lme + nrzt

          phipog(l) = cp5*sqrts(lme)*sqrts(lme)

          guu(l) = guu(l)
                    +       cp5*(  ru(lme)*ru(lme) + zu(lme)*zu(lme))
                    + phipog(l)*(  ru(lmo)*ru(lmo) + zu(lmo)*zu(lmo))
                    +  shalf(l)*(  ru(lme)*ru(lmo) + zu(lme)*zu(lmo))

          guv(l) = guv(l)
                    +       cp5*(  ru(lme)*rv(lme) + zu(lme)*zv(lme))
                    + phipog(l)*(  ru(lmo)*rv(lmo) + zu(lmo)*zv(lmo))
                    +  shalf(l)*(  ru(lme)*rv(lmo) + rv(lme)*ru(lmo)
                                 + zu(lme)*zv(lmo) + zv(lme)*zu(lmo)) * cp5

          gvv(l) = gvv(l)
                    +       cp5*(  rv(lme)*rv(lme) + zv(lme)*zv(lme))
                    + phipog(l)*(  rv(lmo)*rv(lmo) + zv(lmo)*zv(lmo))
                    +  shalf(l)*(  rv(lme)*rv(lmo) + zv(lme)*zv(lmo))
                    +       cp5*   r1(lme)*r1(lme)
                    + phipog(l)*   r1(lmo)*r1(lmo)
                    +  shalf(l)*   r1(lme)*r1(lmo)
        end do
      end do

      ! PUT LAMBDA DERIVATIVES ON RADIAL HALF-MESH
      do l = nrzt,2,-1
        phipog(l) = phip(l)/gsqrt(l)

        lu(l,0) = cp5*phipog(l)*(lu(l,0)+lu(l-1,0)
                +      shalf(l)*(lu(l,1)+lu(l-1,1)))

        lv(l,0) = cp5*phipog(l)*(lv(l,0)+lv(l-1,0)
                +      shalf(l)*(lv(l,1)+lv(l-1,1)))
      end do

      ! COMPUTE IOTA PROFILE
      call getiota(phipog, guu, guv, lu, lv, wint, iotas, jv, czero, ns, ncurr)

      ! PUT LAMBDA FORCES (=covariant magnetic field components)
      ! ON RADIAL HALF-MESH
      do l = 1,nrzt
        bsubu(l,0) = guu(l)*lv(l,0) + guv(l)*lu(l,0)
        bsubv(l,0) = guv(l)*lv(l,0) + gvv(l)*lu(l,0)
        bsubu(l,1) = shalf(l)*bsubu(l,0)
        bsubv(l,1) = shalf(l)*bsubv(l,0)
      end do

      ! COMPUTE SUM OF KINETIC AND MAGNETIC PRESSURES
      ! AND RETURN IF FINAL OUTPUT LOOP (IEQUI=1)
      ! ON ENTRY, BSQ IS THE KINETIC PRESSURE
      wb = -wp
      do l = 2,nrzt
        bsq(l) = bsq(l) + cp5*(lv(l,0)*bsubu(l,0) + lu(l,0)*bsubv(l,0))
        phipog(l) = phipog(l)*wint(l)
        wb = wb + hs*dnorm*wint(l)*abs(gsqrt(l))*bsq(l)
      end do

      if (iequi.eq.1) &
        return

      ! AVERAGE LAMBDA FORCES ONTO FULL MESH
      ! NOTE: EDGE FORCE IS DOWN BY .5 UNTIL 90 LOOP
      do m = meven,modd
        do l=1,nrzt
          bsubu(l,m) = cp5*(bsubu(l,m) + bsubu(l+1,m))
          bsubv(l,m) = cp5*(bsubv(l,m) + bsubv(l+1,m))
        end do

        ! do 70 l=1,nrzt-1
        !        bsubu(l,m) = cp5*(bsubu(l,m) + bsubu(l+1,m))
        !        bsubv(l,m) = cp5*(bsubv(l,m) + bsubv(l+1,m))
        ! 70     continue

        ! l     = nrzt
        ! bsubu(l,m) = cp5*(bsubu(l,m)               )
        ! bsubv(l,m) = cp5*(bsubv(l,m)               )
        do l = ns,nrzt,ns
          bsubu(l,m) = c2p0*bsubu(l,m)
          bsubv(l,m) = c2p0*bsubv(l,m)
        end do
      end do

      !  COMPUTE R,Z AND LAMBDA PRE-CONDITIONING MATRIX
      ! ELEMENTS AND FORCE NORMS EVERY NS4 STEPS
      if(mod(iter2-iter1,ns4).eq.0)then

        call lamcal(phipog, guu, guv, gvv)

        call precondn(lu, bsq, gsqrt, r12, zs, zu12, zu, zu(1+nrzt), z1(1+nrzt), &
                      shalf, wint, pres, &
                      arm, ard, brm, brd, cr, &
                      ohs, cp25, cp5, c1p0, c1p5, czero, ns, iter2)

        call precondn(lu, bsq, gsqrt, r12, rs, ru12, ru, ru(1+nrzt), r1(1+nrzt), &
                      shalf, wint, pres, &
                      azm, azd, bzm, bzd, cr, &
                      ohs, cp25, cp5, c1p0, c1p5, czero, ns, iter2)

        do l=2,nrzt
          guu(l) = guu(l)*r12(l)**2
        end do

        volume = hs*ssum(ns-1,vp(2),1)

        fnorm = dnorm/(sdot(nrzt,guu,1,wint,1)*(wb/volume)**2)

        fnorm1 = c1p0/sdot(4*mns-ns,xc(ns+1),1,xc(ns+1),1)

        ! COMPUTE CONSTRAINT FORCE SCALING FACTOR (TCON)
        do js=2,ns-1
          ar(js) = czero
          az(js) = czero
        end do

        do js=2,ns-1
          do lk = 1,nrzt,ns
            ar(js) = ar(js) + wint(js+lk-1)*ru0(js+lk-1)**2
            az(js) = az(js) + wint(js+lk-1)*zu0(js+lk-1)**2
          end do
        end do

        do js=2,ns-1
          tcon(js) = min(abs(ard(js,1)/ar(js)),abs(azd(js,1)/az(js)))
        end do
        tcon(ns) = cp5*tcon(ns-1)
      endif

      ! STORE LU * LV COMBINATIONS USED IN FORCES
      do l=2,nrzt
         guu(l) = lv(l,0)*lv(l,0)*gsqrt(l)
         guv(l) = lv(l,0)*lu(l,0)*gsqrt(l)
         gvv(l) = lu(l,0)*lu(l,0)*gsqrt(l)

         lv(l,0)  = bsq(l)*gsqrt(l)/r12(l)
         lu(l,0)  = bsq(l)*r12(l)
      end do

      return
end
