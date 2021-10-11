subroutine bcovar(bsubu, bsubv, gsqrt, bsq, r12, rs, zs, &
                  ru12, zu12, guu, guv, gvv, phipog, lu, lv)

      use stel_kinds, only: dp
      use name0, only: czero, cp5, c1p0, c2p0
      use name1, only: nsd
      use scalars, only: nrzt, ns, hs, dnorm, iequi, meven, modd, &
                         iter1, iter2, ns4, mns
      use scalefac, only: sqrts, shalf, wint
      use realsp, only: r1, ru, rv, z1, zu, zv, ru0, zu0
      use profs, only: phip, vp
      use fsqu, only: wp, wb, fnorm, fnorm1
      use precond, only: arm, ard, brm, brd,    &
                         azm, azd, bzm, bzd, cr
      use xstuff, only: xc
      use spectra, only: tcon

      implicit none

      ! TODO: more elegant definition of these BLAS functions
      real(kind=dp) :: ddot

      real(kind=dp), intent(out)   :: bsubu (nrzt,0:1) ! --> clmn (~d(lambda)/d( zeta) => force on  zeta-derivative of lambda)
      real(kind=dp), intent(out)   :: bsubv (nrzt,0:1) ! --> blmn (~d(lambda)/d(theta) => force on theta-derivative of lambda)
      real(kind=dp)                :: gsqrt (nrzt)
      real(kind=dp)                :: bsq   (nrzt)
      real(kind=dp)                :: r12   (nrzt)
      real(kind=dp)                :: rs    (nrzt)
      real(kind=dp)                :: zs    (nrzt)
      real(kind=dp)                :: ru12  (nrzt)
      real(kind=dp)                :: zu12  (nrzt)
      real(kind=dp), intent(out)   :: guu   (nrzt)
      real(kind=dp), intent(out)   :: guv   (nrzt)
      real(kind=dp), intent(out)   :: gvv   (nrzt)
      real(kind=dp)                :: phipog(nrzt)
      real(kind=dp), intent(inout) :: lu    (nrzt,0:1)
      real(kind=dp), intent(inout) :: lv    (nrzt,0:1)

      real(kind=dp) :: ar(nsd)
      real(kind=dp) :: az(nsd)

      integer       :: is, js, l, lk, lme, lmo, m
      real(kind=dp) :: volume

      ! INITIALIZATION BLOCK
      do l = 1, nrzt+1 ! NOTE: incl. magic element at end !
         guu(l) = czero
         guv(l) = czero
         gvv(l) = czero
      end do

      ! COMPUTE METRIC ELEMENTS GIJ ON HALF MESH
      ! loop over is leads to summation in radial direction over l-1 and l
      do is = -1, 0
        do l = 2, nrzt
          lme = l + is     ! index for l-th even-m
          lmo = lme + nrzt ! index for l-th  odd-m

          ! temporary storage of (0.5 * s)
          phipog(l) = cp5 * sqrts(lme)*sqrts(lme) ! depending on is, this is s_j or s_(j-1)

          guu(l) = guu(l) + cp5            *( ru(lme)*ru(lme) + zu(lme)*zu(lme)                                                       ) &
                          +       phipog(l)*( ru(lmo)*ru(lmo) + zu(lmo)*zu(lmo)                                                       ) &
                          +        shalf(l)*( ru(lme)*ru(lmo) + zu(lme)*zu(lmo)                                                       )

          guv(l) = guv(l) + cp5            *( ru(lme)*rv(lme) + zu(lme)*zv(lme)                                                       ) &
                          +       phipog(l)*( ru(lmo)*rv(lmo) + zu(lmo)*zv(lmo)                                                       ) &
                          + cp5 *  shalf(l)*( ru(lme)*rv(lmo) + rv(lme)*ru(lmo) + zu(lme)*zv(lmo) + zv(lme)*zu(lmo)                   )

          gvv(l) = gvv(l) + cp5            *( rv(lme)*rv(lme) + zv(lme)*zv(lme)                                     + r1(lme)*r1(lme) ) &
                          +       phipog(l)*( rv(lmo)*rv(lmo) + zv(lmo)*zv(lmo)                                     + r1(lmo)*r1(lmo) ) &
                          +        shalf(l)*( rv(lme)*rv(lmo) + zv(lme)*zv(lmo)                                     + r1(lme)*r1(lmo) )
        end do
      end do

      ! PUT LAMBDA DERIVATIVES ON RADIAL HALF-MESH
      ! This also re-scales lu, lv to include a factor phip/gsqrt == phipog
      ! for the zero-current algorithm that needs to be kept in mind for later on!
      do l = nrzt,2,-1
        phipog(l) = phip(l)/gsqrt(l)

        lu(l,0) = phipog(l) * cp5*( lu(l,0)+lu(l-1,0) + shalf(l) * (lu(l,1)+lu(l-1,1)) )
        lv(l,0) = phipog(l) * cp5*( lv(l,0)+lv(l-1,0) + shalf(l) * (lv(l,1)+lv(l-1,1)) )
      end do

      ! COMPUTE IOTA PROFILE (consistent with prescribed current profile if ncurr != 0)
      call getiota(phipog, guu, guv, wint, lu, lv)

      ! PUT LAMBDA FORCES ON RADIAL HALF-MESH.
      ! (The lambda forces are the covariant magnetic field components.)
      do l = 1, nrzt
        bsubu(l,0) = guu(l)*lv(l,0) + guv(l)*lu(l,0)
        bsubv(l,0) = guv(l)*lv(l,0) + gvv(l)*lu(l,0)
        bsubu(l,1) = shalf(l)*bsubu(l,0)
        bsubv(l,1) = shalf(l)*bsubv(l,0)
      end do

      ! COMPUTE SUM OF KINETIC AND MAGNETIC PRESSURES
      ! AND RETURN IF FINAL OUTPUT LOOP (IEQUI=1)
      ! ON ENTRY, BSQ IS THE KINETIC PRESSURE ( =pres on each surface; see pressure() )
      wb = -wp
      do l = 2,nrzt
        bsq(l)    = bsq(l) + cp5 * (lv(l,0)*bsubu(l,0) + lu(l,0)*bsubv(l,0))
        wb        = wb + hs*dnorm*wint(l)*abs(gsqrt(l))*bsq(l)

        phipog(l) = phipog(l)*wint(l)
      end do

      if (iequi.eq.1) &
        return

      ! AVERAGE LAMBDA FORCES ONTO FULL MESH
      ! NOTE: EDGE FORCE IS DOWN BY .5 UNTIL 90 LOOP
      do m = meven, modd
        do l = 1, nrzt
          bsubu(l,m) = cp5*(bsubu(l,m) + bsubu(l+1,m))
          bsubv(l,m) = cp5*(bsubv(l,m) + bsubv(l+1,m))
        end do

        ! 90 loop: scale up edge components by 2.0
        do l = ns, nrzt, ns
          bsubu(l,m) = c2p0*bsubu(l,m)
          bsubv(l,m) = c2p0*bsubv(l,m)
        end do
      end do

      ! COMPUTE R,Z AND LAMBDA PRE-CONDITIONING MATRIX ELEMENTS
      ! AND FORCE NORMS EVERY NS4(=25) STEPS
      if (mod(iter2-iter1, ns4) .eq. 0) then

        call lamcal(phipog, guu, guv, gvv)

        call precondn(lu, bsq, gsqrt, r12,                        &
                      zs, zu12, zu, zu(1+nrzt), z1(1+nrzt), wint, &
                      arm, ard, brm, brd, cr)

        call precondn(lu, bsq, gsqrt, r12,                        &
                      rs, ru12, ru, ru(1+nrzt), r1(1+nrzt), wint, &
                      azm, azd, bzm, bzd, cr)

        do l = 2, nrzt
          guu(l) = guu(l) * r12(l)**2
        end do

        volume = hs * sum(vp(2:ns))
        fnorm  = dnorm/(ddot(nrzt, guu, 1, wint, 1) * (wb/volume)**2)

        ! leave out forces on axis (?) --> start at ns+1, 4*mns-ns elements
        fnorm1 = c1p0/ddot(4*mns-ns, xc(ns+1), 1, xc(ns+1), 1)

        ! COMPUTE CONSTRAINT FORCE SCALING FACTOR (TCON)
        ! ar, az: flux-surface averages of (dR/du)^2, (dZ/du)^2 (without Jacobian...)
        do js = 2, ns-1
          ar(js) = czero
          az(js) = czero
        end do
        do js = 2, ns-1
          do lk = 1, nrzt, ns
            ar(js) = ar(js) + wint(js+lk-1) * ru0(js+lk-1)**2
            az(js) = az(js) + wint(js+lk-1) * zu0(js+lk-1)**2
          end do
        end do

        ! TODO: plot tcon profile over iterations
        do js = 2, ns-1
          tcon(js) = min(abs(ard(js,1)/ar(js)), abs(azd(js,1)/az(js)))
        end do
        tcon(ns) = cp5*tcon(ns-1)

      endif ! update preconditioner

      ! STORE LU * LV COMBINATIONS USED IN FORCES
      do l = 2, nrzt
         guu(l) = lv(l,0)*lv(l,0)*gsqrt(l)
         guv(l) = lv(l,0)*lu(l,0)*gsqrt(l)
         gvv(l) = lu(l,0)*lu(l,0)*gsqrt(l)

         lv(l,0)  = bsq(l)*gsqrt(l)/r12(l)
         lu(l,0)  = bsq(l)*r12(l)
      end do

      return
end
