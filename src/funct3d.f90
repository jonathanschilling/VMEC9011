subroutine funct3d

      use stel_kinds, only: dp, p4
      use name0, only: czero, cp25, cp5, c1p0, c1p5
      use name1, only: nvac, nznt, nrztd, mnmax
      use realsp, only: workb, &
                        r1, ru, rv, z1, zu, zv, rcon, zcon, &
                        ru0, zu0, rcon0, zcon0, gcon
      use rforces, only: worka, &
                         armn, brmn, crmn, azmn, bzmn, czmn, &
                         blmn, clmn
      use scalars, only: nrzt, mns, ns, neqs, iter1, iter2, &
                         hs, ohs, irst,  meven, modd, iequi, dnorm, &
                         twopi, voli, isigng
      use xstuff, only: xc, gc
      use mnarray, only: xrz3, xrz4
      use scalefac, only: scalxc, sqrts, shalf, wint
      use fsqu, only: wp, fsqr, fsqz
      use profs, only: mass, vp, pres
      use inputdat, only: gam, nvacskip
      use extfld, only: ivac, ivac2
      use magfield, only: curpol, rbtor, ctor, bsqvac, bsqsav, rbsq, dbsq
      use time, only: timer

      implicit none

      ! TODO: more elegant definition of these BLAS functions
      real(kind=dp) :: ddot
      external dcopy

      real(kind=dp) :: rmnc(mnmax)
      real(kind=dp) :: zmns(mnmax)
      real(kind=dp) :: lmns(mnmax)
      real(kind=dp) :: xm(mnmax) ! for handing over to vacuum()
      real(kind=dp) :: xn(mnmax) ! for handing over to vacuum()

      real(kind=dp) :: rax(nznt)
      real(kind=dp) :: zax(nznt)

      integer       :: lodd, l, lk
      real(kind=dp) :: bz0
      real(kind=p4) :: t4
      real(kind=dp) :: timeon, timeoff

      ! let the madness begin ...
      real(kind=dp), pointer :: guu(:)
      real(kind=dp), pointer :: guv(:)
      real(kind=dp), pointer :: gvv(:)
      real(kind=dp), pointer :: lu(:)
      real(kind=dp), pointer :: lv(:)

      guu => workb(1+13*nrztd:14*nrztd) ! rcon(1+nrztd)
      guv => workb(1+15*nrztd:16*nrztd) ! zcon(1+nrztd)
      gvv => workb(1+ 6*nrztd: 7*nrztd) ! z1
      lu  => worka(1+10*nrztd:12*nrztd) ! czmn
      lv  => worka(1+ 4*nrztd: 6*nrztd) ! crmn

      lodd = 1+nrzt

      ! EXTRAPOLATE M>2 MODES AT JS = 2
      call extrap(xc,          xc(1+  mns), &
                  xc(1+2*mns), xc(1+3*mns), &
                  xc(1+4*mns), xc(1+5*mns), &
                  xrz3, xrz4)

      ! temporary re-use of gc for scaled xc
      do l = 1,neqs
        gc(l) = xc(l) * scalxc(l)
      enddo

      ! INVERSE FOURIER TRANSFORM TO S,THETA,ZETA SPACE
      ! R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
      do l = 1, nrzt
        lu(l) = c1p0 !!! This makes lu == (1 + d(lambda)/du)
        lv(l) = czero
        lu(l+nrzt) = czero
        lv(l+nrzt) = czero
      enddo
      ! note that lv coming out of totzsp is actually -d(lambda)/dv
      call totzsp(gc, gc(1+mns), gc(1+2*mns), gc(1+3*mns), gc(1+4*mns), gc(1+5*mns), &
                  r1, ru,        rv,          z1,          zu,          zv,          &
                  lu, lv,        rcon,        zcon, &
                  worka,worka,worka)

      ! COMPUTE CONSTRAINT FORCE (GCON)
      do l = 1,nrzt
        ru0(l)  = ru  (l) + ru  (l+nrzt)*sqrts(l)
        zu0(l)  = zu  (l) + zu  (l+nrzt)*sqrts(l)
        rcon(l) = rcon(l) + rcon(l+nrzt)*sqrts(l)
        zcon(l) = zcon(l) + zcon(l+nrzt)*sqrts(l)

        ! init gcon to zero
        gcon(l) = czero
        ! azmn: temporary storage for not-yet-dealiased gcon
        ! rcon0, zcon0 are initialized to 0 in vsetup --> ok to use in first iteration
        azmn(l) = (rcon(l)-rcon0(l))*ru0(l) + (zcon(l)-zcon0(l))*zu0(l)
      enddo

      do l = 1,2*mns
        gc(l) = czero
      enddo

      ! first iteration: save rcon into rcon0, zcon into zcon0
      if (iter2.eq.1 .and. nvac.eq.0) then
        call dcopy(nrzt, rcon, 1, rcon0, 1)
        call dcopy(nrzt, zcon, 1, zcon0, 1)
      endif

      ! only put non-zero values into gcon if iter2>1
      ! --> in the first iteration, gcon stays zero
      ! --> no constraint force in first iteration
      if (iter2.gt.1) &
        call alias(gcon, azmn, worka,worka,worka, gc, gc(1+mns))

      ! COMPUTE S AND THETA DERIVATIVE OF R AND Z AND JACOBIAN ON HALF-GRID
      call jacobian(r1, ru, z1, zu,                            &
      !             zu12        ru12        zs          rs
                    armn,       azmn,       brmn,       bzmn,  &
                    azmn(lodd), armn(lodd), brmn(lodd)        )
      !             gsqrt       r12         tau
      if (irst.eq.2 .and. iequi.eq.0) then
        ! Jacobian changes sign somewhere in volume
        ! --> flux surfaces cross, need to restart with smaller timestep
        return
      end if

      ! COMPUTE PRESSURE AND VOLUME ON HALF-GRID
      ! see above for re-use of arrays; additionally:
      !             gsqrt       bsq
      call pressure(azmn(lodd), bzmn(lodd), wint)
      if (iter2.eq.1) then
        voli = (twopi**2)*hs*sum(vp(2:ns))
      end if

      ! COMPUTE COVARIANT COMPONENTS OF B, MAGNETIC PRESSURE,
      ! AND METRIC ELEMENTS ON HALF-GRID
      do lk = 1,nznt
        rax(lk) = r1(1+ns*(lk-1))
        zax(lk) = z1(1+ns*(lk-1))
      enddo

      ! backup z1 in gco (since z1 gets overwritten in bcovar for iequi=1?)
      if (iequi.eq.1) &
        call dcopy(nrzt, z1, 1, gcon, 1)

      ! see above for re-use of arrays; additionally:
      !           bsubu  bsubv  gsqrt       bsq         r12         rs    zs
      call bcovar(clmn,  blmn,  azmn(lodd), bzmn(lodd), armn(lodd), bzmn, brmn, &
                  azmn,  armn,  guu,        guv,        gvv,        brmn(lodd), lu, lv)
      !           ru12   zu12                                       phipog

      ! restore z1 from gcon (since z1 got overwritten in bcovar for iequi=1 ?)
      if (iequi.eq.1) &
        call dcopy(nrzt, gcon, 1, z1, 1)

      ! bsubv@LCFS --> <B_zeta> is being computed here
      bz0 = ddot(nznt, blmn(ns), ns, wint(ns), ns)

      ! COMPUTE VACUUM MAGNETIC PRESSURE AT PLASMA EDGE
      ! NOTE: FOR FREE BOUNDARY RUNS, THE PLASMA VOLUME CAN
      ! BE INCREASED BY DECREASING CURPOL AT FIXED PHIPS.  THE
      ! VALUE FOR CURPOL GIVEN BELOW SHOULD BE USED AS AN
      ! INITIAL GUESS, WHICH CAN BE CHANGED BY READING IN
      ! CURPOL FROM DATA STATEMENT
      if (nvac.ne.0.and.iter2.gt.1) then
        call second(t4)
        timeon = real(t4)

        ivac2=mod(iter2-iter1, nvacskip)

        if ((fsqr+fsqz).le.1.e-2) &
          ivac=ivac+1

        if (ivac.eq.1) &
          ivac2 = 0

        rbtor = curpol

        if (ivac.eq.1 .and. rbtor.eq.0.) &
          rbtor = twopi*dnorm*bz0

        if (ivac.ge.1) then
          ctor = isigng*twopi*dnorm*(                          &
                   c1p5*ddot(nznt, clmn(ns  ), ns, wint(ns), ns)   &
                 - cp5 *ddot(nznt, clmn(ns-1), ns, wint(ns), ns) )

          call convert(rmnc, zmns, lmns, xm, xn, ns,          &
                       xc,          xc(1+  mns), xc(1+2*mns), &
                       xc(1+3*mns), xc(1+4*mns), xc(1+5*mns) )

          call vacuum(rmnc, zmns, xm, xn, ctor, rbtor, bsqvac, rax, zax)
        endif

        do lk = 1, nznt
          !                   bsq
          bsqsav(lk,3) = c1p5*bzmn(ns*lk+nrzt) - cp5*bzmn(ns*lk-1+nrzt)
          rbsq(lk)     = bsqvac(lk) * ohs*(r1(ns*lk) + r1(ns*lk+nrzt))
          dbsq(lk)     = abs(bsqvac(lk)-bsqsav(lk,3))
        enddo

        if(ivac.eq.1) then
          ! copy initial plasma and vacuum field into bsqsav(:,1:2)
          call dcopy(nznt, bzmn(ns+nrzt:ns+2*nrzt), ns, bsqsav(1,1), 1)
          call dcopy(nznt, bsqvac,                   1, bsqsav(1,2), 1)
        endif

        call second(t4)
        timeoff = real(t4)
        timer(1) = timer(1) + (timeoff-timeon)

      endif ! vacuum/free-boundary contribution

      if(iequi.eq.1)then
        ! COMPUTE REMAINING COVARIANT COMPONENT OF B (BSUBS),
        ! CYLINDRICAL COMPONENTS OF B (BR, BPHI, BZ), AND
        ! AVERAGE EQUILIBRIUM PROPERTIES AT END OF RUN

        !        r12         rs    zs    ru12  zu12  bsubs
        call bss(armn(lodd), bzmn, brmn, azmn, armn, crmn(lodd), &
                 lu, lv, rcon, czmn(lodd), zcon)
        !                br    bphi        bz

        ! [bcovar] bsubu bsubv bsq
        ! [eqfor]  bu    bv    bsq         rmag  zmag
        call eqfor(clmn, blmn, bzmn(lodd), xc,   xc(1+2*mns))

        !          bsq         gsqrt       bsubu bsubv
        call wrout(bzmn(lodd), azmn(lodd), clmn, blmn,        &
                   crmn(lodd), rcon, czmn(lodd), zcon, lu, lv)
        !          bsubs       br    bphi        bz

        return
      endif ! output file quantities

      ! COMPUTE MHD FORCES ON INTEGER-MESH
      call forces(guu, guv, gvv)

      ! FOURIER-TRANSFORM MHD FORCES TO (M,N)-SPACE
      do l = 1, neqs
        gc(l) = czero
      enddo
      call tomnsp(gc, gc(1+mns), gc(1+2*mns), gc(1+3*mns), gc(1+4*mns), gc(1+5*mns), &
                  armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn, rcon, zcon, workb,workb,workb)
      do l = 1, neqs
        ! TODO: Can one maybe omit this if the odd-m force components are not scaled by sqrt(s) ?
        ! --> comment in newer VMEC versions; are the odd-m force components different there?
        ! The hnr22_kis test case converges in (3 ... 5) % less iterations if this is commented out!
        gc(l) = gc(l) * scalxc(l)
      enddo

      ! COMPUTE FORCE RESIDUALS
      call residue(gc, gc(1+2*mns), gc(1+4*mns), bz0, workb)

      return

end
