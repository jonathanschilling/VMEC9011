subroutine funct3d

      use stel_kinds, only: dp
      use name0, only: czero, cp25, cp5, c1p0, c1p5
      use name1, only: nvac
      use realsp, only: r1, ru, rv, z1, zu, zv, &
                        ru0, zu0, rcon, zcon, rcon0, zcon0, gcon
      use rforces, only: armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn
      use scalars, only: nrzt, nznt, mns, ns, neqs, iter1, iter2, &
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

      real(kind=dp) :: xm(mnmax) ! for handing over to vacuum()
      real(kind=dp) :: xn(mnmax) ! for handing over to vacuum()
      real(kind=dp) :: rmnc(mnmax)
      real(kind=dp) :: zmns(mnmax)
      real(kind=dp) :: lmns(mnmax)
      real(kind=dp) :: rax(nznt)
      real(kind=dp) :: zax(nznt)
      real(kind=dp) :: lu(2*nrztd)
      real(kind=dp) :: lv(2*nrztd)
      integer       :: lodd, l, lk
      real(kind=dp) :: bz0
      real(kind=dp) :: timeon, timeoff

      real(kind=dp), pointer :: worka(12*nrztd)
      real(kind=dp), pointer :: workb(12*nrztd)
      real(kind=dp), pointer :: guu(nrztd)
      real(kind=dp), pointer :: guv(nrztd)
      real(kind=dp), pointer :: gvv(nrztd)
      real(kind=dp), pointer :: czmn(*)
      real(kind=dp), pointer :: crmn(*)

      worka => armn
      workb => r1
      guu   => rcon(1+nrztd)
      guv   => zcon(1+nrztd)
      gvv   => z1
      czmn  => lu
      crmn  => lv






      lodd = 1+nrzt

      ! EXTRAPOLATE M>2 MODES AT JS = 2
      call extrap(xc,xc(1+mns),xc(1+2*mns),xc(1+3*mns),xc(1+4*mns),xc(1+5*mns),xrz3,xrz4,ns)

      ! temporary re-use of gc for scaled xc
      do l = 1,neqs
        gc(l) = xc(l) * scalxc(l)
      enddo

      ! INVERSE FOURIER TRANSFORM TO S,THETA,ZETA SPACE
      ! R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
      do l = 1,nrzt
        lu(l) = c1p0
        lv(l) = czero
        lu(l+nrzt) = czero
        lv(l+nrzt) = czero
      enddo
      call totzsp(gc,gc(1+mns),gc(1+2*mns),gc(1+3*mns),gc(1+4*mns),gc(1+5*mns), &
                  r1,ru,       rv,         z1,         zu,         zv,          &
                  lu,lv,rcon,zcon, &
                  worka,worka,worka, workb)

      ! COMPUTE CONSTRAINT FORCE (GCON)
      do l = 1,nrzt
        ru0(l)  = ru(l) + ru(l+nrzt)*sqrts(l)
        zu0(l)  = zu(l) + zu(l+nrzt)*sqrts(l)
        rcon(l) = rcon(l) + rcon(l+nrzt)*sqrts(l)
        zcon(l) = zcon(l) + zcon(l+nrzt)*sqrts(l)

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
        call dcopy(nrzt,rcon,1,rcon0,1)
        call dcopy(nrzt,zcon,1,zcon0,1)
      endif

      if (iter2.gt.1) &
        call alias(gcon, azmn, worka, worka, worka, gc, gc(1+mns))

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
      call pressure(azmn(lodd), bzmn(lodd))
      if (iter2.eq.1) then
        voli = (twopi**2)*hs*dsum(ns-1,vp(2),1)
      end if

      ! COMPUTE COVARIANT COMPONENTS OF B, MAGNETIC PRESSURE,
      ! AND METRIC ELEMENTS ON HALF-GRID
      do lk = 1,nznt
        rax(lk) = r1(1+ns*(lk-1))
        zax(lk) = z1(1+ns*(lk-1))
      enddo

      if (iequi.eq.1) &
        call dcopy(nrzt, z1, 1, gcon, 1)

      ! see above for re-use of arrays; additionally:
      !           bsubu  bsubv  gsqrt       bsq         r12         rs    zs
      call bcovar(clmn,  blmn,  azmn(lodd), bzmn(lodd), armn(lodd), bzmn, brmn, &
                  azmn,  armn,  guu,        guv,        gvv,        brmn(lodd), lu, lv)
      !           ru12   zu12                                       phipog

      if (iequi.eq.1) &
        call dcopy(nrzt, gcon, 1, z1, 1)

      bz0 = sdot(nznt,blmn(ns),ns,wint(ns),ns)

      ! COMPUTE VACUUM MAGNETIC PRESSURE AT PLASMA EDGE
      ! NOTE: FOR FREE BOUNDARY RUNS, THE PLASMA VOLUME CAN
      ! BE INCREASED BY DECREASING CURPOL AT FIXED PHIPS.  THE
      ! VALUE FOR CURPOL GIVEN BELOW SHOULD BE USED AS AN
      ! INITIAL GUESS, WHICH CAN BE CHANGED BY READING IN
      ! CURPOL FROM DATA STATEMENT
      if (nvac.ne.0.and.iter2.gt.1) then
        call second(timeon)

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
                   c1p5*sdot(nznt,clmn(ns  ),ns,wint(ns),ns)   &
                 - cp5 *sdot(nznt,clmn(ns-1),ns,wint(ns),ns) )

          call convert(rmnc, zmns, lmns, xm, xn, ns,          &
                       xc,          xc(1+  mns), xc(1+2*mns), &
                       xc(1+3*mns), xc(1+4*mns), xc(1+5*mns) )

          call vacuum(rmnc, zmns, xm, xn, ctor, rbtor, bsqvac, rax, zax)
        endif

        do lk=1,nznt
          bsqsav(lk,3) = c1p5*bzmn(ns*lk+nrzt) - cp5*bzmn(ns*lk-1+nrzt)
          rbsq(lk)     = bsqvac(lk) * ohs*(r1(ns*lk) + r1(ns*lk+nrzt))
          dbsq(lk)     = abs(bsqvac(lk)-bsqsav(lk,3))
        enddo

        if(ivac.eq.1) then
          ! copy initial plasma and vacuum field into bsqsav(:,1:2)
          call scopy(nznt,bzmn(ns+nrzt),ns,bsqsav(1,1),1)
          call scopy(nznt,bsqvac,        1,bsqsav(1,2),1)
        endif

        call second(timeoff)
        timer(1) = timer(1) + (timeoff-timeon)
      endif ! vacuum contribution

      if(iequi.eq.1)then
        ! COMPUTE REMAINING COVARIANT COMPONENT OF B (BSUBS),
        ! CYLINDRICAL COMPONENTS OF B (BR, BPHI, BZ), AND
        ! AVERAGE EQUILIBRIUM PROPERTIES AT END OF RUN

        call bss(armn(lodd), bzmn, brmn, azmn, armn, shalf, crmn(lodd),
                 lu, lv, rcon, czmn(lodd), zcon, cp25, cp5, nrzt)

        call eqfor(clmn, blmn, bzmn(lodd), xc, xc(1+2*mns))

        call wrout(bzmn(lodd), azmn(lodd), clmn, blmn, &
                   crmn(lodd), rcon, czmn(lodd), zcon, lu, lv)

        return
      endif

      ! COMPUTE MHD FORCES ON INTEGER-MESH
      call forces(rbsq,guu,guv,gvv,sqrts,shalf,ohs,cp25,cp5,
                  czero,nrzt,ns,ivac)

      ! FOURIER-TRANSFORM MHD FORCES TO (M,N)-SPACE
      do l = 1,neqs
        gc(l) = czero
      enddo
      call tomnsp(gc,gc(1+mns),gc(1+2*mns),gc(1+3*mns),gc(1+4*mns),gc(1+5*mns), &
                  armn,brmn,crmn, azmn,bzmn,czmn,&
                  blmn,clmn, rcon,zcon, workb,workb,workb)

      do l = 1,neqs
        gc(l) = gc(l) * scalxc(l)
      enddo

      ! COMPUTE FORCE RESIDUALS
      call residue(gc,gc(1+2*mns),gc(1+4*mns),bz0,workb)

      return

end
