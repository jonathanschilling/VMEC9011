        subroutine funct3d
      include 'name1'
      include 'name0'
      include 'name2'
      include 'name3'
      include 'name4'
        real guu(nrztd),guv(nrztd),gvv(nrztd),lu(2*nrztd),lv(2*nrztd),
     >  rmnc(mnmax),zmns(mnmax),lmns(mnmax),xm(mnmax),
     >  xn(mnmax),rax(nznt),zax(nznt),worka(12*nrztd),workb(12*nrztd)
        equivalence(worka,armn),(workb,r1),(guu,rcon(1+nrztd)),
     >  (guv,zcon(1+nrztd)),(gvv,z1),(czmn,lu),(crmn,lv)
        lodd = 1+nrzt
*************
*                 EXTRAPOLATE M>2 MODES AT JS = 2
*************
        call extrap(xc,xc(1+mns),xc(1+2*mns),xc(1+3*mns),
     >  xc(1+4*mns),xc(1+5*mns),xrz3,xrz4,ns)
        do 5 l = 1,neqs
 5      gc(l) = xc(l) * scalxc(l)
*************
*                 INVERSE FOURIER TRANSFORM TO S,THETA,ZETA SPACE
*                 R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
*************
CDIR$ IVDEP
        do 10 l = 1,nrzt
        lu(l) = c1p0
        lv(l) = czero
        lu(l+nrzt) = czero
 10     lv(l+nrzt) = czero
        call totzsp(gc,gc(1+mns),gc(1+2*mns),gc(1+3*mns),gc(1+4*mns),
     >  gc(1+5*mns),r1,ru,rv,z1,zu,zv,lu,lv,rcon,zcon,
     >  worka,worka,worka,workb)
*************
*                 COMPUTE CONSTRAINT FORCE (GCON)
*************
CDIR$ IVDEP
        do 20 l = 1,nrzt
        ru0(l)  = ru(l) + ru(l+nrzt)*sqrts(l)
        zu0(l)  = zu(l) + zu(l+nrzt)*sqrts(l)
        rcon(l) = rcon(l) + rcon(l+nrzt)*sqrts(l)
        zcon(l) = zcon(l) + zcon(l+nrzt)*sqrts(l)
        gcon(l) = czero
 20     azmn(l) = (rcon(l)-rcon0(l))*ru0(l) + (zcon(l)-zcon0(l))*zu0(l)
        do 25 l = 1,2*mns
 25     gc(l) = czero
        if(iter2.eq.1.and.nvac.eq.0)call scopy(nrzt,rcon,1,rcon0,1)
        if(iter2.eq.1.and.nvac.eq.0)call scopy(nrzt,zcon,1,zcon0,1)
        if(iter2.gt.1)
     >  call alias(gcon,azmn,worka,worka,worka,gc,gc(1+mns))
*************
*                 COMPUTE S AND THETA DERIVATIVE OF R AND Z
*                 AND JACOBIAN ON HALF-GRID
*************
        call jacobian(r1,ru,z1,zu,armn,azmn,brmn,bzmn,azmn(lodd),
     >  armn(lodd),brmn(lodd),wint,shalf,ohs,cp25,cp5,czero,
     >  nrzt,nznt,irst,meven,modd)
        if(irst.eq.2.and.iequi.eq.0)return
*************
*                 COMPUTE PRESSURE AND VOLUME ON HALF-GRID
*************
        call pressure(azmn(lodd),bzmn(lodd),wint,wp,dnorm,
     >  mass,vp,pres,gam,nznt,ns)
        if(iter2.eq.1)voli = (twopi**2)*hs*ssum(ns-1,vp(2),1)
*************
*                 COMPUTE COVARIANT COMPONENTS OF B, MAGNETIC
*                 PRESSURE, AND METRIC ELEMENTS ON HALF-GRID
*************
        do 30 lk = 1,nznt
        rax(lk) = r1(1+ns*(lk-1))
 30     zax(lk) = z1(1+ns*(lk-1))
        if( iequi.eq.1 )call scopy(nrzt,z1,1,gcon,1)
        call bcovar(clmn,blmn,azmn(lodd),bzmn(lodd),armn(lodd),bzmn,
     >  brmn,azmn,armn,guu,guv,gvv,brmn(lodd),lu,lv)
        if( iequi.eq.1 )call scopy(nrzt,gcon,1,z1,1)
        bz0 = sdot(nznt,blmn(ns),ns,wint(ns),ns)
*************
*                 COMPUTE VACUUM MAGNETIC PRESSURE AT PLASMA EDGE
*                 NOTE: FOR FREE BOUNDARY RUNS, THE PLASMA VOLUME CAN
*                 BE INCREASED BY DECREASING CURPOL AT FIXED PHIPS.  THE
*                 VALUE FOR CURPOL GIVEN BELOW SHOULD BE USED AS AN
*                 INITIAL GUESS, WHICH CAN BE CHANGED BY READING IN
*                 CURPOL FROM DATA STATEMENT
*************
        if(nvac.eq.0.or.iter2.le.1)goto 55
        call second(timeon)
        ivac2=mod(iter2-iter1,nvacskip)
        if( (fsqr+fsqz).le.1.e-2 )ivac=ivac+1
        if( ivac.eq.1 )ivac2 = 0
        rbtor = curpol
        if(ivac.eq.1.and.rbtor.eq.0.)rbtor = twopi*dnorm*bz0
        if( ivac.ge.1 )then
        ctor = isigng*twopi*dnorm*(
     >         c1p5*sdot(nznt,clmn(ns),ns,wint(ns),ns)
     >       -  cp5*sdot(nznt,clmn(ns-1),ns,wint(ns),ns) )
        call convert(rmnc,zmns,lmns,xm,xn,ns,xc,xc(1+mns),
     >  xc(1+2*mns),xc(1+3*mns),xc(1+4*mns),xc(1+5*mns))
        call vacuum(rmnc,zmns,xm,xn,ctor,rbtor,bsqvac,rax,zax)
          endif
        do 40 lk=1,nznt
        bsqsav(lk,3) = c1p5*bzmn(ns*lk+nrzt) - cp5*bzmn(ns*lk-1+nrzt)
        rbsq(lk) = bsqvac(lk)*ohs*(r1(ns*lk) + r1(ns*lk+nrzt))
 40     dbsq(lk)=abs(bsqvac(lk)-bsqsav(lk,3))
        if(ivac.eq.1)call scopy(nznt,bzmn(ns+nrzt),ns,bsqsav(1,1),1)
        if(ivac.eq.1)call scopy(nznt,bsqvac,1,bsqsav(1,2),1)
        call second(timeoff)
        timer(1) = timer(1) + (timeoff-timeon)
*************
*                 COMPUTE REMAINING COVARIANT COMPONENT OF B (BSUBS),
*                 CYLINDRICAL COMPONENTS OF B (BR, BPHI, BZ), AND
*                 AVERAGE EQUILIBRIUM PROPERTIES AT END OF RUN
*************
 55       if(iequi.eq.1)then
        call bss(armn(lodd),bzmn,brmn,azmn,armn,shalf,crmn(lodd),
     >  lu,lv,rcon,czmn(lodd),zcon,cp25,cp5,nrzt)
        call eqfor(clmn,blmn,bzmn(lodd),xc,xc(1+2*mns))
        call wrout(bzmn(lodd),azmn(lodd),clmn,blmn,crmn(lodd),rcon,
     >  czmn(lodd),zcon,lu,lv)
        return
          endif
*************
*                 COMPUTE MHD FORCES ON INTEGER-MESH
*************
        call forces(rbsq,guu,guv,gvv,sqrts,shalf,ohs,cp25,cp5,
     >  czero,nrzt,ns,ivac)
*************
*                 FOURIER-TRANSFORM MHD FORCES TO (M,N)-SPACE
*************
        do 75 l = 1,neqs
 75     gc(l) = czero
        call tomnsp(gc,gc(1+mns),gc(1+2*mns),gc(1+3*mns),gc(1+4*mns),
     >  gc(1+5*mns),armn,brmn,crmn,azmn,bzmn,czmn,blmn,clmn,
     >  rcon,zcon,workb,workb,workb)
        do 90 l = 1,neqs
 90     gc(l) = gc(l) * scalxc(l)
*************
*                 COMPUTE FORCE RESIDUALS
*************
        call residue(gc,gc(1+2*mns),gc(1+4*mns),bz0,workb)
        return
        end
