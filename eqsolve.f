
        subroutine eqsolve(nsval,intflag,ierflag)
      include 'name1'
      include 'name0'
      include 'name2'
        real(4) :: t4
        data w1,r01,res1/0.,0.,1./
************************************************************************
*                 INDEX OF LOCAL VARIABLES
*
*         hs      radial mesh size increment
*         iequi   counter used to call EQFOR at end of run
*         ijacob  counter for number of times jacobian changes sign
*         irst    counter monitoring sign of jacobian; resets R, Z, and
*                 Lambda when jacobian changes sign and decreases time s
*         isigng  sign of Jacobian : must be =1 (right-handed) or =-1
*(left-handed)
*         iterj   stores position in main iteration loop (j=1,2)
*         itfsq   counter for storing FSQ into FSQT for plotting
*         ivac    counts number of free-boundary iterations
*         ndamp   number of iterations over which damping is averaged
*         meven   parity selection label for even poloidal modes of R an
*         modd    parity selection label for odd poloidal modes of R and
*         rb      boundary coefficient array for R (rcc,rss)
*         zb      boundary coefficient array for Z (zcs,zsc)
*         gc      stacked array of R, Z, Lambda Spectral force coefficie
*         xc      stacked array of scaled R, Z, Lambda Fourier coefficie
*         STACKING ORDER:
*           1:mns   => rmncc Fourier coefficients
*           1+mns,2*mns => rmnss Fourier coefficients
*           1+2*mns,3*mns => zmncs Fourier coefficients
*           1+3*mns,4*mns => zmnsc Fourier coefficients
*           1+4*mns,5*mns => lmncs Fourier coefficients
*           1+5*mns,neqs => lmnsc Fourier coefficients
************************************************************************
*************
*                 INITIALIZE MESH-DEPENDENT SCALARS
*************
        ns = nsval
        hs = c1p0/real(ns-1)
        ohs = c1p0/hs
        mns = ns*mnd
        neqs = 6*mns
        nrzt = nznt*ns
        iter2 = 1
        irst = 1
*       DELT = C1P1
*sgg
*       DELT=0.9
        DELT=0.6
        ijacob = -1
        print   *,'NS = ',ns,' NO. FOURIER MODES = ',mnmax
        write(3,*)'NS = ',ns,' NO. FOURIER MODES = ',mnmax
*************
*                 COMPUTE INITIAL R, Z AND MAGNETIC FLUX PROFILES
*                 AND STORE XC, XCDOT FOR POSSIBLE RESTART
*************
        call profil3d(xc,xc(1+2*mns),intflag)
 10     call restart
        ijacob = ijacob+1
        iter1 = iter2
*************
*                 FORCE ITERATION LOOP
*************
        do 30 miter = iter1,niter
*************
*                 ADVANCE FOURIER AMPLITUDES OF R, Z, AND LAMBDA
*************
        iter2 = miter
        call evolve(ierflag)
        if(ierflag.ne.0)return
*************
*                 COMPUTE ABSOLUTE STOPPING CRITERION
*************
        r00 = xc(1)*mscale(0)*nscale(0)
        w0 = wb+wp/(gam-c1p0)
        wdota = abs(w0-w1)/w0
        r0dot = abs(r00-r01)/r00
        r01 = r00
        w1 = w0
        if( (fsqr.le.ftol.and.fsqz.le.ftol.and.fsql.le.ftol
     >      .and.miter-iter1.gt.10)
     >      .or.(ns.lt.nsd.and.(fsqr+fsqz+fsql).le.c2pm8))goto 40
        if(ijacob.ge.100.or.miter.eq.niter)goto 40
************SGG-RESTZEIT-a****************
*       iendsgg=itime(0)
*       IF(IENDSGG.LT.20)GOTO 40
************SGG-RESTZEIT-e****************
        if(ivac.eq.1)then
           print 110, miter,c1p0/bscale
           write( 3,110)miter,c1p0/bscale
           ivac = ivac+1
        endif
        if(mod(miter,(niter/100+1)).ne.0.or.(ns.lt.nsd))goto 15
        itfsq=itfsq+1
        fsqt(itfsq) = fsqr + fsqz
        wdot(itfsq)=max(wdota,c1pm13)
*************
*                 TIME STEP CONTROL
*************
 15     if(miter.eq.iter1)res0=fsq
        res0=min(res0,fsq)
        if( (fsq.le.res0) .and. ((iter2-iter1).gt.nsd/2) )then
          call restart
        else if(fsq.gt.c100p*res1.and.miter.gt.iter1)then
          irst = 2
        else if(fsq.gt.(c1p4*res1).and.mod(miter-iter1,5).eq.0.and.
     >  (miter-iter1).gt.2*ns4)then
          irst = 3
        endif
        if(irst.ne.1)goto 10
        res1=res0
*************
*                 PRINTOUT EVERY NSTEP ITERATIONS
*************
 20     if(mod(miter,nstep).ne.0.and.miter.gt.1)goto 30
        call printout(iter2,w1,r00)
 30     continue
 40     call printout(iter2,w0,r00)
        call second(t4)
        timer(0) = real(t4)
        if(ijacob.ge.100)ierflag = 2
        if(ns.lt.nsd)write(3,100)timer(0),ijacob
        write(3,60)wdota,r0dot
 60     format(/,' d(ln W)/dt = ',1pe10.3,' d(ln R0)/dt = ',1pe10.3,/)
 100    format(/,' COMPUTATIONAL TIME FOR INITIAL INTERPOLATION = ',
     >  1pe10.2,' SECONDS ',' JACOBIAN RESETS = ',i4)
 110    format(/,'  VACUUM PRESSURE TURNED ON AT ',i4, ' ITERATIONS',/,
     >  '  MAGNETIC FIELD IN PLASMA IS SCALED AT END BY ',1pe10.2,
     >  ' (VAC. UNITS)',/)
        return
        end
