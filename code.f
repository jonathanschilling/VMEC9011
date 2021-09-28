        program vmec
************************************************************************
*                               D * I * S * C * L * A * I * M * E * R
*
*       You are using a BETA version of the program VMEC, which is curre
*       under development by S. P. Hirshman at the Fusion Energy Divisio
*       Oak Ridge National Laboratory.  Please report any problems or co
*       to him.  As a BETA version, this program is subject to change
*       and improvement without notice.
*
*       THIS PROGRAM - VMEC (Variational Moments Equilibrium Code)  -
*       SOLVES THREE-DIMENSIONAL MHD EQUILIBRIUM EQUATIONS USING
*       FOURIER SPECTRAL (MOMENTS) METHODS. A CYLINDRICAL COORDINATE
*       REPRESENTATION IS USED (R-Z COORDINATES). THE POLOIDAL
*       ANGLE VARIABLE IS RENORMALIZED THROUGH THE STREAM FUNCTION
*       LAMBDA, WHICH IS SELF-CONSISTENTLY DETERMINED AND DIFFERENCED
*       VARIATIONALLY ON THE FULL-RADIAL MESH. THE POLOIDAL ANGLE IS
*       DETERMINED BY MINIMIZING <M> = m**2 S(m) , WHERE S(m) =
*       Rm**2 + Zm**2 . AN EVEN-ODD DECOMPOSITION IN THE POLOIDAL MODE
*       NO. OF R,Z, AND LAMDA IS USED TO IMPROVE RADIAL RESOLUTION.
*       A FREE-BOUNDARY OPTION IS AVAILABLE (FOR nvac > 0), WITH A
*       USER-SUPPLIED SUBROUTINE "VACUUM" NEEDED TO COMPUTE THE PLASMA
*       BOUNDARY VALUE OF B**2.
*
*       Added features since last edition
*       1.  Implemented preconditioning algorithm for R,Z
*       2.  The physical (unpreconditioned) residuals are used
*           to determine the level of convergence
*       3.  The original (MOMCON) scaling of lambda is used, i.e.,
*           Bsupu = phip*(iota - lamda[sub]v)/sqrt(g). This is needed to
*           maintain consistency with the time-stepper for arbitrary PHI
*
*       WRITTEN BY S. P. HIRSHMAN (8/28/85 - REVISED 3/1/86) BASED ON
*       1. S. P. Hirshman and J. C. Whitson, Phys. Fluids 26, 3553 (1983
*       2. S. P. Hirshman and H. K. Meier, Phys. Fluids 28, 1387 (1985).
*       3. S. P. Hirshman and D. K. Lee, Comp. Phys. Comm. 39, 161 (1986
************************************************************************
      include 'name1'
        character*80 werror(0:6)
        character*24 date,mach,timeloc
C       CALL LINK(' //')
C       OPEN(UNIT=7,FILE='INDATA',STATUS='OLD',ERR=902)
C       OPEN(UNIT=3,FILE='THREED1',STATUS='NEW',ERR=901)
C       OPEN(UNIT=8,FILE='WOUT',STATUS='UNKNOWN',ERR=903)
        goto 904
 901    print *,' THREED1 FILE ALREADY EXISTS: RENAME OR DELETE IT'
        stop
 902    print *,' INDATA FILE IS MISSING'
        stop
 903    print *,' WOUT FILE ALREADY EXIST: RENAME OR DELETE IT'
        stop
 904    continue
************************************************************************
*                 PARAMETER STATEMENT DIMENSIONS AND RECOMMENDED VALUES
*
*         mpol    number of poloidal modes, starting at m=0 (note: mpol
*         nmax    mode number of largest toroidal mode (normed to nfp),
*                 with -nmax <= n <= nmax
*         ntheta  number of theta mesh points (at least 2*mpol+4)
*                 (note: in fft routines in VACUUM, ntheta and nzeta
*                 must be even with no prime factors greater than 5.)
*         nzeta   number of zeta mesh points [at least 2*nmax+4 if nmax>
*                 or =1 if nmax = 0 (axisymmetric case)]
*         nvac    turns on fixed (=0) or free(=1) boundary option
************
*                 INDEX OF LOCAL VARIABLES
************
*         ierflag specifies error condition if nonzero
*         intflag = 0,    compute xc array from boundary data
*                 = nsin, compute xc array by interpolation
************************************************************************

c       call timedate(timeloc,date,mach)
        WRITE(3,*)'THIS IS THE PRECONDITIONED VMEC.FULL CODE: VMEC9011
     >                                      '
************
*               INITIALIZE ALL VARIABLES SO VMEC CAN BE CALLED
*               AS A SUBROUTINE.
************
        ierflag = 0
        intflag = 0
        call vsetup
************
*               READ INPUT FILE INDATA AND STORE IN INPUTDAT COMMON BLOC
************
        call readin(nsin,ierflag)
        if(ierflag.ne.0)goto 10
************
*                 COMPUTE INVARIANT ARRAYS (IF VMEC IS USED AS A
*                 SUBROUTINE, CALL FIX_ARAY ONLY ONCE AT SETUP)
************
        call fixaray(ierflag)
        if(ierflag.ne.0)goto 10
************
*                 MAKE INITIAL CALL TO EQUILIBRIUM SOLVER
************
        if(nsin.lt.2.or.nsin.ge.nsd)nsin = nsd
        call eqsolve(nsin,intflag,ierflag)
************
*                 PERFORM COARSE-TO-FINE MESH INTERPOLATION
************
        if(nsin.lt.nsd.and.ierflag.eq.0)then
          intflag = nsin
          call eqsolve(nsd,intflag,ierflag)
        end if
************
*                 WRITE OUTPUT TO FILES THREED1, WOUT
************
        if(ierflag.eq.0.or.ierflag.eq.1.or.ierflag.eq.5)call output
 10     werror(0)=' EXECUTION TERMINATED NORMALLY '
        werror(1)=' INITIAL JACOBIAN CHANGED SIGN (NEED A BETTER GUESS)'
        werror(2)=' MORE THAN 100 JACOBIAN ITERATIONS (DECREASE DELT)'
        werror(3)=' m IN BOUNDARY ARRAY EXCEEDS mpol-1 (INCREASE MPOL)'
        werror(4)=' n IN BOUNDARY ARRAYS OUTSIDE nmax,(-nmax) RANGE'
        werror(5)=' ASSUMED SIGN OF JACOBIAN (ISIGNG) IS WRONG'
        print *, werror(ierflag)
        write(3,*)werror(ierflag)
        call exit
        end

************************************************************************
*               INITIALIZE COMMON BLOCK CONSTANTS
************************************************************************
        block data
      include 'name1'
      include 'name0'
      include 'name2'
        parameter(mup2=nmax+nmax1,mlo3=mup2+1,
     >  zerod=0.0e+00,oned=1.0e+00,cp707d=7.071e-01)
        common /bounds/ mupper(3),mlower(3)

        data mupper/nmax,mup2,mnd1/, mlower/0,nmax1,mlo3/
        DATA NDAMP/10/, NS4/25/, ISIGNG/-1/,
     >  meven/0/, modd/1/, jmin1/1,1,mpol1*2/, jmin2/1,2,mpol1*3/,
     >  jlam/2,3,mpol1*3/,phips(1)/zerod/,mscale/cp707d,mpol1*oned/,
     >  nscale/cp707d,nmax1*oned/,cp15/1.5e-01/,
     >  cp25/2.5e-01/,cp5/5.0e-01/,cp9/9.0e-01/,c2pm8/2.e-08/,
     >  cp96/9.6e-01/,c1p0/oned/,c1p1/1.1e+00/,c1p4/1.4e+00/,
     >  c1p5/1.5e+00/,c2p0/2.0e+00/,c1pm13/1.0e-13/,
     >  c3p0/3.0e+00/,c8p0/8.0e+00/,c100p/1.00e+02/,czero/zerod/
        end

        subroutine eqsolve(nsval,intflag,ierflag)
      include 'name1'
      include 'name0'
      include 'name2'
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
c          print 110, miter,c1p0/bscale
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
        call second(timer(0))
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

        subroutine output
      include 'name1'
      include 'name0'
      include 'name2'

        iequi = 1
        call funct3d
*       PRINT 10, IJACOB
        write(3,10)ijacob
 10     format(/,'  NUMBER OF JACOBIAN RESETS = ',i4)
*       PRINT 20, TIMER(0),TIMER(1)
        write(3,20)timer(0),timer(1)
 20     format (/,
     >  '  TOTAL COMPUTATIONAL TIME :       ',1pe10.2,' SECONDS',/
     >  '  TIME IN VACUUM LOOP :            ',1pe10.2,' SECONDS',/)
        return
        end

        subroutine vsetup
      include 'name1'
      include 'name0'
      include 'name2'
      include 'name3'

        do 10 l = 1,nrztd
        rcon0(l) = czero
 10     zcon0(l) = czero
        do 20 i = 0,10
 20     timer(i) = czero
        iequi = 0
        itfsq = 0
        bscale = c1p0
        delbsq = c1p0
        ivac   = 0
        return
        end

        subroutine readin(nsin,ierflag)
      include 'name1'
      include 'name0'
      include 'name2'
        real rc(0:mpol1,-nmax:nmax),zs(0:mpol1,-nmax:nmax),
     >  rmag(-nmax:nmax),zmag(-nmax:nmax)
************************************************************************
*                 INPUT DATA AND RECOMMENDED VALUES
*
*         ftol    value of fsq = fr**2 + fz**2 at which iteration ends
*         gam     value of compressibility index (gam=0 => pressure pres
*         ncurr   flux conserving (=0) or no net toroidal current (=1)
*         nfp     number of toroidal field periods
*         niter   number of iterations (used to terminate run)
*         nsin    number of radial intervals used to obtain initial gues
*                 for interpolation onto finer grid; <2, no interpolatio
*         ns      number of radial mesh points (>2)
*         nstep   number of timesteps between printouts on screen
*         rc      boundary coefficients of cos(m*theta-n*zeta) for R
*         zs      boundary coefficients of sin(m*theta-n*zeta) for Z
*
*               INITIALIZATION BLOCK (NEEDED IF VMEC CALLED REPETITIVELY
*
************************************************************************
        do 5 n = -nmax,nmax
        rmag(n) = czero
        zmag(n) = czero
        do 5 m = 0,mpol1
        rc(m,n) = czero
 5      zs(m,n) = czero
        do 10 n = 0,nmax
        raxis(n) = czero
        zaxis(n) = czero
        do 10 m = 0,mpol1
        rb(n,m,1) = czero
        rb(n,m,2) = czero
        zb(n,m,1) = czero
 10     zb(n,m,2) = czero
************************************************************************
*               READ IN DATA FROM INDATA FILE
************************************************************************
        read(5,*)
        read(5,*)nfp,mnbound,ncurr,nsin,niter,nstep,ftol
        read(5,*)
        read(5,*)gam,phiedge,curpol,curtor
        write(3,15)nsd,ntheta1,nzeta,mpol,nmax,nfp,gam,
     >  phiedge,curpol,curtor
 15     format(/,' COMPUTATION PARAMETERS: (u = theta, v = zeta)'//,
     >  '     ns     nu     nv     mu     mv'//,5i7,//,
     >  ' CONFIGURATION PARAMETERS:'//,'    nfp      gamma',
     >  '    phiedge     curpol     curtor',//,i7,1p4e11.3)
        nvacskip = nfp
        write(3,20)ncurr,niter,nsin,nstep,nvacskip,ftol
 20     format(/' CONTROL PARAMETERS:'//,'  ncurr  niter',
     >  '   nsin  nstep  nvacskip      ftol',//,4i7,i10,1pe10.2,/)
        write(3,25)
 25     format(' MASS PROFILE EXPANSION COEFFICIENTS (am):',/)
        read(5,*)
        read(5,*)(am(i),i=0,5)
        read(5,*)
        write(3,40)(am(i),i=0,5)
        write(3,45)
        read(5,*)(ai(i),i=0,5)
        read(5,*)
        write(3,40)(ai(i),i=0,5)
        WRITE(3,46)
        READ(5,*)(AC(I),I=1,5)
        READ(5,*)
        WRITE(3,40)(AC(I),I=1,5)
 40     format(1p6e12.3,/)
 45     format(' IOTA PROFILE EXPANSION COEFFICIENTS (ai):',/)
 46     format(' CURRENT PROFILE EXPANSION COEFFICIENTS (ac):',/)
*       DO 50 I=1,MNBOUND
        DO 50 I=1, MPOL1*(2*NMAX+1)
        READ(5,*,END=51)M,N,RC(M,N),ZS(M,N)
        write(*,*)M,N,RC(M,N),ZS(M,N)
        if(m.gt.mpol1.or.m.lt.0) ierflag = 3
        if(n.gt.nmax.or.n.lt.(-nmax)) ierflag = 4
        if(ierflag.ne.0)return
        if(m.eq.0)backspace 5
        if(m.eq.0) read(5,*)m,n,rc(m,n),zs(m,n),rmag(n),zmag(n)
 50     continue
 51     WRITE(3,55)
*A_NORMIERUNG SGG                                                               
        IF((RC(1,0)+ZS(1,0)).EQ.0)STOP'KANN FOURKOEF. NICHT NORMIEREN'          
        FOURNORM=2./(abs(RC(1,0))+abs(ZS(1,0)))  
*       FOURNORM=1.
         DO 47 I=0,MPOL1                                                         
            DO 47 J=-NMAX,NMAX                                                   
               RC(I,J)=RC(I,J)*FOURNORM                                          
               ZS(I,J)=ZS(I,J)*FOURNORM                                          
  47     CONTINUE                                                                
         DO 48 J=-NMAX,NMAX                                                      
           RMAG(J)=RMAG(J)*FOURNORM                                              
           ZMAG(J)=ZMAG(J)*FOURNORM                                              
  48     CONTINUE                                                                
*E_NORMIERUNG SGG                                                               
 55     format
     >  (/,'   mb  nb     rbc         zbs        raxis       zaxis',/)
        do 60 m=0,mpol1
        do 60 n=-nmax,nmax
        n1 = iabs(n)
        isgn = 1
        if(n .lt. 0)isgn = -1
        if(m.eq.0)raxis(n1)=raxis(n1) + rmag(n)
        if(m.eq.0)zaxis(n1)=zaxis(n1) - zmag(n)*isgn
        rb(n1,m,1) = rb(n1,m,1) + rc(m,n)
        rb(n1,m,2) = rb(n1,m,2) + isgn*rc(m,n)
        zb(n1,m,1) = zb(n1,m,1) - isgn*zs(m,n)
        zb(n1,m,2) = zb(n1,m,2) + zs(m,n)
        if(n.eq.0 .or. m.eq.0)rb(n1,m,2) = czero
        if(n.eq.0)zb(n1,m,1) = czero
        if(m.eq.0)zb(n1,m,2) = czero
        if(rc(m,n).eq.czero.and.zs(m,n).eq.czero)goto 60
        if(m.eq.0)write(3,65)m,n,rc(m,n),zs(m,n),rmag(n),zmag(n)
        if(m.ne.0)write(3,65)m,n,rc(m,n),zs(m,n)
 60     continue
 65     format(i5,i4,1p4e12.4)
        WRITE(3,70)
 70     format(/,' LEGEND: FSQR, FSQZ = Normalized Physical Force ',
     >  'Residuals;    fsqr, fsqz = Preconditioned Force Residuals'/)
        return
        end

        subroutine fixaray(ierflag)
      include 'name1'
      include 'name0'
      include 'name2'
*************
*                 INDEX OF LOCAL VARIABLES
*************
*         mscale  array for norming theta-trig functions (internal use o
*         nscale  array for norming zeta -trig functions (internal use o
*************
*                 COMPUTE TRIGONOMETRIC FUNCTION ARRAYS
*************
        twopi = c8p0*atan(c1p0)
        dnorm = c1p0/real(nzeta*(ntheta2-1))

        do 10 i = 1,ntheta2
        argi = twopi*real(i-1)/real(ntheta1)
        do 10 m = 0,mpol1
        arg = argi*real(m)
        cosmui(i,m) = cos(arg)*mscale(m)
        if(i.eq.1.or.i.eq.ntheta2)cosmui(i,m) = cp5*cosmui(i,m)
        sinmu (i,m) = sin(arg)*mscale(m)
        cosmu (i,m) = cos(arg)*mscale(m)
        sinmum(i,m) =-sinmu(i,m)*real(m)
        cosmum(i,m) = cosmu(i,m)*real(m)
        cosmumi(i,m)= cosmui(i,m)*real(m)
 10     continue

        do 20 j = 1,nzeta
        argj = twopi*real(j-1)/real(nzeta)
        do 20 n = 0,nmax
        arg = argj*real(n)
        cosnv (j,n) = cos(arg)*nscale(n)
        sinnv (j,n) = sin(arg)*nscale(n)
        cosnvn(j,n) = cosnv(j,n)*real(n*nfp)
        sinnvn(j,n) =-sinnv(j,n)*real(n*nfp)
 20     continue

        do 30 m = 0,mpol1
        power = -cp5*real(m)
        xmpq(m,1) = real(m*(m-1))
        xmpq(m,2) = real(m**4)
        xmpq(m,3) = real(m**5)
        t1 = -cp5*real(isigng)*dnorm/real((1+m)**4)
        do 30 n = 0,nmax
        faccon(n,m) = t1/(mscale(m)*nscale(n))**2
        if(m.eq.0.or.n.eq.0)faccon(n,m) = cp5*faccon(n,m)
        xrz3(n,m) = c2p0**(power + c1p0)
 30     xrz4(n,m) =-c3p0**power
*************
*                 CHECK SIGN OF JACOBIAN (SHOULD BE SAME AS ISIGNG)
*************
        rtest = ssum(nmax1,rb(0,1,1),1)
        ztest = ssum(nmax1,zb(0,1,2),1)
        if( (rtest*ztest*real(isigng)) .ge. czero)ierflag = 5
        return
        end

        subroutine profil3d(rmn,zmn,intflag)
      include 'name1'
      include 'name0'
      include 'name2'
        real rmn(ns,0:nmax,0:mpol1,2),zmn(ns,0:nmax,0:mpol1,2)
        pmass(x) = am(0)+x*(am(1)+x*(am(2)+x*(am(3)+x*(am(4)+x*am(5)))))
        piota(x) = ai(0)+x*(ai(1)+x*(ai(2)+x*(ai(3)+x*(ai(4)+x*ai(5)))))
        pcurr(x) =       x*(ac(1)+x*(ac(2)+x*(ac(3)+x*(ac(4)+x*ac(5)))))
************************************************************************
*                 INDEX OF LOCAL VARIABLES
*
*         ai      vector of coefficients in phi-series for iota
*         am      vector of coefficients in phi-series for mass
*         ac      vector of coefficients in phi-series for toroidal curr
*         iotas   rotational transform , on half radial mesh
*         jv   (-)toroidal current inside flux surface (vanishes like s)
*         mass    mass profile on half-grid
*         phip    radial derivative of phi/(2*pi) on half-grid
*         phiedge value of real toroidal flux at plasma edge (s=1)
*         phips   same as phip , one-dimensional array
*         pressure pressure profile on half-grid, mass/phip**gam
*         shalf   sqrt(s) ,two-dimensional array on half-grid
*         sqrts   sqrt(s), two-dimensional array on full-grid
*         wint    two-dimensional array for normalizing angle integratio
************************************************************************
*                 VERIFY SIGN OF CURPOL CONSISTENT WITH TOROIDAL FLUX
*************
        if( phiedge*curpol .lt. czero )then
        print *,'CHANGING SIGN OF PHIEDGE '
        write(3,*)'CHANGING SIGN OF PHIEDGE '
        phiedge = -phiedge
        endif
**************
*                 COMPUTE PHIP, MASS (PRESSURE), IOTA PROFILES ON HALF-G
**************
        torcur = real(isigng)*curtor/(dnorm*twopi)
        do 10 js=2,ns
        phij = hs*(js - c1p5)
        phips(js) = real(isigng) * phiedge / twopi
        iotas(js) = piota(phij)
        jv(js) = -real(isigng) * torcur * pcurr(phij)
 10     mass(js) = pmass(phij)*(abs(phips(js))*rb(0,0,1))**gam
        do 15 js=1,ns
        do 15 lk=1,nznt
        loff = js + ns*(lk-1)
        shalf(loff) = sqrt(hs*abs(js - c1p5))
        sqrts(loff) = sqrt(hs*(js-1))
 15     phip(loff) = phips(js)
        lk = 0
        do 20 lt = 1,ntheta2
        do 20 lz = 1,nzeta
        lk = lk+1
        do 20 js = 2,ns
 20     wint(js+ns*(lk-1)) = cosmui(lt,0)/mscale(0)
        do 25 lk = 1,nznt
        shalf(1+ns+ns*(lk-1)) = c1p0
 25     wint(1+ns*(lk-1)) = czero
*************
*                 INITIALIZE XCDOT=0 AND SAVE COARSE-MESH, SCALED
*                 FOURIER COEFFICIENTS IN GC FOR INTERPOLATION
*************
        do 30 l = 1,neqs
        gc(l) = xc(l)
        xc(l) = czero
 30     xcdot(l) = czero
        do 35 l = 1,6*mns
 35     gc(l) = gc(l) * scalxc(l)
*************
*                 COMPUTE INITIAL R AND Z FOURIER COEFFICIENTS,
*                 FROM SCALED BOUNDARY VALUES, AND SCALXC ARRAY
*                 (1/SQRTS FACTOR FOR ODD M'S)
*************
        do 40 ntype = 1,2
        do 40 m = 0,mpol1
        do 40 n = 0,nmax
        t1 = c1p0/(mscale(m)*nscale(n))
        do 40 js = 1,ns
        sm0 = c1p0 - (hs*(js-1))
        l = js + ns*(n + nmax1*m) + (ntype-1)*mns
        if(mod(m,2).eq.1)scalxc(l) = c1p0/max(sqrts(js),sqrts(2))
        if(mod(m,2).eq.0)scalxc(l) = c1p0
        if(m.eq.0.and.ntype.eq.1)then
                rmn(js,n,m,1) = (rb(n,m,1)+(raxis(n)-rb(n,m,1))*sm0)*t1
                zmn(js,n,m,1) = (zb(n,m,1)+(zaxis(n)-zb(n,m,1))*sm0)*t1
        else if(m.eq.0 .and. ntype.eq.2)then
                rmn(js,n,m,2) = czero
                zmn(js,n,m,2) = czero
        else if(m.ne.0)then
                facj = t1*sqrts(js)**m
                rmn(js,n,m,ntype) = rb(n,m,ntype)*facj
                zmn(js,n,m,ntype) = zb(n,m,ntype)*facj
        endif
 40     continue
        do 50 l = 1,4*mns
 50     scalxc(l+2*mns) = scalxc(l)
*************
*                 INTERPOLATE FROM COARSE TO FINE RADIAL GRID
*************
        if(intflag.ne.0)call interp(xc,gc,scalxc,intflag)
        return
        end

        subroutine evolve(ierflag)
      include 'name1'
      include 'name0'
      include 'name2'
*************
*                 COMPUTE MHD FORCES
*************
        call funct3d
        if(iter2.ne.1.or.irst.ne.2)goto 10
        ierflag = 1
        return
*************
*                 COMPUTE DAMPING PARAMETER (DTAU) AND EVOLVE
*                 R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
*************
 10     if(iter2.ne.iter1)goto 15
        ndamp1 = min0(ndamp,15)
        do 35 i = 1,ndamp1
 35     otau(i) = cp15/delt
 15     fsq1 = fsqr1 + fsqz1 + fsql1
        if(iter2.gt.iter1)dtau = min(abs(log(fsq/fsq1)),cp15)
        fsq = fsq1
        if(iter2.le.1)return
        do 25 i = 1,ndamp1-1
 25     otau(i) = otau(i+1)
        if(iter2.gt.iter1)otau(ndamp1) = dtau/delt
        otav = ssum(ndamp1,otau,1)/ndamp1
        dtau = delt*otav
        b1 = c1p0 - cp5*dtau
        fac = c1p0/(c1p0 + cp5*dtau)
        do 20 l = 1,neqs
        xcdot(l) = fac*(xcdot(l)*b1 + delt*gc(l))
 20     xc(l) = xc(l) + xcdot(l)*delt
        return
        end

        subroutine restart
      include 'name1'
      include 'name0'
      include 'name2'
        goto (10,20,20)irst
 10     call scopy(neqs,xc,1,xstore,1)
        return
 20     do 30 l = 1,neqs
        xcdot(l) = czero
 30     xc(l) = xstore(l)
        delt = delt*(cp96*(irst-2) + cp9*(3-irst))
        irst = 1
        return
        end

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

        subroutine extrap(rmncc,rmnss,zmncs,zmnsc,lmncs,lmnsc,x3,x4,ns)
      include 'name1'
        real rmncc(ns,0:mnd1),rmnss(ns,0:mnd1),x3(0:mnd1),
     >       zmncs(ns,0:mnd1),zmnsc(ns,0:mnd1),x4(0:mnd1),
     >       lmncs(ns,0:mnd1),lmnsc(ns,0:mnd1)

        do 10 mn = 2*nmax1,mnd1
        rmncc(2,mn) = x3(mn)*rmncc(3,mn) + x4(mn)*rmncc(4,mn)
        rmnss(2,mn) = x3(mn)*rmnss(3,mn) + x4(mn)*rmnss(4,mn)
        zmncs(2,mn) = x3(mn)*zmncs(3,mn) + x4(mn)*zmncs(4,mn)
 10     zmnsc(2,mn) = x3(mn)*zmnsc(3,mn) + x4(mn)*zmnsc(4,mn)
        do 20 mn = nmax1,mnd1
        lmncs(2,mn) = x3(mn)*lmncs(3,mn) + x4(mn)*lmncs(4,mn)
 20     lmnsc(2,mn) = x3(mn)*lmnsc(3,mn) + x4(mn)*lmnsc(4,mn)
        do 30 n = 1,nmax
 30     lmncs(1,n) = x3(0)*lmncs(2,n) - lmncs(3,n)
        return
        end

        subroutine alias(gcon,zcon,work1,work2,work3,gcs,gsc)
      include 'name1'
      include 'name0'
      include 'name2'
        real gcon(ns*nzeta,ntheta2),zcon(ns*nzeta,ntheta2),
     >  gcs(ns,0:nmax,0:mpol1),gsc(ns,0:nmax,0:mpol1),
     >  work3(ns,nzeta,4),work1(4*ns*nzeta),work2(ns*nzeta,4)

*************
*                 BEGIN DE-ALIASING (TRUNCATION OF GCON IN FOURIER-SPACE
*************
        do 60 m = 1,mpol1-1
        do 20 l = 1,4*ns*nzeta
 20     work1(l) = czero
        do 25 i = 1,ntheta2
        do 25 jk = 1,ns*nzeta
        work2(jk,01) = work2(jk,01) + zcon(jk,i)*cosmui(i,m)
 25     work2(jk,02) = work2(jk,02) + zcon(jk,i)*sinmu (i,m)
        do 30 n = 0,nmax
        fm = faccon(n,m)
        do 30 k = 1,nzeta
        do 30 js= 2,ns
        gcs(js,n,m) =gcs(js,n,m) +fm*tcon(js)*work3(js,k,01)*sinnv(k,n)
 30     gsc(js,n,m) =gsc(js,n,m) +fm*tcon(js)*work3(js,k,02)*cosnv(k,n)
*************
*                 RECONSTRUCT DE-ALIASED GCON
*************
        do 40 n = 0,nmax
        do 40 k = 1,nzeta
        do 40 js= 2,ns
        work3(js,k,03) = work3(js,k,03) + gcs(js,n,m)*sinnv(k,n)
 40     work3(js,k,04) = work3(js,k,04) + gsc(js,n,m)*cosnv(k,n)
        do 50 i = 1,ntheta2
        do 50 jk= 1,ns*nzeta
 50     gcon(jk,i) = gcon(jk,i) + work2(jk,03)*cosmu(i,m)
     >                          + work2(jk,04)*sinmu(i,m)
 60     continue
        return
        end

        subroutine totzsp(rmncc,rmnss,zmncs,zmnsc,lmncs,lmnsc,
     >  r1,ru,rv,z1,zu,zv,lu,lv,rcon,zcon,
     >  work1,work2,work3,realsp)
      include 'name1'
      include 'name0'
      include 'name2'
        real rmncc(ns,0:nmax,0:mpol1),rmnss(ns,0:nmax,0:mpol1),
     >       zmncs(ns,0:nmax,0:mpol1),zmnsc(ns,0:nmax,0:mpol1),
     >       lmncs(ns,0:nmax,0:mpol1),lmnsc(ns,0:nmax,0:mpol1),
     >  r1(ns*nzeta,ntheta2,0:1),ru(ns*nzeta,ntheta2,0:1),
     >  rv(ns*nzeta,ntheta2,0:1),z1(ns*nzeta,ntheta2,0:1),
     >  zu(ns*nzeta,ntheta2,0:1),zv(ns*nzeta,ntheta2,0:1),
     >  lu(ns*nzeta,ntheta2,0:1),lv(ns*nzeta,ntheta2,0:1),
     >  rcon(ns*nzeta,ntheta2,0:1),zcon(ns*nzeta,ntheta2,0:1)
        real realsp(16*nrztd),work3(ns,nzeta,12),work1(ns*nzeta*12)
        real work2(ns*nzeta,12)
*************
*                 THIS ROUTINE ASSUMES THE FOLLOWING STACKING OF R, Z,
*                 LAMBDA ARRAYS:
*                 rmncc(ns,0:nmax,0:mpol1),rmnss,zmncs,zmncc,lmncs,lmnsc
*************
*                 INITIALIZATION BLOCK
*************
        do 10 l = 1,16*nrztd
 10     realsp(l) = czero
*************
*                 EXTRAPOLATION AT JS=1 FOR M=1 MODES
*************
        do 30 n = 0,nmax
        rmncc(1,n,1) = c2p0*rmncc(2,n,1) - rmncc(3,n,1)
        rmnss(1,n,1) = c2p0*rmnss(2,n,1) - rmnss(3,n,1)
        zmncs(1,n,1) = c2p0*zmncs(2,n,1) - zmncs(3,n,1)
        zmnsc(1,n,1) = c2p0*zmnsc(2,n,1) - zmnsc(3,n,1)
        lmncs(1,n,1) = c2p0*lmncs(2,n,1) - lmncs(3,n,1)
 30     lmnsc(1,n,1) = c2p0*lmnsc(2,n,1) - lmnsc(3,n,1)
*************
*                 COMPUTE R, Z, AND LAMBDA IN REAL SPACE
*                 BEGIN INVERSE TRANSFORM IN N-ZETA
*************
        do 70 m = 0,mpol1
        mparity = mod(m,2)
        do 40 l = 1,12*ns*nzeta
 40     work1(l) = czero
        do 50 n = 0,nmax
        do 50 k = 1,nzeta
CDIR$ IVDEP
        do 50 js= jmin1(m),ns
        work3(js,k,1) = work3(js,k,1) + rmncc(js,n,m)*cosnv (k,n)
        work3(js,k,2) = work3(js,k,2) + rmnss(js,n,m)*sinnv (k,n)
        work3(js,k,3) = work3(js,k,3) + rmncc(js,n,m)*sinnvn(k,n)
        work3(js,k,4) = work3(js,k,4) + rmnss(js,n,m)*cosnvn(k,n)
        work3(js,k,5) = work3(js,k,5) + zmncs(js,n,m)*sinnv (k,n)
        work3(js,k,6) = work3(js,k,6) + zmnsc(js,n,m)*cosnv (k,n)
        work3(js,k,7) = work3(js,k,7) + zmncs(js,n,m)*cosnvn(k,n)
        work3(js,k,8) = work3(js,k,8) + zmnsc(js,n,m)*sinnvn(k,n)
        work3(js,k,9) = work3(js,k,9) + lmncs(js,n,m)*sinnv (k,n)
        work3(js,k,10)= work3(js,k,10)+ lmnsc(js,n,m)*cosnv (k,n)
        work3(js,k,11)= work3(js,k,11)+ lmncs(js,n,m)*cosnvn(k,n)
 50     work3(js,k,12)= work3(js,k,12)+ lmnsc(js,n,m)*sinnvn(k,n)
*************
*                 INVERSE TRANSFORM IN M-THETA
*************
        do 60 i = 1,ntheta2
        cosmux = xmpq(m,1)*cosmu(i,m)
        sinmux = xmpq(m,1)*sinmu(i,m)
        do 60 jk = 1, nzeta*ns
        r1(jk,i,mparity) = r1(jk,i,mparity) +
     >      work2(jk,01)*cosmu (i,m) + work2(jk,02)*sinmu (i,m)
        ru(jk,i,mparity) = ru(jk,i,mparity) +
     >      work2(jk,02)*cosmum(i,m) + work2(jk,01)*sinmum(i,m)
        rv(jk,i,mparity) = rv(jk,i,mparity) +
     >      work2(jk,03)*cosmu (i,m) + work2(jk,04)*sinmu (i,m)
        z1(jk,i,mparity) = z1(jk,i,mparity) +
     >      work2(jk,05)*cosmu (i,m) + work2(jk,06)*sinmu (i,m)
        zu(jk,i,mparity) = zu(jk,i,mparity) +
     >      work2(jk,06)*cosmum(i,m) + work2(jk,05)*sinmum(i,m)
        zv(jk,i,mparity) = zv(jk,i,mparity) +
     >      work2(jk,07)*cosmu (i,m) + work2(jk,08)*sinmu (i,m)
        lu(jk,i,mparity) = lu(jk,i,mparity) +
     >      work2(jk,10)*cosmum(i,m) + work2(jk,09)*sinmum(i,m)
        lv(jk,i,mparity) = lv(jk,i,mparity) -
     >     (work2(jk,11)*cosmu (i,m) + work2(jk,12)*sinmu (i,m))
        rcon(jk,i,mparity) = rcon(jk,i,mparity) +
     >      work2(jk,01)*cosmux      + work2(jk,02)*sinmux
 60     zcon(jk,i,mparity) = zcon(jk,i,mparity) +
     >      work2(jk,05)*cosmux      + work2(jk,06)*sinmux
 70     continue
        return
        end

        subroutine jacobian(r1,ru,z1,zu,zu12,ru12,zs,rs,gsqrt,r12,
     >  tau,wint,shalf,ohs,cp25,cp5,czero,nrzt,nznt,irst,meven,modd)
        real r1(nrzt,0:1),ru(nrzt,0:1),z1(nrzt,0:1),zu(nrzt,0:1),rs(*),
     >  zs(*),r12(*),gsqrt(*),ru12(*),zu12(*),shalf(*),tau(*),wint(*)
        implicit real*8 (a-h,o-z)
*************
*                 (RS, ZS)=(R, Z) SUB S, (RU12, ZU12)=(R, Z) SUB THETA(=
*                 AND GSQRT=SQRT(G) ARE DIFFERENCED ON HALF MESH
*************
        do 10 l=2,nrzt
        lm = l-1
        ru12(l) = cp5*(ru(l,meven) + ru(lm,meven) + shalf(l)*
     >                (ru(l,modd ) + ru(lm,modd )))
        zu12(l) = cp5*(zu(l,meven) + zu(lm,meven) + shalf(l)*
     >                (zu(l,modd ) + zu(lm,modd )))
        rs(l)   = ohs*(r1(l,meven) - r1(lm,meven) + shalf(l)*
     >                (r1(l,modd ) - r1(lm,modd )))
        zs(l)   = ohs*(z1(l,meven) - z1(lm,meven) + shalf(l)*
     >                (z1(l,modd ) - z1(lm,modd )))
        r12(l)  = cp5*(r1(l,meven) + r1(lm,meven) + shalf(l)*
     >                (r1(l,modd ) + r1(lm,modd )))
        gsqrt(l)= r12(l) * ( ru12(l)*zs(l) - rs(l)*zu12(l) + cp25*
     >  ( ru(l,modd )*z1(l,modd) + ru(lm,modd )*z1(lm,modd)
     >  - zu(l,modd )*r1(l,modd) - zu(lm,modd )*r1(lm,modd)
     >  +(ru(l,meven)*z1(l,modd) + ru(lm,meven)*z1(lm,modd)
     >  - zu(l,meven)*r1(l,modd) - zu(lm,meven)*r1(lm,modd))/shalf(l)) )
 10     tau(l) = wint(l)*gsqrt(l)
*************
*                 TEST FOR SIGN CHANGE IN JACOBIAN
*************
        taumax = czero
        taumin = czero
        do 20 l=2,nrzt
        taumax = amax0(tau(l),taumax)
 20     taumin = amin0(tau(l),taumin)
        if(taumax*taumin.lt.czero)irst=2
        return
        end

        subroutine pressure(gsqrt,bsq,wint,wp,dnorm,mass,
     >  vp,pres,gam,nznt,ns)
        real mass(*),gsqrt(ns,*),bsq(ns,*),vp(*),pres(*),wint(ns,*)

        do 10 js = 2,ns
 10     vp(js) = sdot(nznt,gsqrt(js,1),ns,wint(js,1),ns)
        do 20 js = 2,ns
        vp(js) = dnorm*abs(vp(js))
 20     pres(js) = mass(js)/vp(js)**gam
        wp = sdot(ns,vp,1,pres,1)/real(ns-1)
        do 30 js = 2,ns
        do 30 lk = 1,nznt
 30     bsq(js,lk) = pres(js)
        return
        end

        subroutine bcovar(bsubu,bsubv,gsqrt,bsq,r12,rs,zs,ru12,zu12,
     >  guu,guv,gvv,phipog,lu,lv)
      include 'name1'
      include 'name0'
      include 'name2'
      include 'name3'
        common/precond/ard(nsd1,2),arm(nsd1,2),brd(nsd1,2),brm(nsd1,2),
     >        cr(nsd1),azd(nsd1,2),azm(nsd1,2),bzd(nsd1,2),bzm(nsd1,2)
        real bsubu(nrzt,0:1),bsubv(nrzt,0:1),ar(nsd),az(nsd),gsqrt(*),
     >  phipog(*),bsq(*),r12(*),ru12(*),zu12(*),rs(*),zs(*),guu(*),
     >  guv(*),gvv(*),lu(nrzt,0:1),lv(nrzt,0:1)
*************
*                 INITIALIZATION BLOCK
*************
        do 10 l = 1,nrzt+1
           guu(l) = czero
           guv(l) = czero
           gvv(l) = czero
 10     continue
*************
*                 COMPUTE METRIC ELEMENTS GIJ ON HALF MESH
*************
        do 20 is = -1,0
           do 20 l=2,nrzt
           lme = l + is
           lmo = lme + nrzt
           phipog(l) = cp5*sqrts(lme)*sqrts(lme)
           guu(l) = guu(l) + cp5*(ru(lme)*ru(lme) + zu(lme)*zu(lme))
     >               + phipog(l)*(ru(lmo)*ru(lmo) + zu(lmo)*zu(lmo))
     >               + shalf(l)*(ru(lme)*ru(lmo) + zu(lme)*zu(lmo))
           guv(l) = guv(l) + cp5*(ru(lme)*rv(lme) + zu(lme)*zv(lme))
     >               + phipog(l)*(ru(lmo)*rv(lmo) + zu(lmo)*zv(lmo))
     >            + cp5*shalf(l)*(ru(lme)*rv(lmo) + rv(lme)*ru(lmo)
     >                         + zu(lme)*zv(lmo) + zv(lme)*zu(lmo))
           gvv(l) = gvv(l) + cp5*(rv(lme)*rv(lme) + zv(lme)*zv(lme))
     >               + phipog(l)*(rv(lmo)*rv(lmo) + zv(lmo)*zv(lmo))
     >               + shalf(l)*(rv(lme)*rv(lmo) + zv(lme)*zv(lmo))
     >                      + cp5*r1(lme)*r1(lme) + r1(lmo)*r1(lmo)
     >                      * phipog(l)  + shalf(l)*r1(lme)*r1(lmo)
 20     continue
*************
*                 PUT LAMBDA DERIVATIVES ON RADIAL HALF-MESH
*************
CDIR$ IVDEP
        do 30 l = nrzt,2,-1
        phipog(l) = phip(l)/gsqrt(l)
        lu(l,0) = cp5*phipog(l)*(lu(l,0)+lu(l-1,0)
     >          +     shalf(l)*(lu(l,1)+lu(l-1,1)))
 30     lv(l,0) = cp5*phipog(l)*(lv(l,0)+lv(l-1,0)
     >          +     shalf(l)*(lv(l,1)+lv(l-1,1)))
*************
*                 COMPUTE IOTA PROFILE
*************
        call getiota(phipog,guu,guv,lu,lv,wint,iotas,jv,
     >  czero,ns,ncurr)
*************
*                 PUT LAMBDA FORCES ON RADIAL HALF-MESH
*************
        do 40 l = 1,nrzt
        bsubu(l,0) = guu(l)*lv(l,0) + guv(l)*lu(l,0)
        bsubv(l,0) = guv(l)*lv(l,0) + gvv(l)*lu(l,0)
        bsubu(l,1) = shalf(l)*bsubu(l,0)
 40     bsubv(l,1) = shalf(l)*bsubv(l,0)
*************
*                 COMPUTE SUM OF KINETIC AND MAGNETIC PRESSURES
*                 AND RETURN IF FINAL OUTPUT LOOP (IEQUI=1)
*                 ON ENTRY, BSQ IS THE KINETIC PRESSURE
*************
        wb = -wp
        do 50 l = 2,nrzt
        bsq(l) = bsq(l) + cp5*(lv(l,0)*bsubu(l,0) + lu(l,0)*bsubv(l,0))
        phipog(l) = phipog(l)*wint(l)
 50     wb = wb + hs*dnorm*wint(l)*abs(gsqrt(l))*bsq(l)
        if(iequi.eq.1)return
*************
*                 AVERAGE LAMBDA FORCES ONTO FULL MESH
*                 NOTE: EDGE FORCE IS DOWN BY .5 UNTIL 90 LOOP
*************
        do 60 m = meven,modd
CDIR$ IVDEP
        do 70 l=1,nrzt
          bsubu(l,m) = cp5*(bsubu(l,m) + bsubu(l+1,m))
 70       bsubv(l,m) = cp5*(bsubv(l,m) + bsubv(l+1,m))
*do 70 l=1,nrzt-1
*       bsubu(l,m) = cp5*(bsubu(l,m) + bsubu(l+1,m))
*       bsubv(l,m) = cp5*(bsubv(l,m) + bsubv(l+1,m))
*70     continue
*        l     = nrzt
*       bsubu(l,m) = cp5*(bsubu(l,m)               )
*       bsubv(l,m) = cp5*(bsubv(l,m)               )
        do 90 l = ns,nrzt,ns
        bsubu(l,m) = c2p0*bsubu(l,m)
 90     bsubv(l,m) = c2p0*bsubv(l,m)
 60     continue
*************
*                 COMPUTE R,Z AND LAMBDA PRE-CONDITIONING MATRIX
*                 ELEMENTS AND FORCE NORMS EVERY NS4 STEPS
*************
        if(mod(iter2-iter1,ns4).eq.0)then
        call lamcal(phipog,guu,guv,gvv)
        call precondn(lu,bsq,gsqrt,r12,zs,zu12,zu,zu(1+nrzt),
     >  z1(1+nrzt),shalf,wint,pres,arm,ard,brm,brd,cr,ohs,
     >  cp25,cp5,c1p0,c1p5,czero,ns,iter2)
        call precondn(lu,bsq,gsqrt,r12,rs,ru12,ru,ru(1+nrzt),
     >  r1(1+nrzt),shalf,wint,pres,azm,azd,bzm,bzd,cr,ohs,
     >  cp25,cp5,c1p0,c1p5,czero,ns,iter2)
        do 110 l=2,nrzt
 110    guu(l) = guu(l)*r12(l)**2
        volume = hs*ssum(ns-1,vp(2),1)
        fnorm = dnorm/(sdot(nrzt,guu,1,wint,1)*(wb/volume)**2)
        fnorm1 = c1p0/sdot(4*mns-ns,xc(ns+1),1,xc(ns+1),1)
*************
*                 COMPUTE CONSTRAINT FORCE SCALING FACTOR (TCON)
*************
        do 120 js=2,ns-1
        ar(js) = czero
 120    az(js) = czero
        do 130 js=2,ns-1
        do 130 lk = 1,nrzt,ns
        ar(js) = ar(js) + wint(js+lk-1)*ru0(js+lk-1)**2
 130    az(js) = az(js) + wint(js+lk-1)*zu0(js+lk-1)**2
        do 140 js=2,ns-1
 140    tcon(js) = min(abs(ard(js,1)/ar(js)),abs(azd(js,1)/az(js)))
        tcon(ns) = cp5*tcon(ns-1)
        endif
*************
*                 STORE LU * LV COMBINATIONS USED IN FORCES
*************
        do 150 l=2,nrzt
           guu(l) = lv(l,0)*lv(l,0)*gsqrt(l)
           guv(l) = lv(l,0)*lu(l,0)*gsqrt(l)
           gvv(l) = lu(l,0)*lu(l,0)*gsqrt(l)
           lv(l,0)  = bsq(l)*gsqrt(l)/r12(l)
           lu(l,0)  = bsq(l)*r12(l)
 150    continue
        return
        end

        subroutine getiota(phipog,guu,guv,lu,lv,wint,iotas,jv,
     >  czero,ns,ncurr)
      include 'name1'
        real guu(ns,*),guv(ns,*),lu(ns,*),lv(ns,*),
     >  wint(ns,*),phipog(ns,*),iotas(*),jv(*)

        if(ncurr .eq. 0)goto 30
        do 10 js = 2,ns
        top = jv(js)
        bot = czero
        do 20 lk = 1,nznt
        top = top - wint(js,lk)*(guu(js,lk)*lv(js,lk)
     >      +                    guv(js,lk)*lu(js,lk))
 20     bot = bot + wint(js,lk)*phipog(js,lk)*guu(js,lk)
 10     iotas(js) = top/bot
*************
*                 ADD IOTA TO LAMBDA (ON HALF MESH NOW)
*************
 30     do 40 js = 2,ns
        do 40 lk = 1,nznt
 40     lv(js,lk) = lv(js,lk) + phipog(js,lk)*iotas(js)
        return
        end

        subroutine lamcal(phipog,guu,guv,gvv)
      include 'name1'
      include 'name0'
      include 'name2'
        real phipog(*),guu(*),guv(*),gvv(*),
     >  blam(nsd1),clam(nsd1),dlam(nsd1)
        data blam,clam,dlam/nsd1*0.,nsd1*0.,nsd1*0./

        do 10 js = 2,ns
        blam(js) = sdot(nznt,guu(js),ns,phipog(js),ns)
        dlam(js) = sdot(nznt,guv(js),ns,phipog(js),ns)*c2p0*real(nfp)
 10     clam(js) = sdot(nznt,gvv(js),ns,phipog(js),ns)
        do 20 m = 0,mpol1
        do 20 n = 0,nmax
        lmn = ns*(n + nmax1*m)
        tnn = real( (n*nfp)**2 )
        if( m.eq.0 .and. n.eq.0 )tnn = -c1p0
        do 30 js = jlam(m),ns
 30     faclam(js+lmn) = -c2p0*c2p0/( (blam(js)+blam(js+1))*tnn
     >  + sign((dlam(js)+dlam(js+1)),blam(js))*real(m*n)
     >  + (clam(js) + clam(js+1))*real(m*m))
 20     faclam(ns+lmn) = cp5*faclam(ns+lmn)
CDIR$ IVDEP
        do 40 l = 1,mns
 40     faclam(l+mns) = faclam(l)
        return
        end

        subroutine precondn(lu,bsq,gsqrt,r12,xs,xu12,xue,xuo,xodd,shalf,
     >  wint,pres,axm,axd,bxm,bxd,cx,ohs,cp25,cp5,c1p0,
     >  c1p5,czero,ns,iter2)
      include 'name1'
        real bsq(ns,*),gsqrt(ns,*),r12(ns,*),lu(ns,*),xu12(ns,*),
     >  xue(ns,*),xuo(ns,*),xodd(ns,*),xs(ns,*),wint(ns,*),shalf(*),
     >  pres(*),ax(nsd1,4),bx(nsd1,4),sm(nsd1),sp(0:nsd1),ptau(nznt),
     >  axm(nsd1,2),axd(nsd1,2),bxm(nsd1,2),bxd(nsd1,2),cx(*)
        save sm,sp

        if(iter2.le.1)then
        do 5 js=2,ns
        sm(js) = sqrt( (js - c1p5)/(js - c1p0) )
 5      sp(js) = sqrt( (js - cp5)/(js - c1p0) )
        sm(1) = czero
        sp(0) = czero
        sp(1) = sm(2)
        endif
*************
*                 COMPUTE PRECONDITIONING MATRIX ELEMENTS FOR R,Z
*                 FORCE (ALL ARE MULTIPLIED BY 0.5)
*************
        do 10 i = 1,4
        do 10 js = 1,ns+1
        ax(js,i) = czero
 10     bx(js,i) = czero
        do 15 js = 1,ns+1
 15     cx(js) = czero
        do 20 js = 2,ns
*************
*                 COMPUTE DOMINANT (1/DELTA-S)**2 PRECONDITIONING
*                 MATRIX ELEMENTS
*************
        do 30 lk = 1,nznt
        ptau(lk) = r12(js,lk)**2*(bsq(js,lk)-pres(js))
     >         * wint(js,lk)/gsqrt(js,lk)
        t1 = xu12(js,lk)*ohs
        t2 = cp25*(xue(js,lk)/shalf(js) + xuo(js,lk))/shalf(js)
        t3 = cp25*(xue(js-1,lk)/shalf(js) + xuo(js-1,lk))/shalf(js)
        ax(js,1) = ax(js,1) + ptau(lk)*t1*t1
        ax(js,2) = ax(js,2) + ptau(lk)*(-t1+t3)*(t1+t2)
        ax(js,3) = ax(js,3) + ptau(lk)*(t1+t2)**2
 30     ax(js,4) = ax(js,4) + ptau(lk)*(-t1+t3)**2
*************
*                 COMPUTE ORDER M**2 PRECONDITIONING MATRIX ELEMENTS
*************
        do 35 lk=1,nznt
        t1 = cp5*(xs(js,lk) + cp5*xodd(js,lk)/shalf(js))
        t2 = cp5*(xs(js,lk) + cp5*xodd(js-1,lk)/shalf(js))
        bx(js,1) = bx(js,1) + ptau(lk)*t1*t2
        bx(js,2) = bx(js,2) + ptau(lk)*t1**2
        bx(js,3) = bx(js,3) + ptau(lk)*t2**2
 35     cx(js) = cx(js) + cp25*lu(js,lk)**2*gsqrt(js,lk)*wint(js,lk)
 20     continue
        do 50 js = 1,ns
        axm(js,1) =-ax(js,1)
        axd(js,1) = ax(js,1) + ax(js+1,1)
        axm(js,2) = ax(js,2) * sm(js) * sp(js-1)
        axd(js,2) = ax(js,3)*sm(js)**2 + ax(js+1,4)*sp(js)**2
        bxm(js,1) = bx(js,1)
        bxm(js,2) = bx(js,1) * sm(js) * sp(js-1)
        bxd(js,1) = bx(js,2) + bx(js+1,3)
        bxd(js,2) = bx(js,2)*sm(js)**2 + bx(js+1,3)*sp(js)**2
 50     cx(js) = cx(js) + cx(js+1)
        return
        end

        subroutine forces(rbsq,guu,guv,gvv,sqrts,shalf,ohs,cp25,cp5,
     >  czero,nrzt,ns,ivac)
      include 'name1'
      include 'name3'
      include 'name4'
        real rbsq(*),guu(*),guv(*),gvv(*),guus(nrztd),guvs(nrztd),
     >  shalf(*),sqrts(*),bsqr(nrztd),gvvs(nrztd)
        equivalence (czmn(1+nrztd),guvs),(crmn(1+nrztd),guus)
*************
*                 ON ENTRY, ARMN=ZU,BRMN=ZS,AZMN=RU,BZMN=RS,CZMN=R*BSQ
*                 IT IS ESSENTIAL THAT CRMN,CZMN AT j=1 ARE ZERO INITIAL
*                 IN LOOPS, L (L1) INDEX REPRESENTS EVEN (ODD) COMPONENT
*************
        do 5 l = 1,nrzt+1,ns
        czmn(l) = czero
 5      crmn(l) = czero
        do 10 l=1,nrzt+1
        guus(l)  = guu(l)*shalf(l)
        guvs(l)  = guv(l)*shalf(l)
        gvvs(l)  = gvv(l)*shalf(l)
        armn(l)  = ohs*armn(l)*czmn(l)
        azmn(l)  =-ohs*azmn(l)*czmn(l)
        brmn(l)  = brmn(l)*czmn(l)
        bzmn(l)  =-bzmn(l)*czmn(l)
 10     bsqr(l)  = czmn(l)/shalf(l)
CDIR$ IVDEP
        do 15 l = 2,nrzt+1
        l1 = l+nrzt
        armn(l1) = armn(l)*shalf(l)
        azmn(l1) = azmn(l)*shalf(l)
        brmn(l1) = brmn(l)*shalf(l)
 15     bzmn(l1) = bzmn(l)*shalf(l)
*************
*                 CONSTRUCT CYLINDRICAL FORCE KERNELS (L=EVEN,
*                 L1=ODD COMPONENTS)
*************
CDIR$ IVDEP
        do 20 l = 1,nrzt
        l1 = l+nrzt
        bsqr(l) = cp25*(bsqr(l) + bsqr(l+1))
        czmn(l) = cp25*(czmn(l) + czmn(l+1))
        guu(l)  = cp5*(guu(l)  +  guu(l+1))
        guus(l) = cp5*(guus(l) + guus(l+1))
        guv(l)  = cp5*(guv(l)  +  guv(l+1))
        guvs(l) = cp5*(guvs(l) + guvs(l+1))
        gvv(l)  = cp5*(gvv(l)  +  gvv(l+1))
        gvvs(l) = cp5*(gvvs(l) + gvvs(l+1))
        armn(l)  = armn(l+1) - armn(l) + cp5*(crmn(l) + crmn(l+1))
     >   - (gvv(l)*r1(l) + gvvs(l)*r1(l1))
        crmn(l) = cp5*(crmn(l)*shalf(l) + crmn(l+1)*shalf(l+1))
        azmn(l)  = azmn(l+1) - azmn(l)
        brmn(l)  = cp5*(brmn(l) + brmn(l+1)) + z1(l1)*bsqr(l)
     >  - (guu(l)*ru(l) +guus(l)*ru(l1) +guv(l)*rv(l) +guvs(l)*rv(l1))
 20     bzmn(l)  = cp5*(bzmn(l) + bzmn(l+1)) - r1(l1)*bsqr(l)
     >  - (guu(l)*zu(l) +guus(l)*zu(l1) +guv(l)*zv(l) +guvs(l)*zv(l1))
CDIR$ IVDEP
        do 30 l=1,nrzt
        l1 = l+nrzt
        s2 = sqrts(l)**2
        guus2    = guu(l)*s2
        guvs2    = guv(l)*s2
        gvvs2    = gvv(l)*s2
        armn(l1) = armn(l1+1) -armn(l1) -(zu(l1)*czmn(l) +zu(l)*bsqr(l))
     >           + crmn(l) - (gvvs(l)*r1(l) + gvvs2*r1(l1))
        azmn(l1) = azmn(l1+1) -azmn(l1) + ru(l1)*czmn(l) +ru(l)*bsqr(l)
        brmn(l1) = cp5*(brmn(l1) + brmn(l1+1)) + z1(l1)*czmn(l)
     >  - (guus(l)*ru(l) + guus2*ru(l1) + guvs(l)*rv(l) + guvs2*rv(l1))
        bzmn(l1) = cp5*(bzmn(l1) + bzmn(l1+1)) - r1(l1)*czmn(l)
     >  - (guus(l)*zu(l) + guus2*zu(l1) + guvs(l)*zv(l) + guvs2*zv(l1))
        crmn(l)  = guv(l)  * ru(l) + guvs(l) * ru(l1)
     >           + gvv(l)  * rv(l) + gvvs(l) * rv(l1)
        crmn(l1) = guvs(l) * ru(l) + guvs2 * ru(l1)
     >           + gvvs(l) * rv(l) + gvvs2 * rv(l1)
        czmn(l)  = guv(l)  * zu(l) + guvs(l) * zu(l1)
     >           + gvv(l)  * zv(l) + gvvs(l) * zv(l1)
 30     czmn(l1) = guvs(l) * zu(l) + guvs2 * zu(l1)
     >           + gvvs(l) * zv(l) + gvvs2 * zv(l1)
*************
*                 ASSIGN EDGE FORCES (JS = NS) FOR FREE BOUNDARY
*                 CALCULATION
*************
        if( ivac.ge.1 )then
        do 40 m = 0,1
        do 40 lk = 1,nznt
        l = ns*lk
        l1 = l + m*nrzt
        armn(l1) = armn(l1) + zu0(l)*rbsq(lk)
 40     azmn(l1) = azmn(l1) - ru0(l)*rbsq(lk)
        endif
*************
*                 COMPUTE CONSTRAINT FORCE KERNELS
*************
CDIR$ IVDEP
        do 50 l = 1,nrzt
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
 50     zcon(l1)= zcon(l) * sqrts(l)
        return
        end

        subroutine tomnsp(frcc,frss,fzcs,fzsc,flcs,flsc,armn,brmn,
     >  crmn,azmn,bzmn,czmn,blmn,clmn,arcon,azcon,work1,work2,work3)
      include 'name1'
      include 'name0'
      include 'name2'
        real frcc(ns,0:nmax,0:mpol1), frss(ns,0:nmax,0:mpol1),
     >  fzcs(ns,0:nmax,0:mpol1),fzsc(ns,0:nmax,0:mpol1),
     >  flcs(ns,0:nmax,0:mpol1),flsc(ns,0:nmax,0:mpol1),
     >  armn(ns*nzeta,ntheta2,0:1),brmn(ns*nzeta,ntheta2,0:1),
     >  crmn(ns*nzeta,ntheta2,0:1),azmn(ns*nzeta,ntheta2,0:1),
     >  bzmn(ns*nzeta,ntheta2,0:1),czmn(ns*nzeta,ntheta2,0:1),
     >  blmn(ns*nzeta,ntheta2,0:1),clmn(ns*nzeta,ntheta2,0:1),
     >  arcon(ns*nzeta,ntheta2,0:1),azcon(ns*nzeta,ntheta2,0:1),
     >  work3(ns,nzeta,12),work2(ns*nzeta,12),work1(ns*nzeta*12)

        jmax = ns
        if( ivac.lt.1 )jmax = ns-1
*************
*                 BEGIN INVERSE FOURIER TRANSFORM
*                 DO THETA (U) INTEGRATION FIRST
*************
        do 40 m = 0,mpol1
        mparity = mod(m,2)
        do 10 l = 1,12*ns*nzeta
 10     work1(l) = czero
        do 20 i = 1,ntheta2
        do 30 jk= 1,ns*nzeta
        temp1 = armn(jk,i,mparity) + xmpq(m,1)*arcon(jk,i,mparity)
        temp3 = azmn(jk,i,mparity) + xmpq(m,1)*azcon(jk,i,mparity)
        work2(jk,01) = work2(jk,01) + temp1               *cosmui(i,m)
     >                               + brmn(jk,i,mparity)*sinmum(i,m)
        work2(jk,02) = work2(jk,02) - crmn(jk,i,mparity)*cosmui(i,m)
        work2(jk,03) = work2(jk,03) + temp1               *sinmu (i,m)
     >                               + brmn(jk,i,mparity)*cosmumi(i,m)
        work2(jk,04) = work2(jk,04) - crmn(jk,i,mparity)*sinmu (i,m)
        work2(jk,05) = work2(jk,05) + temp3               *cosmui(i,m)
     >                               + bzmn(jk,i,mparity)*sinmum(i,m)
        work2(jk,06) = work2(jk,06) - czmn(jk,i,mparity)*cosmui(i,m)
        work2(jk,07) = work2(jk,07) + temp3               *sinmu (i,m)
     >                               + bzmn(jk,i,mparity)*cosmumi(i,m)
        work2(jk,08) = work2(jk,08) - czmn(jk,i,mparity)*sinmu (i,m)
        work2(jk,09) = work2(jk,09) + blmn(jk,i,mparity)*sinmum(i,m)
        work2(jk,10) = work2(jk,10) - clmn(jk,i,mparity)*cosmui(i,m)
        work2(jk,11) = work2(jk,11) + blmn(jk,i,mparity)*cosmumi(i,m)
 30     work2(jk,12) = work2(jk,12) - clmn(jk,i,mparity)*sinmu (i,m)
 20     continue
*************
*                 NEXT, DO ZETA (V) INTEGRATION
*************
        do 45 n = 0,nmax
        do 45 k = 1,nzeta
        do 50 js= jmin2(m),jmax
        frcc(js,n,m) = frcc(js,n,m) + work3(js,k,01)*cosnv (k,n)
     >                              + work3(js,k,02)*sinnvn(k,n)
        frss(js,n,m) = frss(js,n,m) + work3(js,k,03)*sinnv (k,n)
     >                              + work3(js,k,04)*cosnvn(k,n)
        fzcs(js,n,m) = fzcs(js,n,m) + work3(js,k,05)*sinnv (k,n)
     >                              + work3(js,k,06)*cosnvn(k,n)
 50     fzsc(js,n,m) = fzsc(js,n,m) + work3(js,k,07)*cosnv (k,n)
     >                              + work3(js,k,08)*sinnvn(k,n)
        do 55 js= jlam(m),ns
        flcs(js,n,m) = flcs(js,n,m) + work3(js,k,09)*sinnv (k,n)
     >                              + work3(js,k,10)*cosnvn(k,n)
 55     flsc(js,n,m) = flsc(js,n,m) + work3(js,k,11)*cosnv (k,n)
     >                              + work3(js,k,12)*sinnvn(k,n)
 45     continue
 40     continue
        return
        end

        subroutine residue(gcr,gcz,gcl,bz0,work)
      include 'name1'
      include 'name0'
      include 'name2'
        common/precond/ard(nsd1,2),arm(nsd1,2),brd(nsd1,2),brm(nsd1,2),
     >        cr(nsd1),azd(nsd1,2),azm(nsd1,2),bzd(nsd1,2),bzm(nsd1,2)
        real gcr(ns,0:nmax,0:mpol1,2),gcz(ns,0:nmax,0:mpol1,2),
     >  gcl(ns,0:nmax,0:mpol1,2),work(mns,6)
*************
*                 IMPOSE M=1 MODE CONSTRAINT (ON Z CS, NTYPE = 1)
*************
        if(fsqz.lt.1.e-8.and.(iter1.ne.1))then
        do 10 n = 1,nmax
        do 10 js= 2,ns
 10     gcz(js,n,1,1) = czero
        endif
*************
*                 CONSTRUCT INVARIANT RESIDUALS
*************
        bz0 = c2p0*hs/bz0**2
        fsql = bz0*sdot(2*mns,gcl,1,gcl,1)
        call getfsq(gcr,gcz,fsqr,fsqz,fnorm,meven)
        fedge = fnorm*(sdot(mnd2,gcr(ns,0,0,1),ns,gcr(ns,0,0,1),ns)
     >        +        sdot(mnd2,gcz(ns,0,0,1),ns,gcz(ns,0,0,1),ns))
*************
*                 PERFORM PRECONDITIONING AND COMPUTE RESIDUES
*************
        call scalfor(gcr,arm,brm,ard,brd,cr,work,work(1,2),
     >  work(1,3),work(1,4),work(1,5))
        call scalfor(gcz,azm,bzm,azd,bzd,cr,work,work(1,2),
     >  work(1,3),work(1,4),work(1,5))
*************
*                 REDUCE R-Z FORCES IF INITIALLY VERY LARGE
*************
        fac = c1p0/(c1p0 + (fsqr+fsqz))
        if(iter1.gt.1)call sscal(4*mns,fac,gc,1)
        call getfsq(gcr,gcz,fsqr1,fsqz1,fnorm1,modd)
        if(iter1.eq.1)call sscal(4*mns,fac,gc,1)
*************
*                 CONSTRUCT PRECONDITIONED (SCALED) LAMBDA FORCES
*************
        do 30 l=1,2*mns
 30     gc(l+4*mns)=faclam(l)*gc(l+4*mns)
        fsql1 = hs*sdot(2*mns,gcl,1,gcl,1)
        return
        end

        subroutine scalfor(gcx,axm,bxm,axd,bxd,cx,bx,dx,ax,gm,alf)
      include 'name1'
      include 'name0'
      include 'name2'
        real gcx(ns,0:nmax,0:mpol1,2),ax(ns,0:nmax,0:mpol1),gm(*),
     >  bx(ns,0:nmax,0:mpol1),dx(ns,0:nmax,0:mpol1),alf(*),
     >  axm(nsd1,2),axd(nsd1,2),cx(*),bxm(nsd1,2),bxd(nsd1,2)

        jmax = ns
        if( ivac.lt.1 )jmax = ns-1
        do 10 m = 0,mpol1
        mp = mod(m,2)+1
        do 10 n = 0,nmax
        do 20 js=jmin2(m),jmax
        ax(js,n,m) = axm(js+1,mp) + bxm(js+1,mp)*m**2
        bx(js,n,m) = axm(js  ,mp) + bxm(js  ,mp)*m**2
 20     dx(js,n,m) = axd(js  ,mp) + bxd(js  ,mp)*m**2 +cx(js)*(n*nfp)**2
        if(m.eq.1)dx(2,n,1) = dx(2,n,1) + bx(2,n,1)
        dx(ns,n,m) = c1p1*dx(ns,n,m)
        bx(ns,n,m) = cp9*bx(ns,n,m)
 10     continue
        do 30 ntype=1,2
        call trid(ax,dx,bx,gcx(1,0,0,ntype),gm,alf,jmin2,jmax,ns)
 30     continue
        return
        end

        subroutine trid(a,d,b,c,gam,alf,jmin,nn,ns)
      include 'name1'
*************
*                 SOLVES B(I)*X(I-1)+D(I)*X(I)+A(I)*X(I+1)=C(I), I=IN,NN
*                 AND RETURNS ANSWER IN C(I)
*************
        real a(ns,0:mnd1),b(ns,0:mnd1),c(ns,0:mnd1),d(ns,0:mnd1),
     >  alf(0:mnd1,ns),gam(0:mnd1,ns)
        integer jmin(0:3)
        common /bounds/ mupper(3),mlower(3)
*************
*                 SEPARATE M=0 (IMODES=1), M=1 (IMODES=2), M>1 (IMODES=3
*************
        do 100 imodes = 1,3
        in = jmin(imodes-1)
        in1 = in + 1
        do 10 mn = mlower(imodes),mupper(imodes)
 10     gam(mn,in) = d(in,mn)
        do 20 i=in1,nn
        do 20 mn = mlower(imodes),mupper(imodes)
        alf(mn,i-1) = a(i-1,mn)/gam(mn,i-1)
 20     gam(mn,i)   = d(i,mn) - b(i,mn)*alf(mn,i-1)
        do 30 mn = mlower(imodes),mupper(imodes)
 30     c(in,mn) = c(in,mn)/gam(mn,in)
        do 40 i=in1,nn
        do 40 mn = mlower(imodes),mupper(imodes)
 40     c(i,mn) = (c(i,mn) - b(i,mn)*c(i-1,mn))/gam(mn,i)
        n2 = nn + in
        do 50 i=in1,nn
        i1 = n2 -i
        do 50 mn = mlower(imodes),mupper(imodes)
 50     c(i1,mn) = c(i1,mn) - alf(mn,i1)*c(i1+1,mn)
 100    continue
        return
        end

        subroutine getfsq(gcr,gcz,gnormr,gnormz,gnorm,mprecon)
      include 'name1'
      include 'name0'
      include 'name2'
        real gcr(ns,0:mnd1,2),gcz(ns,0:mnd1,2)
        gnormr = czero
        gnormz = czero
        jsmax = (ns-1) + mprecon
        do 10 mn = 0,mnd1
        do 10 js = 1,jsmax
        gnormr = gnormr + gnorm*(gcr(js,mn,1)**2 + gcr(js,mn,2)**2)
 10     gnormz = gnormz + gnorm*(gcz(js,mn,1)**2 + gcz(js,mn,2)**2)
        return
        end

        subroutine interp(xnew,xold,scale,nsin)
      include 'name1'
      include 'name0'
      include 'name2'
        real xnew(ns,0:nmax,0:mpol1,6),xold(nsin,0:nmax,0:mpol1,6)
        real scale(ns,0:nmax,0:mpol1)
        hsold = c1p0/real(nsin-1)
*************
*                 INTERPOLATE R,Z AND LAMBDA ON FULL GRID
*************
        do 30 ntype = 1,6
        do 20 n = 0,nmax
 20     xold(1,n,1,ntype) = c2p0*xold(2,n,1,ntype) - xold(3,n,1,ntype)
        do 10 js=1,ns
        sj=(js-1)*hs
        js1 = 1 + int(sj/hsold + 1.e-10)
        js2 = min0(js1+1,nsin)
        s1=(js1-1)*hsold
        do 10 mn = 0,mnd1
        xint=(sj-s1)/hsold
 10     xnew(js,mn,0,ntype) = ((c1p0-xint)*xold(js1,mn,0,ntype)
     >                      + xint*xold(js2,mn,0,ntype))/scale(js,mn,0)
        do 40 n = 0,nmax
 40     xnew(1,n,1,ntype) = czero
 30     continue
        return
        end

        subroutine bss(r12,rs,zs,ru12,zu12,shalf,bsubs,lu,lv,
     >  br,bphi,bz,cp25,cp5,nrzt)
      include 'name1'
      include 'name3'
        real r12(*),rs(*),zs(*),ru12(*),zu12(*),shalf(*),
     >  bsubs(*),br(*),bphi(*),bz(*),lu(*),lv(*)
        do 10  l=2,nrzt
        lodd = l + nrzt
        rv12 =      cp5*(   rv(l)+   rv(l-1)
     >       +shalf(l)*(rv(lodd)+rv(lodd-1)))
        zv12 =      cp5*(   zv(l)+   zv(l-1)
     >       +shalf(l)*(zv(lodd)+zv(lodd-1)))
        gsu= rs(l)*ru12(l) + zs(l)*zu12(l)
     >        +cp25*(r1(lodd)*ru(lodd)+r1(lodd-1)*ru(lodd-1)
     >             +z1(lodd)*zu(lodd)+z1(lodd-1)*zu(lodd-1)
     >        +(    r1(lodd)*ru(l)   +r1(lodd-1)*ru(l-1)
     >        +     z1(lodd)*zu(l)   +z1(lodd-1)*zu(l-1))/shalf(l))
        gsv= rs(l)*rv12    + zs(l)*zv12
     >        +cp25*(r1(lodd)*rv(lodd)+r1(lodd-1)*rv(lodd-1)
     >             +z1(lodd)*zv(lodd)+z1(lodd-1)*zv(lodd-1)
     >        +(    r1(lodd)*rv(l)   +r1(lodd-1)*rv(l-1)
     >        +     z1(lodd)*zv(l)   +z1(lodd-1)*zv(l-1))/shalf(l))
        br(l)   = lv(l)*ru12(l) + lu(l)*rv12
        bphi(l) =                 lu(l)*r12(l)
        bz(l)   = lv(l)*zu12(l) + lu(l)*zv12
 10     bsubs(l)= lv(l)*gsu     + lu(l)*gsv
        return
        end

        subroutine wrout(bsq,gsqrt,bsubu,bsubv,bsubs,br,bphi,bz,lu,lv)
      include 'name1'
      include 'name0'
      include 'name2'
      include 'name3'
        real
     >  bmod(nznt),bsq(ns,*),gsqrt(ns,*),bsubu(ns,*),bsubv(ns,*),
     >  bsubs(ns,*),br(*),bphi(*),bz(*),phi(nsd),rmnc(mnmax),
     >  zmns(mnmax),lmns(mnmax),xm(mnmax),xn(mnmax),bmodmn(mnmax),
     >  bmodmn1(mnmax),lu(ns,*),lv(ns,*)
*************
*                 THIS SUBROUTINE CREATES THE TEXT FILE WOUT WHICH
*                 CONTAINS THE CYLINDRICAL COORDINATE SPECTRAL
*                 COEFFICIENTS RMNC,ZMNS,LMNS (FULL MESH)
*************
C       WRITE(8,701)GAM,NFP,NS,MPOL,NMAX,MNMAX,ITFSQ,NITER/100+1
C 701   FORMAT(F10.3,7I6)
        WRITE(8    )GAM,NFP,NS,MPOL,NMAX,MNMAX,ITFSQ,NITER/100+1
        do 20 js=1,ns
        call convert(rmnc,zmns,lmns,xm,xn,js,xc,xc(1+mns),
     >  xc(1+2*mns),xc(1+3*mns),xc(1+4*mns),xc(1+5*mns))
        mn = 0
        do 25 lk=1,nznt
 25     bmod(lk)=sqrt(c2p0*abs(bsq(js,lk)/bscale**2-pres(js)))
        do 20 m = 0,mpol1
        nmin0 = -nmax
        if(m.eq.0)nmin0 = 0
        do 20 n = nmin0,nmax
        mn = mn+1
        dmult = c2p0*dnorm/(mscale(m)*nscale(iabs(n)))
        if(m.eq.0.and.n.eq.0)dmult=cp5*dmult
        gmn = czero
        bmn = czero
        bsubumn = czero
        bsubvmn = czero
        bsubsmn = czero
        bsupumn = czero
        bsupvmn = czero
        if(js.eq.1)goto 22
        do 200 j = 1,ntheta2
        do 200 k = 1,nzeta
        lk = k + nzeta*(j-1)
        if(n.ge.0)then
        tcosi = dmult*(cosmui(j,m)*cosnv(k,n)+sinmu(j,m) *sinnv(k,n))
        tsini = dmult*(sinmu (j,m)*cosnv(k,n)-cosmui(j,m)*sinnv(k,n))
        else
        tcosi = dmult*(cosmui(j,m)*cosnv(k,-n)-sinmu(j,m) *sinnv(k,-n))
        tsini = dmult*(sinmu (j,m)*cosnv(k,-n)+cosmui(j,m)*sinnv(k,-n))
        endif
        bmn = bmn + tcosi*bmod(lk)
        gmn = gmn + tcosi*gsqrt(js,lk)
        bsubumn = bsubumn + tcosi*bsubu(js,lk)
        bsubvmn = bsubvmn + tcosi*bsubv(js,lk)
        bsubsmn = bsubsmn + tsini*bsubs(js,lk)
        bsupumn = bsupumn + tcosi*lv(js,lk)
        bsupvmn = bsupvmn + tcosi*lu(js,lk)
 200    continue
        if(js.eq.ns/2)bmodmn(mn) = bmn
        if(js.eq.ns)bmodmn1(mn) = bmn
        goto 20
C22     WRITE(8,702)NINT(XM(MN)),NINT(XN(MN))
 22     WRITE(8    )NINT(XM(MN)),NINT(XN(MN))
C20     WRITE(8,703)RMNC(MN),ZMNS(MN),LMNS(MN),
 20     WRITE(8    )RMNC(MN),ZMNS(MN),LMNS(MN),
     >  bmn,gmn,bsubumn/bscale,bsubvmn/bscale,bsubsmn/bscale,
     >  bsupumn/bscale,bsupvmn/bscale
 702    format(2i10)
 703    format(5e20.13)
        phi(1) = czero
        do 30 js = 2,ns
 30     phi(js) = twopi*hs*ssum(js-1,phips(2),1)
        fac = abs(bscale)**(gam-c2p0)
C       WRITE(8,703)(IOTAS(JS),MASS(JS)*FAC,PRES(JS),PHIPS(JS),BUCO(JS),
        WRITE(8    )(IOTAS(JS),MASS(JS)*FAC,PRES(JS),PHIPS(JS),BUCO(JS),
     >  bvco(js),phi(js),vp(js),ju(js),jv(js),specw(js),js=2,ns)
        WRITE(8    )(FSQT(I),WDOT(I),I=1,100)
C       WRITE(8,703)(FSQT(I),WDOT(I),I=1,100)
        if(nvac.eq.0)return
        write(3,60)rbtor/bscale,ctor/bscale,bscale
 60     format(/,' 2*PI*R*BTOR = ',1pe16.8,' NET TOROIDAL CURRENT = ',
     >  1pe16.8,' BSCALE = ',1pe16.8)
        do 80 iprint=1,2
        if(iprint.eq.1)write(3,70)
        if(iprint.eq.2)write(3,75)
        ntskip=1+ntheta1/12
        nzskip=1+nzeta/6
        do 80 l=1,nzeta,nzskip
        zeta = real(360*(l-1))/real(nzeta)
        do 80 k=1,ntheta2,ntskip
        lk=l+nzeta*(k-1)
        nl = ns*lk
        if(iprint.eq.1)write(3,90)zeta,r1(nl)+r1(nl+nrzt),
     >  z1(nl)+z1(nl+nrzt),(bsqsav(lk,n)/bscale**2,n=1,3),
     >  bsqvac(lk)/bscale**2
        if(iprint.eq.2)write(3,95)zeta,r1(nl)+r1(nl+nrzt),
     >  z1(nl)+z1(nl+nrzt),(1.5*br(nl) - 0.5*br(nl-1))/bscale,
     >  (1.5*bphi(nl) - 0.5*bphi(nl-1))/bscale,
     >  (1.5*bz(nl) - 0.5*bz(nl-1))/bscale,brv(lk)/bscale,
     >  bphiv(lk)/bscale,bzv(lk)/bscale
 80     continue
 70     format(/,4x,'ZETA',8x,' Rb ',8x,' Zb ',6x,
     >  'BSQMHDI',5x,'BSQVACI',5x,'BSQMHDF',5x,'BSQVACF',/)
 75     format(/,4x,'ZETA',8x,' Rb ',8x,' Zb ',6x,
     >  'BR',8x,'BPHI',6x,'BZ',8x,'BRv',7x,'BPHIv',5x,'BZv',/)
 90     format(1pe10.2,1p6e12.4)
 95     format(1pe10.2,1p2e12.4,1p6e10.2)
        write(3,100)
 100    format(//,3x,'mb',2x,'nb',9x,'rbc',9x,'zbs',3x,'|B|(s=.5)',
     >  3x,'|B|(s=1.)',6x,'mb',2x,'nb',9x,'rbc',9x,'zbs',3x,
     >  '|B|(s=.5)',3x,'|B|(s=1.)'/)
        do 110 mn=1,mnmax,2
 110    write(3,115)nint(xm(mn)),nint(xn(mn)/nfp),rmnc(mn),zmns(mn),
     >  bmodmn(mn),bmodmn1(mn),nint(xm(mn+1)),nint(xn(mn+1)/nfp),
     >  rmnc(mn+1),zmns(mn+1),bmodmn(mn+1),bmodmn1(mn+1)
 115    format(i5,i4,1p4e12.4,3x,i5,i4,1p4e12.4)
        write(3,120)
 120    format(/,3x,'mf',2x,'nf',5x,'potvacs',6x,'mf',2x,'nf',5x,
     >  'potvacs',6x,'mf',2x,'nf',5x,'potvacs'/)
        do 130 mn=1,mpmax,3
 130    write(3,135)nint(xmpot(mn)),nint(xnpot(mn)/nfp),
     >  potvac(mn)/bscale,nint(xmpot(mn+1)),nint(xnpot(mn+1)/nfp),
     >  potvac(mn+1)/bscale,nint(xmpot(mn+2)),nint(xnpot(mn+2)/nfp),
     >  potvac(mn+2)/bscale
 135    format(i5,i4,1pe12.4,3x,i5,i4,1pe12.4,3x,i5,i4,1pe12.4)
        return
        end

        subroutine convert(rmnc,zmns,lmns,xm,xn,js,rmncc,rmnss,
     >  zmncs,zmnsc,lmncs,lmnsc)
      include 'name1'
      include 'name0'
      include 'name2'
        real rmnc(*),zmns(*),lmns(*),xm(*),xn(*),
     >  rmncc(ns,0:nmax,0:mpol1),rmnss(ns,0:nmax,0:mpol1),
     >  zmncs(ns,0:nmax,0:mpol1),zmnsc(ns,0:nmax,0:mpol1),
     >  lmncs(ns,0:nmax,0:mpol1),lmnsc(ns,0:nmax,0:mpol1)
*************
*                 CONVERTS INTERNAL MODE REPRESENTATION TO STANDARD
*                 FORM FOR OUTPUT (COEFFICIENTS OF cos(mu-nv), sin(mu-nv
*************
        mn = 0
        do 10 m = 0,mpol1
        nmin0 = -nmax
        if(m.eq.0)nmin0 = 0
        do 10 n = nmin0,nmax
        n1 = iabs(n)
        t1 = mscale(m)*nscale(n1)
        mn = mn + 1
        xm(mn) = real(m)
        xn(mn) = real(n*nfp)
        if( m.eq.0 )then
                rmnc(mn) = t1*rmncc(js,n,m)
                zmns(mn) =-t1*zmncs(js,n,m)
                lmns(mn) =-t1*lmncs(js,n,m)
        else if( n.eq.0 )then
                rmnc(mn) = t1*rmncc(js,n,m)
                zmns(mn) = t1*zmnsc(js,n,m)
                lmns(mn) = t1*lmnsc(js,n,m)
        else if( js.gt.1 )then
                sign = real(n/n1)
                rmnc(mn) = cp5*t1*(rmncc(js,n1,m) + sign*rmnss(js,n1,m))
                zmns(mn) = cp5*t1*(zmnsc(js,n1,m) - sign*zmncs(js,n1,m))
                lmns(mn) = cp5*t1*(lmnsc(js,n1,m) - sign*lmncs(js,n1,m))
        else if( js.eq.1 )then
                rmnc(mn) = czero
                zmns(mn) = czero
                lmns(mn) = czero
        endif
 10     continue
        return
        end

        subroutine eqfor(bu,bv,bsq,rmag,zmag)
      include 'name1'
      include 'name0'
      include 'name2'
        real bu(ns,1),bv(ns,1),bsq(ns,1),rmag(ns,0:nmax),zmag(ns,0:nmax)

        write(3,5)
  5     format(/,'     S       EQUIF      PHI       CHI     <JTHETA>  ',
     >  '   IOTA     <JZETA>   <BSUBV>    P''/V''    V''/PHIP   <M>',/)
        betaxis =
     >   c1p5*pres(2)/(dnorm*sdot(nznt,bsq(2,1),ns,wint(2),ns)-pres(2))
     >  - cp5*pres(3)/(dnorm*sdot(nznt,bsq(3,1),ns,wint(2),ns)-pres(3))
        phi1 = czero
        chi1 = czero
        do 10 i = 2,ns
        buco(i) = dnorm*sdot(nznt,bu(i,1),ns,wint(2),ns)/bscale
        bvco(i) = dnorm*sdot(nznt,bv(i,1),ns,wint(2),ns)/bscale
        pres(i) = pres(i)/bscale**2
 10     phips(i) = phips(i)/bscale
        do 20 js = 2,ns-1
        t0 = cp5*(vp(js+1)+vp(js))
        aiotaf = cp5*(iotas(js) + iotas(js+1))
        jv(js) =  isigng*ohs*(buco(js+1)-buco(js))/t0
        ju(js) = -isigng*ohs*(bvco(js+1)-bvco(js))/t0
        t1 = jv(js)*aiotaf
        t2 = ohs*(pres(js+1)-pres(js))/t0
        t3 = cp5*(vp(js+1)/phips(js+1) + vp(js)/phips(js))
        phi1 = phi1 + phips(js)/ohs
        chi1 = chi1 + phips(js)*iotas(js)/ohs
        equif(js) = (ju(js)-t1-t2*t3)/(abs(t1)+abs(ju(js))+abs(t2*t3))
        es = (js-1)/real(ns-1)
 20     write(3,30)es,equif(js),isigng*twopi*phi1,isigng*twopi*chi1,
     >  ju(js),aiotaf,jv(js),bvco(js),t2,t3,specw(js)
*30     FORMAT(1P9E10.2,1P1E11.3,0PF7.3)
 30     FORMAT(1P5E10.2,1p1E11.3,1P3E10.2,1P1E11.3,0PF7.3)
        volf = (twopi**2)*ssum(ns,vp,1)/real(ns-1)
        write(3,40)voli,volf,betaxis
 40     format(/,'  INITIAL VOLUME  =  ',1pe20.10,'  FINAL VOLUME  =  ',
     >  1pe20.10,/'  BETA ON AXIS    =  ',2x,1pe12.4,/)
        write(3,50)
 50     format(2x,'MAGNETIC AXIS COEFFICIENTS'/,
     >  '    n     rmag       zmag',/)
        do 60 n = 0,nmax
        t1 = mscale(0)*nscale(n)
 60     write(3,70)n,t1*rmag(1,n),-t1*zmag(1,n)
 70     format(i5,1p2e12.4)
        return
        end

        subroutine spectrum(rcc,rss,zcs,zsc)
      include 'name1'
      include 'name0'
      include 'name2'
        real rcc(ns,0:nmax,0:mpol1),rss(ns,0:nmax,0:mpol1),
     >  zsc(ns,0:nmax,0:mpol1),zcs(ns,0:nmax,0:mpol1)
        real t1(nsd),dnumer(nsd),denom(nsd)
        do 5 js=2,ns
        dnumer(js) = czero
 5      denom (js) = czero
        do 10 n = 0,nmax
        do 10 m = 1,mpol1
        scale = (mscale(m)*nscale(n))**2
        do 15 js=2,ns
 15     t1(js) = (rcc(js,n,m)**2 + rss(js,n,m)**2
     >          + zcs(js,n,m)**2 + zsc(js,n,m)**2)*scale
        do 25 js=2,ns
        dnumer(js) = dnumer(js) + t1(js)*xmpq(m,3)
 25     denom (js) = denom (js) + t1(js)*xmpq(m,2)
 10     continue
        do 30 js=2,ns
 30     specw(js) = dnumer(js)/denom(js)
        return
        end

        subroutine printout(i,w0,r00)
      include 'name1'
      include 'name0'
      include 'name2'
        betav = wp/wb
        w = w0*twopi*twopi/nfp
        avm = czero
        den = czero
        specw(1) = c1p0
        call spectrum(xc,xc(1+mns),xc(1+2*mns),xc(1+3*mns))
        do 10 j=2,ns
        den = den+phips(j)
 10     avm = avm+phips(j)*(specw(j)+specw(j-1))
        avm = cp5*avm/den
        if( ivac.ge.1 )delbsq=sdot(nznt,dbsq,1,wint(2),ns)/
     >                 sdot(nznt,bsqsav(1,3),1,wint(2),ns)
c       if(i.eq.1.and.nvac.ne.0)print 20
        if(i.eq.1.and.nvac.ne.0)write(3,15)
c       if(i.eq.1.and.nvac.eq.0)print 30
        if(i.eq.1.and.nvac.eq.0)write(3,25)
 15     format(/,' ITER    FSQR      FSQZ      FSQL      fsqr      ',
     >  'fsqz      DELT     R00(0)        WMHD        BETA     ',
     >  '<M>   DEL-BSQ   FEDGE',/)
 20     format(/,' ITER    FSQR      FSQZ      FSQL      ',
     >  'R00(0)     WMHD      DEL-BSQ',/)
 25     format(/,' ITER    FSQR      FSQZ      FSQL      fsqr      ',
     >  'fsqz      DELT     R00(0)        WMHD        BETA     <M>',/)
 30     format(/,' ITER    FSQR      FSQZ      FSQL      fsqr      ',
     >  'fsqz     R00(0)     WMHD',/)
        if(nvac.ne.0)goto 35
*       PRINT 45, I,FSQR,FSQZ,FSQL,FSQR1,FSQZ1,R00,W
        write(3,40)i,fsqr,fsqz,fsql,fsqr1,fsqz1,delt,r00,w,betav,avm
        return
*35     PRINT 50, I,FSQR,FSQZ,FSQL,R00,W,DELBSQ
 35     CONTINUE
        write(3,40)i,fsqr,fsqz,fsql,fsqr1,fsqz1,delt,r00,w,
     >  betav,avm,delbsq,fedge
 40     format(i5,1p6e10.2,1pe11.4,1pe15.8,1pe9.2,0pf7.3,1p2e9.2)
 45     format(i5,1p5e10.2,1pe10.3,1p2e11.4)
 50     format(i5,1p3e10.2,1pe10.3,1p2e11.4)
        return
        end
