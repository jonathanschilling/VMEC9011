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
