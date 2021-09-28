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
