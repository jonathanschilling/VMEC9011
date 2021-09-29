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
