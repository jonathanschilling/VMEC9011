
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
