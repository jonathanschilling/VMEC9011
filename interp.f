
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
