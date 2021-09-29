
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
