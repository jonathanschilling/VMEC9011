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
