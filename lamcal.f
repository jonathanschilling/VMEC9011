
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
