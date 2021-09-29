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
