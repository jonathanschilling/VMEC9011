
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
