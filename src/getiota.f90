        subroutine getiota(phipog,guu,guv,lu,lv,wint,iotas,jv,
     >  czero,ns,ncurr)
      include 'name1'
        real guu(ns,*),guv(ns,*),lu(ns,*),lv(ns,*),
     >  wint(ns,*),phipog(ns,*),iotas(*),jv(*)

        if(ncurr .eq. 0)goto 30
        do 10 js = 2,ns
        top = jv(js)
        bot = czero
        do 20 lk = 1,nznt
        top = top - wint(js,lk)*(guu(js,lk)*lv(js,lk)
     >      +                    guv(js,lk)*lu(js,lk))
 20     bot = bot + wint(js,lk)*phipog(js,lk)*guu(js,lk)
 10     iotas(js) = top/bot
*************
*                 ADD IOTA TO LAMBDA (ON HALF MESH NOW)
*************
 30     do 40 js = 2,ns
        do 40 lk = 1,nznt
 40     lv(js,lk) = lv(js,lk) + phipog(js,lk)*iotas(js)
        return
        end
