        subroutine getfsq(gcr,gcz,gnormr,gnormz,gnorm,mprecon)
      include 'name1'
      include 'name0'
      include 'name2'
        real gcr(ns,0:mnd1,2),gcz(ns,0:mnd1,2)
        gnormr = czero
        gnormz = czero
        jsmax = (ns-1) + mprecon
        do 10 mn = 0,mnd1
        do 10 js = 1,jsmax
        gnormr = gnormr + gnorm*(gcr(js,mn,1)**2 + gcr(js,mn,2)**2)
 10     gnormz = gnormz + gnorm*(gcz(js,mn,1)**2 + gcz(js,mn,2)**2)
        return
        end
