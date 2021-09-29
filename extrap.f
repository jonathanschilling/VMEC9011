
        subroutine extrap(rmncc,rmnss,zmncs,zmnsc,lmncs,lmnsc,x3,x4,ns)
      include 'name1'
        real rmncc(ns,0:mnd1),rmnss(ns,0:mnd1),x3(0:mnd1),
     >       zmncs(ns,0:mnd1),zmnsc(ns,0:mnd1),x4(0:mnd1),
     >       lmncs(ns,0:mnd1),lmnsc(ns,0:mnd1)

        do 10 mn = 2*nmax1,mnd1
        rmncc(2,mn) = x3(mn)*rmncc(3,mn) + x4(mn)*rmncc(4,mn)
        rmnss(2,mn) = x3(mn)*rmnss(3,mn) + x4(mn)*rmnss(4,mn)
        zmncs(2,mn) = x3(mn)*zmncs(3,mn) + x4(mn)*zmncs(4,mn)
 10     zmnsc(2,mn) = x3(mn)*zmnsc(3,mn) + x4(mn)*zmnsc(4,mn)
        do 20 mn = nmax1,mnd1
        lmncs(2,mn) = x3(mn)*lmncs(3,mn) + x4(mn)*lmncs(4,mn)
 20     lmnsc(2,mn) = x3(mn)*lmnsc(3,mn) + x4(mn)*lmnsc(4,mn)
        do 30 n = 1,nmax
 30     lmncs(1,n) = x3(0)*lmncs(2,n) - lmncs(3,n)
        return
        end
