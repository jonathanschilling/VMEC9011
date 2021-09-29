        subroutine vsetup
      include 'name1'
      include 'name0'
      include 'name2'
      include 'name3'

        do 10 l = 1,nrztd
        rcon0(l) = czero
 10     zcon0(l) = czero
        do 20 i = 0,10
 20     timer(i) = czero
        iequi = 0
        itfsq = 0
        bscale = c1p0
        delbsq = c1p0
        ivac   = 0
        print *, "vsetup done"
        return
        end
