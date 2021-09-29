
        subroutine restart
      include 'name1'
      include 'name0'
      include 'name2'
        goto (10,20,20)irst
 10     call scopy(neqs,xc,1,xstore,1)
        return
 20     do 30 l = 1,neqs
        xcdot(l) = czero
 30     xc(l) = xstore(l)
        delt = delt*(cp96*(irst-2) + cp9*(3-irst))
        irst = 1
        return
        end
