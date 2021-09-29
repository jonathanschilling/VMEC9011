
        subroutine output
      include 'name1'
      include 'name0'
      include 'name2'

        iequi = 1
        call funct3d
*       PRINT 10, IJACOB
        write(3,10)ijacob
 10     format(/,'  NUMBER OF JACOBIAN RESETS = ',i4)
*       PRINT 20, TIMER(0),TIMER(1)
        write(3,20)timer(0),timer(1)
 20     format (/,
     >  '  TOTAL COMPUTATIONAL TIME :       ',1pe10.2,' SECONDS',/
     >  '  TIME IN VACUUM LOOP :            ',1pe10.2,' SECONDS',/)
        return
        end
