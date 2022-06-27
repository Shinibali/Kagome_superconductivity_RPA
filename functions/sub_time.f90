SUBROUTINE time()
 IMPLICIT NONE

    character(8)  :: date
    character(10)  :: hhmmss
    character(5)  :: zone
    integer,dimension(8) :: values
    ! using keyword arguments
    call date_and_time(date,hhmmss,zone,values)
    write(*,*) 'Time:',values(5),'-',values(6),'-',values(7)

RETURN
END subroutine
