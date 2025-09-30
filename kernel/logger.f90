module logger
!
!   Simple timing logger
!
    use numtype
    implicit none

    public :: logger_start, logger_stop, get_total_time

    real(dp) :: total_time = 0.0_dp
    real(dp) :: t_start

contains


    subroutine logger_start()
    !
    ! start logging
    !
        call cpu_time(t_start)
    end subroutine logger_start


    subroutine logger_stop()
    !
    !   end logging
    !
        real(dp) :: t_end, elapsed
        call cpu_time(t_end)
        elapsed = t_end - t_start
        total_time = total_time + elapsed
        if (isDebug .eq. 1) print*, "Elapsed time (s): ", elapsed
        if (isDebug .eq. 1) print*, "Total accumulated time (s): ", total_time
    end subroutine logger_stop

    function get_total_time() result(t)
    !
    !   return total compute time
    !
        real(dp) :: t
        t = total_time
    end function get_total_time


    subroutine logger_reset()
    !
    !   reset total compute time
    !
        total_time = 0._dp

    end subroutine logger_reset

end module logger
