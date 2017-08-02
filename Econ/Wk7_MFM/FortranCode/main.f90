program main

    use ifport
    
    use vars
    use solve
    implicit none

    integer :: ec
	real    :: time_0, time_T
	real(8) :: dist(10000)
	character(24) :: start_time

    ! Program body -------------------------------------------------------------
    start_time = fdate()
	call load_params()
    call initialize()

    ! Compute model ------------------------------------------------------------
    call cpu_time(time_0)
    print "(a24/,a,i0,a,a)", fdate(), "Compute ", T, " periods w ", param_str

    call compute_prep()
	call save_results()

    dist(1) = 0d0
    periods = 0
    do while (periods < T)
        call compute(1, dist(periods+1))
		call save_results()
        periods = periods + 1
    end do
    call compute_free()

    call cpu_time(time_T)

    ! Simulate model -----------------------------------------------------------
    if (SimLen > 0) then
        print "(a)", param_str
		print "(a,i0,a\)", "Simulate ", SimLen, " periods"

        call simulate(0.5d0)
        call histogram("sim-hist.txt")
    end if

    call cpu_time(time_T)
    print "(a,a24,a,f0.3,a)", "", fdate(), ": Total time = ", (time_T-time_0)/60, " min"

    !ec = system("D:/Programs/gnuplot/bin/gnuplot.exe D:/Dropbox/macrofin-modeling/code-bs/plots_ps.txt")
    
    ! Free memory --------------------------------------------------------------    
    call deinitialize()

end