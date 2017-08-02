module vars

    implicit none
    real(8), parameter :: pi = 3.14159265359d0

    ! Model parameters (loaded from params.txt, see description in file)
    integer :: nu
    real(8) :: beta, betah, betaf, gammah, gammaf, kmin, dloss, dprob, interval, eta
    real(8) :: omega_min, omega_max, omega_ext
    real(8) :: mgammainvh, mgammainvf, mgammah, mgammaf, nuinv
    real(8) :: Ah, Af, delta, sd_z

    integer :: Ny = 3, Nz = 2
    integer :: M, Nw, periods, T
    integer :: burnin = 20, SimLen, seed = 1543264423

    integer :: dw                           ! Current exog state
    integer :: z, zw                        ! Current exog state
    real(8) :: omega, AD                    ! Current endo state
    real(8) :: prc(2)                       ! Parameters for funcs and subs
    real(8) :: prm(3)			            ! Parameters for funcs and subs

    real(8) :: tol = 1d-12					! Tolerance for eq. system
    real(8) :: eps = 1d-15					! "practical" zero

    integer :: wIter = 1000					! Maximum number of iterations
    real(8) :: tol_w = 1d-14				! Tolerance for omega transition eq.
    real(8) :: w_bnd(2)

    character(58) :: param_str              ! Parametrization

    real(8), pointer :: yGrid(:), wGrid(:), Zs(:), pz(:)
    real(8), pointer :: sol_t0(:,:,:), pol_t0(:,:,:), prc_t0(:,:,:), err_t0(:,:,:)
    real(8), pointer :: sol_t1(:,:,:), pol_t1(:,:,:), prc_t1(:,:,:)
    real(8), pointer :: trans(:,:,:), ushock(:)
    real(8), pointer :: brk_w(:), spl_qk(:,:,:), spl_wp(:,:,:,:), spl_kh(:,:,:)

contains

    subroutine initialize()
    use rnun_int
        integer :: i, ik, status
        real(8) :: wstep

        mgammainvh = -1d0/gammah
        mgammainvf = -1d0/gammaf
        mgammah = -gammah
        mgammaf = -gammaf

        nuinv = 1d0/nu

		M = (Nw-1)/2+1
        if (dprob > eps) then; Nz = Nz+1; end if

        allocate(Zs(Nz), pz(Nz), &
            yGrid(Ny), wGrid(Nw), &
            brk_w(Nw+4), spl_qk(4,Nw+4,Ny), spl_wp(4,Nw+4,Ny,Nz), spl_kh(4,Nw+4,Ny), &
			sol_t0(4,Nw,Ny), pol_t0(4,Nw,Ny), prc_t0(2,Nw,Ny), err_t0(4,Nw,Ny), &
			sol_t1(4,Nw,Ny), pol_t1(4,Nw,Ny), prc_t1(2,Nw,Ny), &
            trans(Nz,Nw,Ny), ushock(SimLen), &
			stat = status)

		call random_seed(seed)
        call random_number(ushock)

        ! Capital productivity process
        if (dprob > eps) then
            Zs = [-sd_z-dloss, -sd_z, sd_z]
            pz = [dprob, 0.5d0*(1d0-dprob), 0.5d0*(1d0-dprob)]
        else
            Zs = [-sd_z, sd_z]
            pz = [0.5d0, 0.5d0]
        end if

		interval = omega_max-omega_min
		wstep = interval/(Nw-1)
		wGrid = [ (dble(i)*wstep, i = 0, Nw-1) ] + omega_min       ! Create grid for omega
		!wGrid = [ (1d0-cos(pi*(dble(i)-0.5d0)/Nw), i = 0, Nw-1) ]*(omega_max-omega_min)/2 + omega_min       ! Create grid for omega
		w_bnd = [ omega_min-omega_ext, omega_max+omega_ext ]   ! Bounds for omega

        yGrid = [ Af+Zs(1),(Ah+Af)/2,Ah+Zs(Nz) ]

        do ik = 1, Ny
            pol_t0(1,:,ik) = wGrid                  ! Consumption share of H
            pol_t0(2,:,ik) = wGrid                  ! Capital owndership of H
            pol_t0(3,:,ik) = 1d0-wGrid              ! Capital owndership of F
            pol_t0(4,:,ik) = 0d0                    ! Bond investment of H

            sol_t0(1,:,ik) = log(wGrid)-log(1d0-wGrid) ! Consumption share of H
            sol_t0(2,:,ik) = wGrid**nuinv           ! Capital owndership of H
            sol_t0(3,:,ik) = max(1d0-wGrid-kmin,0d0)**nuinv  ! Capital owndership of F
            sol_t0(4,:,ik) = 0d0                    ! Bond investment of H

            prc_t0(1,:,ik) = 0d0
            prc_t0(2,:,ik) = (betah+betaf)/2d0
        end do
    end subroutine

    subroutine deinitialize()
        integer :: ec

        deallocate(wGrid, Zs, pz, &
			sol_t0, pol_t0, prc_t0, err_t0, &
			sol_t1, pol_t1, prc_t1, &
			trans, brk_w, spl_qk, spl_wp, spl_kh, stat = ec)
    end subroutine

	subroutine load_params()
		open(1, file="params.txt")
		read(1, "(f8.4)") betah
		read(1, "(f8.4)") betaf
		read(1, "(f8.4)") gammah
		read(1, "(f8.4)") gammaf
		read(1, "(f8.4)") Ah
		read(1, "(f8.4)") Af
		read(1, "(f8.4)") kmin
		read(1, "(f8.4)") dloss
		read(1, "(f8.4)") dprob
		read(1, "(f8.4)") sd_z
		read(1, "(f8.4)") omega_min
		read(1, "(f8.4)") omega_max
		read(1, "(f8.4)") omega_ext
		read(1, "(i8)") nu
		read(1, "(f8.4)") eta
		read(1, "(i8)") T
		read(1, "(i8)") Nw
		read(1, "(i8)") SimLen
		close(1)

		sd_z = sd_z/100d0
    end subroutine

    subroutine save_results()
        integer :: i

        open( 1, file = "results\sol-ch.txt")
        open( 2, file = "results\sol-kh.txt")
        open( 3, file = "results\sol-kf.txt")
        open( 4, file = "results\sol-bh.txt")
        open( 5, file = "results\sol-bf.txt")
        open( 6, file = "results\sol-qk.txt")
        open( 7, file = "results\sol-qb.txt")
        open( 8, file = "results\sol-wp.txt")
        do i = 1, Nw
            write( 1, "(f8.4,<Ny*Nz>f10.6)") wGrid(i), pol_t0(1, i, :)
            write( 2, "(f8.4,<Ny*Nz>f10.6)") wGrid(i), pol_t0(2, i, :)
            write( 3, "(f8.4,<Ny*Nz>f10.6)") wGrid(i), pol_t0(3, i, :)
            write( 4, "(f8.4,<Ny*Nz>f10.6)") wGrid(i), pol_t0(4, i, :)
            write( 5, "(f8.4,<Ny*Nz>f10.6)") wGrid(i), -pol_t0(4, i, :)
            write( 6, "(f8.4,<Ny*Nz>f10.6)") wGrid(i), prc_t0(1, i, :)
            write( 7, "(f8.4,<Ny*Nz>f10.6)") wGrid(i), prc_t0(2, i, :)
            write( 8, "(f8.4,<Ny*Nz*Nz>f10.6)") wGrid(i), trans (:, i, :)
		end do
        do i = 1, 9; close(i); end do
    end subroutine

	! Histogram ----------------------------------------------------------------
	subroutine histogram(filename)
		character(*) :: filename

		integer :: i, j, nb
		real(8) :: bx(2), dx
		real(8) :: Ma,Mz
		real(8), allocatable :: X(:), freq(:)

		allocate(X(SimLen))
		open(1, file="results/sim-omega", form="unformatted")
		read(1) X
		close(1)

		nb = floor(2d0*SimLen**0.33d0)
		Ma = minval(X)
		Mz = maxval(X)
		bx = (/ dble(Ma)-1d-12, dble(Mz) /)
		dx = (bx(2)-bx(1))/(nb-1)
		allocate(freq(nb))

		freq = 0
		do i = 1, SimLen
			j = min(max(ceiling((X(i)-bx(1))/dx), 1), nb)
			freq(j) = freq(j)+1
		end do
		freq = freq/SimLen/dx

		open (1, file="results/"//filename);
		do i = 1,nb
			write(1,"(f8.5,f12.6)") bx(1)+dble(i)*dx, freq(i)
		end do
		close(1)
		deallocate(freq, X)
    end subroutine

end
