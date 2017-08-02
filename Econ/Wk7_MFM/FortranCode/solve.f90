module solve

    use vars
    use solver
    use spl

    implicit none
    real(8), external :: zero

    integer, parameter :: Nv = 4

	real(8), pointer, private :: MU(:,:,:), MUh(:), MUf(:), qk2(:), w2(:)
	real(8), pointer, private :: spl_muh(:,:,:), spl_muf(:,:,:)
    real(8), private :: prices(2)
    real, private :: t0, t1

contains

    subroutine eqConditions(Narg, arg, err, flag)
		integer, intent(in) :: Narg
		real(8), intent(in) :: arg(Narg)
		real(8), intent(out) :: err(Narg)
		integer, intent(inout) :: flag

		integer :: iz
		real(8) :: bh, ch, cf, mh, mf, kh, kf, qk, qb, uh, Yw
        
		ch = AD/(1d0 + exp(-eta*arg(1)))    ! H's consumption
		cf = AD - ch                        ! F's consumption
		kh = max(arg(2),0d0)**nu            ! H's portfolio
		kf = max(arg(3),0d0)**nu+kmin       ! H's portfolio
		bh = arg(4)							! H bond position

        mh = max(-arg(2),0d0)**nu           ! H's Lagrange multiplier on kh>=0
        mf = max(-arg(3),0d0)**nu           ! F's Lagrange multiplier on kf>=0

		! Solve for the next period wealth share
		do iz = 1, Nz
			zw = iz
            Yw = (Ah*kh+Af*kf)/(kh+kf) + Zs(iz)

			prm = [ Yw,kh,bh ]
			if (wEvolErr(w_bnd(2)) > 0d0) then
				w2(iz) = wGrid(Nw)
			elseif (wEvolErr(w_bnd(1)) < 0d0) then
				w2(iz) = w_bnd(1)
			else
                w2(iz) = zero( w_bnd(1), w_bnd(2), epsilon(w2), tol_w, wEvolErr )
            end if
			qk2(iz) = prc(1)

            ! Compute next period MUs for asset valuation
			MUh(iz) = csval_y(w2(iz), Yw, yGrid, brk_w, spl_muh(:,:,:))
			MUf(iz) = csval_y(w2(iz), Yw, yGrid, brk_w, spl_muf(:,:,:))
        end do

        uh = ch**mgammah
		qk = (sum(pz*(MUh*(Ah+Zs+qk2))) + mh) / uh   ! Price of H tree (H's EE)
		qb = (sum(pz*(MUh))) / uh                    ! Price of bond

		err(1) = ch + qk*kh + qb*bh - omega*(AD+qk)  ! Budget constraint

		err(2) = (sum(pz*(MUf*(Af+Zs+qk2))) + mf)/qk ! Price of H tree (F's EE)
		err(2) = cf - err(2)**mgammainvf
		err(3) = (sum(pz*(MUf))) / qb                ! Price of bond
		err(3) = cf - err(3)**mgammainvf

		err(4) = 1d0 - kh - kf                       ! Error in mkt clearing for H tree
		prices = [ qk, qb ]                          ! Record eqm prices
        flag = 0
	end subroutine

	function wEvolErr(arg) result(err)
		real(8) :: arg, err
        real(8) :: qk, Yw
 
        Yw = prm(1)
		qk = csval_y(arg, Yw, yGrid, brk_w, spl_qk(:,:,:))
        prc(1) = qk

		err = (qk + Ah + Zs(zw))*prm(2) + prm(3)
		err = err/(Yw + qk) - arg
	end function
 
	function csval_y(w, y, yGrid, brk, spl) result(res)
		real(8) :: w, y, yGrid(:), brk(:), spl(:,:,:), res
        
        integer :: iy, Ny
        real(8) :: res1, res2
 
        Ny = size(yGrid)
        iy = min(max(sum(-1*(yGrid < y)), 1), Ny-1)
		res1 = csval(brk, spl(:,:,iy), w)
		res2 = csval(brk, spl(:,:,iy+1), w)

        if (abs(yGrid(2)-yGrid(1)) > eps) then
            res = res1 + (res2-res1)*(y-yGrid(1))/(yGrid(2)-yGrid(1))
        else
            res = res1
        end if
    end function

    subroutine compute_node(wstate, ystate, guess, sol)
	    integer :: wstate, ystate
	    real(8) :: guess(Nv), sol(Nv)

        integer :: info
        real(8) :: diag(Nv)

	    real(8) :: err(Nv), max_err
 
        omega = wGrid(wstate)           ! Fix wealth share
        AD    = yGrid(ystate)           ! Fix dividend state

        sol = guess
        call hbrd(eqConditions, Nv, sol, err, tol, 1d3*tol, info, diag)
		max_err = maxval(abs(err))

		sol_t0(:,wstate,ystate) = sol
		pol_t0(:,wstate,ystate) = [ 1d0/(1d0+exp(-eta*sol(1))), [0d0,kmin]+max(sol(2:3),0d0)**nu, sol(4)]
		prc_t0(:,wstate,ystate) = prices
		trans (:,wstate,ystate) = w2
		err_t0(:,wstate,ystate) = err
    end subroutine

    subroutine compute(print_flag, dist)
		integer :: print_flag
        real(8) :: dist
 
        integer :: iy, iw
        real(8) :: ch, cf
        real(8) :: sol(nv), sol_m(nv), x0(nv)
 
        if (print_flag > 0) then; print "(i4,a\)", periods+1, ":"; end if
 
        ! Save previous solution
        sol_t1 = sol_t0
        pol_t1 = pol_t0
        prc_t1 = prc_t0
 
        ! Create spline interpolants for policies and prices
        do iy = 1, Ny
			do iw = 1, Nw
				ch = pol_t0(1,iw,iy)*yGrid(iy)
				cf = yGrid(iy) - ch
				MU(1,iw,iy) = betah*ch**mgammah
				MU(2,iw,iy) = betaf*cf**mgammaf
			end do
            call csint_ext(omega_ext, wGrid, MU(1,:,iy), brk_w, spl_muh(:,:,iy))
			call csint_ext(omega_ext, wGrid, MU(2,:,iy), brk_w, spl_muf(:,:,iy))
			call csint_ext(omega_ext, wGrid, prc_t0(1,:,iy), brk_w, spl_qk(:,:,iy))
        end do

        call cpu_time(t0)
        do iy = 1, Ny
	        x0 = sol_t0(:,M,iy)
	        call compute_node(M,iy, x0, sol_M)
 
	        x0 = sol_M
	        do iw = M-1, 1, -1
		        call compute_node(iw,iy, x0, sol)
		        x0 = sol
            end do
 
	        x0 = sol_M
	        do iw = M+1, Nw
		        call compute_node(iw,iy, x0, sol)
		        x0 = sol
	        end do
	        if (print_flag>0) then; print "(a,i0\)", " ", iy; end if
        end do
        call cpu_time(t1)
 
        dist = maxval(abs(prc_t1-prc_t0))
        if (print_flag>0) then; print "(a,f10.7,a,f6.2)", ", dist = ", dist, ", time = ", t1-t0; end if
    end subroutine

    subroutine csint_ext(dxe, x, y, brk, spl)
    ! Create a spline on [xmin-dx,xmax+dx]
        real(8) :: dxe, x(:), y(:), brk(size(x,1)+4), spl(4,size(x,1)+4)
        real(8) :: brk0(size(x,1)), spl0(4,size(x,1))
 
        integer :: kx
        real(8) :: d1, d2, dx, xe(size(x,1)+4), ye(size(x,1)+4), yl(2), yh(2)
 
        call csint(x, y, brk0, spl0)
        
        yl(1) = y(1) - spl0(2,1)*dxe/1 + spl0(3,1)*dxe*dxe/2
        yl(2) = y(1) - spl0(2,1)*dxe/2 + spl0(3,1)*dxe*dxe/8
 
        kx = size(x,1)
        dx = x(kx)-x(kx-1)
        d1 = spl0(2,kx-1) + spl0(3,kx-1)*dx + spl0(4,kx-1)*dx*dx/2
        d2 = spl0(3,kx-1) + spl0(4,kx-1)*dx
        yh(1) = y(kx) + d1*dxe/2 + d2*dxe*dxe/8
        yh(2) = y(kx) + d1*dxe/1 + d2*dxe*dxe/2
 
        xe = [x(1)-dxe, x(1)-dxe/2, x, x(kx)+dxe/2, x(kx)+dxe]
        ye = [yl, y, yh]
        call csint(xe, ye, brk, spl)
    end subroutine

    subroutine simulate(w0)
        real(8) :: w0
        real(8) :: kh, n0, y0
        real(8), pointer :: wSim(:)
        integer :: iy, iz1, it, z0, z1

        allocate(wSim(SimLen))
        do iy = 1, Ny
			call csint_ext(omega_ext, wGrid, pol_t0(2,:,iy), brk_w, spl_kh(:,:,iy))
            do iz1 = 1, Nz
			    call csint_ext(omega_ext, wGrid, trans(iz1,:,iy), brk_w, spl_wp(:,:,iy,iz1))
            end do
        end do

        z0 = 1
        n0 = w0
        y0 = (Af+Ah)/2d0 + Zs(z0)
        do it = 1, SimLen
            if (ushock(it) < 0.5d0) then
                z1 = 1
            else
                z1 = 2
            end if

            wSim(it) = min(max(csval_y(n0, y0, yGrid, brk_w, spl_wp(:,:,:,z1)), 0.001d0), 0.999d0)
            kh = csval_y(n0, y0, yGrid, brk_w, spl_kh(:,:,:))
            y0 = Ah*kh+Af*(1d0-kh)+Zs(z1)
            n0 = wSim(it)
            z0 = z1
        end do

        open (201, file = 'results/sim-omega.txt')
        write(201, '(f12.6)') wSim
        close(201)

        open (201, file = 'results/sim-omega', form='unformatted')
        write(201) wSim
        close(201)
        deallocate(wSim)
    end subroutine

    subroutine compute_prep()
        integer :: ec
 
		allocate(spl_muh(4,Nw+4,Ny), spl_muf(4,Nw+4,Ny), MU(2,Nw,Ny), MUh(Nz), MUf(Nz), &
            qk2(Nz), w2(Nz), stat = ec)
    end subroutine
 
    subroutine compute_free()
        integer :: ec
 
        deallocate(spl_muh, spl_muf, MU, MUh, MUf, qk2, w2, spl_qk, stat = ec)
    end subroutine

end
