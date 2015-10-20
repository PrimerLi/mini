 Complex(kind = 8) function f(gamman, epsilon2, xq)
    implicit none
    complex(kind = 8), intent(in) :: gamman
    real(kind = 8), intent(in) :: epsilon2, xq

    real(kind = 8) :: energyUpper, energyLower, denergy
    integer, parameter :: dp = 8
    integer :: ienergy, Nenergy
    real(kind = 8) :: epsilon1
    real(kind = 8), parameter :: pi = acos(-1.0_dp)
    complex(kind = 8) :: csum

    Nenergy = 1000
    energyUpper = 20.0_dp
    energyLower = -energyUpper
    denergy = (energyUpper - energyLower)/dble(Nenergy)
    csum = dcmplx(0, 0)

    do ienergy = 1, Nenergy
    	epsilon1 = energyLower + (ienergy - 0.5_dp)*denergy
	csum = csum + denergy*1.0_dp/(sqrt(pi*(1.d0 - xq**2)))*1.0_dp/(gamman - epsilon1)*exp(-(epsilon1 - xq*epsilon2)**2/(1.d0 - xq**2))
    enddo

    f = csum

 end function f

 program main
 implicit none
    interface 
    	subroutine inverse(A, A_inv)
    		complex(kind = 8), intent(in) :: A(:, :)
    		complex(kind = 8), intent(out) :: A_inv(:, :)
    	end subroutine inverse
    end interface
 	integer :: narg
	character *100 :: buffer
	character *20 :: sigma_file
	character *20 :: mu_file
	character *20 :: griomfile
	integer, parameter :: dp = 8	
	integer :: Niom, Niom2
	Integer :: NiomMin = 100
	real(kind = 8) :: beta
	integer :: Niom_factor = 5
	integer :: Niom2_factor = 2
	real(kind = 8), dimension(:), allocatable :: xiom
	real(kind = 8), dimension(:), allocatable :: xiom2
	complex(kind = 8) :: ii = dcmplx(0, 1)
    	complex(kind = 8), dimension(:, :, :), allocatable :: sigma
	complex(kind = 8), dimension(:, :, :), allocatable :: sigma_temp
	integer :: norb = 2
	integer :: row, col
	integer :: nxq = 2
	integer :: ixq
	real(kind = 8) :: xqmin
	real(kind = 8) :: xqmax
	real(kind = 8) :: dxq
    	complex(kind = 8), dimension(:, :, :, :), allocatable :: Chi_0_charge, Chi_0_spin, Chi_0_pair
	integer :: Nenergy = 600
	real(kind = 8), dimension(:), allocatable :: energy
	real(kind = 8) :: energyUpper
	real(kind = 8) :: energyLower
	real(kind = 8) :: denergy
    	Integer :: L_dim = 4
	real(kind = 8) :: eps_f, chem
	real(kind = 8) :: Rv = 0.955_dp
	complex(kind = 8), dimension(:, :, :), allocatable :: Chi_int
	real(kind = 8), dimension(:), allocatable :: xqarray
	complex(kind = 8), dimension(:, :), allocatable :: f_z, g_z
	integer :: i, j, mu, nu, ku, lu
	integer :: nw, nw1
	real(kind = 8) :: re, reerror, imag, imagerror
	complex(kind = 8) :: alpha_n, beta_n, gamma_n
	real(kind = 8) :: x1
	complex(kind = 8), dimension(:, :), allocatable :: Z_hlp, Z_hlp1, Z_hlp2, Z_hlp3
	complex(kind = 8), dimension(:, :), allocatable :: Z_hlp1_inv, Z_hlp2_inv, Z_hlp3_inv
	complex(kind = 8) :: czero = dcmplx(0, 0)
    	integer :: ienergy
	complex(kind = 8) :: f
	real(kind = 8), parameter :: pi = acos(-1.0_dp)
    	Complex(kind = 8), dimension(:, :, :), allocatable :: Chi_0_i_charge, Chi_0_i_spin, Chi_0_i_pair
	complex(kind = 8), dimension(:, :, :), allocatable :: g_temp
	complex(kind = 8), dimension(:, :, :), allocatable :: giom
	integer :: Ndim
	complex(kind = 8), dimension(:, :, :, :), allocatable :: Chi_i_charge, Chi_i_spin, Chi_i_pair
	complex(kind = 8) :: Z, Z1, Z2, Z3
	complex(kind = 8), dimension(:, :, :, :), allocatable :: Chi_0_l_charge, Chi_0_l_spin, Chi_0_l_pair
	complex(kind = 8), dimension(:, :), allocatable :: gamma_charge, gamma_spin, gamma_pair
	integer :: k, l
	complex(kind = 8), dimension(:, :, :), allocatable :: Chi_int_charge, Chi_int_spin, Chi_int_pair
	complex(kind = 8), dimension(:, :, :), allocatable :: M_charge, M_spin, M_pair
	complex(kind = 8) :: s 

	narg = command_argument_count()
    	if (narg .ne. 4) then
		write(*, *) "Argument number is wrong."
		write(*, *) "sigma_file is the first argument"
		write(*, *) "beta is the second argument. "
		write(*, *) "mu_file is the third argument. "
		write(*, *) "griomfile is the fourth argument"
		stop
	endif
	call getarg(1, buffer)
    	read(buffer, *) sigma_file
    	call getarg(2, buffer)
	read(buffer, *) beta
	call getarg(3, buffer)
    	read(buffer, *) mu_file
	call getarg(4, buffer)
    	read(buffer, *) griomfile


	Niom = max(int(Niom_factor*beta), NiomMin)
    	Niom2 = max(int(beta*Niom2_factor), 20)
    	Ndim = 2*Niom2
    	allocate(xiom(Niom), xiom2(-Niom+1:Niom))

    	allocate(xqarray(nxq))
    	xqmin = -1.0_dp
	xqmax = -xqmin
	dxq = (xqmax - xqmin)/dble(nxq-1)
    	do ixq = 1, nxq
		xqarray(ixq) = xqmin + (ixq - 1)*dxq
	enddo

	Allocate(Chi_int(nxq, L_dim, L_dim))
    	allocate(f_z(norb, norb), g_z(norb, norb))
    	allocate(sigma(Niom, norb, norb))
    	allocate(sigma_temp(-Niom+1:Niom, norb, norb))
    	allocate(Chi_0_charge(nxq, -Niom2+1:Niom2, L_dim, L_dim))
    	allocate(Chi_0_spin(nxq, -Niom2+1:Niom2, L_dim, L_dim))
    	allocate(Chi_0_pair(nxq, -Niom2+1:Niom2, L_dim, L_dim))
    	allocate(Z_hlp2(L_dim*Ndim, L_dim*Ndim), Z_hlp3(L_dim*Ndim, L_dim*Ndim))
    	allocate(Z_hlp1_inv(L_dim*Ndim, L_dim*Ndim), Z_hlp2_inv(L_dim*Ndim, L_Dim*Ndim), Z_hlp3_inv(L_dim*Ndim, L_dim*Ndim))
    	allocate(energy(Nenergy))
    	allocate(Chi_0_i_charge(-Niom2+1:Niom2, L_dim, L_dim))
    	allocate(Chi_0_i_spin(-Niom2+1:Niom2, L_dim, L_dim))
    	allocate(Chi_0_i_pair(-Niom2+1:Niom2, L_dim, L_dim))
    	allocate(giom(Niom, norb, norb))
    	allocate(g_temp(-Niom+1:Niom, norb, norb))
    	allocate(Chi_i_charge(-Niom2+1:Niom2, -Niom2+1:Niom2, L_dim, L_dim))
    	allocate(Chi_i_spin(-Niom2+1:Niom2, -Niom2+1:Niom2, L_dim, L_dim))
    	allocate(Chi_i_pair(-Niom2+1:Niom2, -Niom2+1:Niom2, L_dim, L_dim))
    	allocate(gamma_charge(L_dim*Ndim, L_dim*Ndim))
    	allocate(gamma_spin(L_dim*Ndim, L_dim*Ndim))
    	allocate(gamma_pair(L_dim*Ndim, L_dim*Ndim))
    	allocate(Chi_int_charge(nxq, L_dim, L_dim), Chi_int_spin(nxq, L_dim, L_dim))
    	allocate(Chi_int_pair(nxq, L_dim, L_dim))
    	allocate(M_charge(nxq, L_dim*Ndim, L_dim*Ndim))
    	allocate(M_spin(nxq, L_dim*Ndim, L_dim*Ndim))
    	allocate(M_pair(nxq, L_dim*Ndim, L_dim*Ndim))

    	energyUpper = 20.0_dp
	energyLower = -energyUpper
    	denergy = (energyUpper - energyLower)/dble(Nenergy)
    	do ienergy = 1, Nenergy
		energy(ienergy) = energyLower + (ienergy - 0.5_dp)*denergy
	enddo
    	
	open(unit = 17, file = sigma_file, action = "read")
	do i = 1, norb
	do j = 1, norb
		read(17, *) row, col
		do nw = 1, Niom
			read(17, *) xiom(nw), re
			read(17, *) imag
			sigma(nw, row, col) = dcmplx(re, imag)
		enddo
	enddo
	enddo
	close(17)

    	open(unit = 23, file = griomfile, action = "read")
	do i = 1, norb
	do j = 1, norb
		read(23, *) row, col
		do nw = 1, Niom
			read(23, *) xiom(nw), re, reerror, imag, imagerror
			giom(nw, row, col) = dcmplx(re, imag)
		enddo
	enddo
	enddo
	close(23)

    	do i = 1, norb
	do j = 1, norb
		do nw  = 1, Niom
			sigma_temp(nw, i, j) = sigma(nw, i, j)
    			sigma_temp(-nw+1, i, j) = conjg(sigma(nw, i, j))
    			g_temp(nw, i, j) = giom(nw, i, j)
    			g_temp(-nw + 1, i, j) = conjg(giom(nw, i, j))
		enddo
	enddo
	enddo

	do nw = -Niom+1, Niom
		xiom2(nw) = dble(2*(nw-1) + 1)*pi/beta
	enddo

    	open(unit = 19, file = mu_file, action = "read")
	read(19, *) chem
	close(19)
    	eps_f = -0.136

	allocate(Z_hlp(L_dim, L_dim))
    	allocate(Z_hlp1(L_dim, L_dim))
    	do ixq = 1, nxq
	do nw = -Niom2 + 1, Niom2
		x1 = xqarray(ixq)
    		alpha_n = ii*xiom2(nw) - (eps_f - chem) - sigma_temp(nw, 2, 2)
    		beta_n = ii*xiom2(nw) + chem - sigma_temp(nw, 1, 1)
    		gamma_n = beta_n - (Rv + sigma_temp(nw, 1, 2))*(Rv + sigma_temp(nw, 2, 1))/alpha_n
		Z_hlp = czero
		Z_hlp1 = czero
		if (.true.) then
		do ienergy = 1, Nenergy
			f_z(1, 1) = 1.0_dp/(gamma_n - energy(ienergy)*xqarray(ixq)) 
    			f_z(1, 2) = (Rv + sigma_temp(nw, 1, 2))/(alpha_n*(gamma_n - xqarray(ixq)*energy(ienergy)))
    			f_z(2, 1) = (Rv + sigma_temp(nw, 2, 1))/(alpha_n*(gamma_n - xqarray(ixq)*energy(ienergy)))
    			f_z(2, 2) = (beta_n -xqarray(ixq)*energy(ienergy))/(alpha_n*(gamma_n - xqarray(ixq)*energy(ienergy)))
    			g_z(1, 1) = 1.0_dp/(gamma_n - energy(ienergy))
    			g_z(1, 2) = (Rv + sigma_temp(nw, 1, 2))/alpha_n*(1.0_dp/(gamma_n - energy(ienergy)))
    			g_z(2, 1) = (Rv + sigma_temp(nw, 2, 1))/alpha_n*(1.0_dp/(gamma_n - energy(ienergy)))
    			g_z(2, 2) = (beta_n - energy(ienergy))/(alpha_n*(gamma_n - energy(ienergy)))
    			i = 0
			do nu = 1, norb
			do mu = 1, norb
				i = i + 1
				j = 0
				do ku = 1, norb
				do lu = 1, norb
					j = j + 1
					Z_hlp(i, j) = Z_hlp(i, j) + denergy*(1.0_dp/sqrt(pi))&
					&*exp(-energy(ienergy)**2)*f_z(nu, ku)*g_z(lu, mu)
    					Z_hlp1(i, j) = Z_hlp1(i, j) + denergy*(1.0_dp/sqrt(pi))*&
					&exp(-energy(ienergy)**2)*conjg(f_z(nu, ku))*g_z(mu, lu)
				enddo
				enddo
			enddo
			enddo
		enddo
		endif
		do i = 1, L_dim
		do j = 1, L_dim
			Chi_0_charge(ixq, nw, i, j) = -2.0_dp*Z_hlp(i, j)
    			Chi_0_spin(ixq, nw, i, j) = -1.0_dp*Z_hlp(i, j)
    			Chi_0_pair(ixq, nw, i, j) = Z_hlp1(i, j)
		enddo
		enddo
	enddo
	enddo
	deallocate(Z_hlp, Z_hlp1)

    	do nw = -Niom2+1, Niom2
	i = 0
	do nu = 1, norb
	do mu = 1, norb
		i = i + 1
		j = 0
		do ku = 1, norb
		do lu = 1, norb
			j = j + 1
			Chi_0_i_charge(nw, i, j) = -2.0_dp*g_temp(nw, nu, ku)*g_temp(nw, lu, mu)
    			Chi_0_i_spin(nw, i, j) = -g_temp(nw, nu, ku)*g_temp(nw, lu, mu)
    			Chi_0_i_pair(nw, i, j) = g_temp(-nw+1, nu, ku)*g_temp(nw, mu, lu)
		enddo
		enddo
	enddo
	enddo
	enddo

	open(unit = 29, file = "Chi_c_average", action = "read")
	open(unit = 31, file = "Chi_m_average", action = "read")
	open(unit = 37, file = "Chi_p_average", action = "read")
	do nu = 1, L_dim
	do mu = 1, L_dim
		read(29, *) row, col
		read(31, *) row, col
		read(37, *) row, col
		do nw = -Niom2 + 1, Niom2
		do nw1 = 1, Niom2
			read(29, *) xiom2(nw), xiom2(nw1), Z1
			Chi_i_charge(nw, nw1, nu, mu) = Z1
    			Chi_i_charge(-nw + 1, -nw1 + 1, nu, mu) = conjg(Z1)
    			read(31, *) xiom2(nw), xiom2(nw1), Z2
			Chi_i_spin(nw, nw1, nu, mu) = Z2
			Chi_i_spin(-nw + 1, -nw1 + 1, nu, mu) = conjg(Z2)
    			read(37, *) xiom2(nw), xiom2(nw1), Z3
			Chi_i_pair(nw, nw1, nu, mu) = Z3
			Chi_i_pair(-nw + 1, -nw1 + 1, nu, mu) = conjg(Z3)
		enddo
		enddo
	enddo
	enddo
	close(37)
	close(31)
	close(29)

    	allocate(Z_hlp1(L_dim*Ndim, L_dim*Ndim))
    	Z_hlp1 = czero
	Z_hlp1_inv = czero
	Z_hlp3 = czero
	Z_hlp3_inv = czero

	do i = 1, L_dim
	do j = 1, L_dim
		do nw = -Niom2 + 1, Niom2
			k = (i - 1)*Ndim + nw + Niom2
			l = (j - 1)*Ndim + nw + Niom2
			Z_hlp1(k, l) = conjg(Chi_0_i_charge(nw, i, j))
		enddo
		do nw = -Niom2+1, Niom2
		do nw1 = -Niom2+1, Niom2
			k = (i - 1)*Ndim + nw + Niom2
			l = (j - 1)*Ndim + nw1 + Niom2
			Z_hlp3(k, l) = Chi_i_charge(nw, nw1, i, j)
		enddo
		enddo
	enddo
	enddo

	call inverse(Z_hlp1, Z_hlp1_inv)
    	call inverse(Z_hlp3, Z_hlp3_inv)
	gamma_charge = Z_hlp1_inv - Z_hlp3_inv

    	
	do i = 1, L_dim
	do j = 1, L_dim
		do nw = -Niom2 + 1, Niom2
			k = (i - 1)*Ndim + nw + Niom2
			l = (j - 1)*Ndim + nw + Niom2
			Z_hlp1(k, l) = conjg(Chi_0_i_spin(nw, i, j))
		enddo
		do nw = -Niom2+1, Niom2
		do nw1 = -Niom2+1, Niom2
			k = (i - 1)*Ndim + nw + Niom2
			l = (j - 1)*Ndim + nw1 + Niom2
			Z_hlp3(k, l) = Chi_i_spin(nw, nw1, i, j)
		enddo
		enddo
	enddo
	enddo

	call inverse(Z_hlp1, Z_hlp1_inv)
    	call inverse(Z_hlp3, Z_hlp3_inv)
	gamma_spin = Z_hlp1_inv - Z_hlp3_inv


	do i = 1, L_dim
	do j = 1, L_dim
		do nw = -Niom2 + 1, Niom2
			k = (i - 1)*Ndim + nw + Niom2
			l = (j - 1)*Ndim + nw + Niom2
			Z_hlp1(k, l) = conjg(Chi_0_i_pair(nw, i, j))
		enddo
		do nw = -Niom2+1, Niom2
		do nw1 = -Niom2+1, Niom2
			k = (i - 1)*Ndim + nw + Niom2
			l = (j - 1)*Ndim + nw1 + Niom2
			Z_hlp3(k, l) = Chi_i_pair(nw, nw1, i, j)
		enddo
		enddo
	enddo
	enddo

	call inverse(Z_hlp1, Z_hlp1_inv)
    	call inverse(Z_hlp3, Z_hlp3_inv)
	gamma_pair = Z_hlp3_inv - Z_hlp1_inv


    	Chi_int_charge = czero
	do ixq = 1, nxq
		Z_hlp2 = czero
		Z_hlp2_inv = czero
		do i = 1, L_dim
		do j = 1, L_dim
			do nw = -Niom2+1, Niom2
				k = (i - 1)*Ndim + nw + Niom2
				l = (j - 1)*Ndim + nw + Niom2
				Z_hlp2(k, l) = conjg(Chi_0_charge(ixq, nw, i, j))
			enddo
		enddo
		enddo
		do i = 1, L_dim*Ndim
		do j = 1, L_dim*Ndim
			s = czero
			do k = 1, L_dim*Ndim
				s = s + gamma_charge(i, k) * Z_hlp2(k, j)
			enddo
			M_charge(ixq, i, j) = s 
		enddo
		enddo
		call inverse(Z_hlp2, Z_hlp2_inv)
    		Z_hlp1 = Z_hlp2_inv - gamma_charge
		call inverse(Z_hlp1, Z_hlp1_inv) !Z_hlp1_inv is Chi_charge
		do i = 1, L_dim
		do j = 1, L_dim
			s = dcmplx(0, 0)
			do nw = -Niom2+1, Niom2
			do nw1 = -Niom2+1, Niom2
				k = (i - 1)*Ndim + nw + Niom2
				l = (j - 1)*Ndim + nw1 + Niom2
				Chi_int_charge(ixq, i, j) = Chi_int_charge(ixq, i, j) + Z_hlp1_inv(k, l)
			enddo
			enddo
		enddo
		enddo
	enddo

	Chi_int_charge = Chi_int_charge/beta

	Chi_int_spin = czero
	do ixq = 1, nxq
		Z_hlp2 = czero
		Z_hlp2_inv = czero
		do i = 1, L_dim
		do j = 1, L_dim
			do nw = -Niom2+1, Niom2
				k = (i - 1)*Ndim + nw + Niom2
				l = (j - 1)*Ndim + nw + Niom2
				Z_hlp2(k, l) = conjg(Chi_0_spin(ixq, nw, i, j))
			enddo
		enddo
		enddo
		do i = 1, L_dim*Ndim
		do j = 1, L_dim*Ndim
			s = czero
			do k = 1, L_dim*Ndim
				s = s + gamma_spin(i, k) * Z_hlp2(k, j)
			enddo
			M_spin(ixq, i, j) = s 
		enddo
		enddo
		call inverse(Z_hlp2, Z_hlp2_inv)
    		Z_hlp1 = Z_hlp2_inv - gamma_spin
		call inverse(Z_hlp1, Z_hlp1_inv) !Z_hlp1_inv is Chi_l
		do i = 1, L_dim
		do j = 1, L_dim
			do nw = -Niom2+1, Niom2
			do nw1 = -Niom2+1, Niom2
				k = (i - 1)*Ndim + nw + Niom2
				l = (j - 1)*Ndim + nw1 + Niom2
				Chi_int_spin(ixq, i, j) = Chi_int_spin(ixq, i, j) + Z_hlp1_inv(k, l)
			enddo
			enddo
		enddo
		enddo
	enddo

	Chi_int_spin = Chi_int_spin/beta

	Chi_int_pair = czero
	do ixq = 1, nxq
		Z_hlp2 = czero
		Z_hlp2_inv = czero
		do i = 1, L_dim
		do j = 1, L_dim
			do nw = -Niom2+1, Niom2
				k = (i - 1)*Ndim + nw + Niom2
				l = (j - 1)*Ndim + nw + Niom2
				Z_hlp2(k, l) = conjg(Chi_0_pair(ixq, nw, i, j))
			enddo
		enddo
		enddo
		do i = 1, L_dim*Ndim
		do j = 1, L_dim*Ndim
			s = czero
			do k = 1, L_dim*Ndim
				s = s + gamma_pair(i, k) * Z_hlp2(k, j)
			enddo
			M_pair(ixq, i, j) = s 
		enddo
		enddo
		call inverse(Z_hlp2, Z_hlp2_inv)
    		Z_hlp1 = Z_hlp2_inv + gamma_pair
		call inverse(Z_hlp1, Z_hlp1_inv) !Z_hlp1_inv is Chi_l
		do i = 1, L_dim
		do j = 1, L_dim
			do nw = -Niom2+1, Niom2
			do nw1 = -Niom2+1, Niom2
				k = (i - 1)*Ndim + nw + Niom2
				l = (j - 1)*Ndim + nw1 + Niom2
				Chi_int_pair(ixq, i, j) = Chi_int_pair(ixq, i, j) + Z_hlp1_inv(k, l)
			enddo
			enddo
		enddo
		enddo
	enddo

	Chi_int_pair = Chi_int_pair/beta

	open(unit = 23, file = "Chi_int", action = "write")
	write(23, *) "Charge susceptibility: "
	do i = 1, L_dim
	do j = 1, L_dim
		write(23, *) i, "    ", j
		do ixq = 1, nxq
			write(23, *) xqarray(3 - ixq), dble(Chi_int_charge(3 - ixq, i, j))
    		enddo
	enddo
	enddo

	write(23, *) "Magnetic susceptibility: "
	do i = 1, L_dim
	do j = 1, L_dim
		write(23, *) i, "    ", j
		do ixq = 1, nxq
			write(23, *) xqarray(3 - ixq), dble(Chi_int_spin(3 - ixq, i, j))
		enddo
	enddo
	enddo

	write(23, *) "Pair susceptibility: "
	do i = 1, L_dim
	do j = 1, L_dim
		write(23, *) i, "    ", j
		do ixq = 1, nxq
			write(23, *) xqarray(3 - ixq), dble(Chi_int_pair(3 - ixq, i, j))
		enddo
	enddo
	enddo
	close(23)

    	Deallocate(Chi_int)
    	deallocate(f_z, g_z)
    	deallocate(sigma, sigma_temp)
    	deallocate(Chi_0_charge, Chi_0_spin, Chi_0_pair)
    	deallocate(xiom, xiom2)
    	deallocate(Z_hlp1, Z_hlp2, Z_hlp3)
    	deallocate(Z_hlp1_inv, Z_hlp2_inv, Z_hlp3_inv)
    	deallocate(energy)
    	deallocate(Chi_0_i_charge, Chi_0_i_spin, Chi_0_i_pair)
    	deallocate(giom, g_temp)
    	deallocate(Chi_i_charge, Chi_i_spin, Chi_i_pair)
    	deallocate(gamma_charge, gamma_spin, gamma_pair)
    	deallocate(Chi_int_charge, Chi_int_spin, Chi_int_pair)
    	deallocate(M_charge, M_spin, M_pair)

 end program main
