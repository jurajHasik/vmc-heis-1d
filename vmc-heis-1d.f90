	module heisenberg

	implicit none

	real*8, parameter :: pi = 4.0d0 * atan(1.0d0) 

	! Holds constant Jastrow factors
	real*8, dimension(:,:), allocatable, save :: JF
	! Holds indices for selecting spin-flip pair
	integer, dimension(:), allocatable, save  :: auxI

	contains

	! spin configuration is given as array of integers s where
	! s(i) = -1 => s^z_i s = -1/2 s
	! s(i) = 1  => s^z_i s =  1/2 s

	! Guiding function - variational wavefunction
	function gf(p, s)
		implicit none

		real*8 :: gf
		real*8, dimension(:)  :: p  ! p - array of params
		integer, dimension(:) :: s  ! s - spin conf given in s_z basis 

		integer :: sgn_M            ! Marshall sign	

		integer :: i, L

		L = size(s)
		sgn_M = (-1)**(sum(s(1:L:2)+1)/2) ! sum over even elements of s

		! create outer product of s x s as (s x s)_ij := s_i*s_j
		! multiply element-wise JF * (s*s)
		gf = sgn_M*exp(0.5d0*p(1)*sum((spread(s,dim=2,ncopies=L)&
			&*spread(s,dim=1,ncopies=L))*JF))
	end function gf

	! Local energy
	! Compute <s|H|gf> = \sum_k c_k<s_k|)|gf> = E(s)gf(p,s)
	! with H of AFM Heisenberg model + PBC (interaction J=1) 
	! H = \sum_i [1/2(s^+_i s^-_{i+1} + s^-_i s^+_{i+1}) 
	!               + s^z_i s^z_{i+1}]
	function el(p, s)
		implicit none

		real*8 :: el
		real*8, dimension(:)  :: p
		integer, dimension(:) :: s  ! s - spin conf given in s_z basis

		integer, dimension(size(s)) :: s_c ! cyclic shift of s by 1
        integer :: i

        s_c = CSHIFT(s,1)
		! Compute diagonal part
		el = 0.25d0*DOT_PRODUCT(s, s_c)

		! Compute off diagonal part
		! a) s^+_i s^-_{i+1} brings same contribution as  
		! b) performing spin flip on bond results in aditional -1 within
		!    Marshall sign factor wrt. to original state s
		! c) Jastrow factor can be expressed as original Jastrow factor
		!    + difference due to spin flip
		! Hence the ratio <s_spin-flip|gf>/<s|gf> = -1 * exp(p(1)/2*dif)
		do i=1, size(s)                       ! loop over bonds
			if((s(i)+s_c(i)) .eq. 0) then     ! spin flip possible
				el = el - 0.5d0*exp(-2.0d0*p(1)*(s(i)*sum(s*JF(i,:)) + &
				&s_c(i)*sum(s_c*JF(i,:)) - 2.0d0*JF(1,2)*s(i)*s_c(i)))
			endif
		enddo

	end function el

	! Computes s^z_i s^z_i+x correlation function
	function szsz(p, s, dist)
		implicit none

		real*8 :: szsz
		real*8, dimension(:)  :: p
		integer, dimension(:) :: s
		integer :: dist

		integer, dimension(size(s)) :: s_c ! cyclic shift of s by $dist

		s_c = CSHIFT(s, dist)
		! Take avg of all pairs shifted by $dist 
		szsz = 0.25d0*DOT_PRODUCT(s, s_c)/dble(size(s))
	end function szsz

	! Computes s^p_i s^p_i+x correlation function
	! Q: Should be automatically zero because the variational wavefunction
	! includes projector on singlet subspace ?
	function spsp(p, s, dist)
		implicit none

		real*8 :: spsp
		real*8, dimension(:)  :: p
		integer, dimension(:) :: s
		integer :: dist

		integer :: i, sgn_m
		integer, dimension(size(s)) :: s_c ! cyclic shift of s by $dist

		s_c = CSHIFT(s, dist)

		spsp = 0.0d0
		do i=1, size(s) ! Loop over all pairs and check if we can apply spsp
			if((s(i)+s_c(i)) == -2) then
				sgn_m = (-1)**((1-mod(i,2))+(1-mod((i+dist),2)))
				spsp = spsp + sgn_m*exp(-2.0d0*p(1)*(s(i)*sum(JF(i,:)*s)&
					&+s_c(i)*sum(JF(i,:)*s_c)-2.0d0*JF(1,dist+1)*s(i)*s_c(i)))
			endif
		enddo
		! Take avg of all pairs shifted by $dist
		spsp = spsp/dble(size(s))
	end function spsp

	! Computes s^p_i s^m_i+x correlation function
	function spsm(p, s, dist)
		implicit none

		real*8 :: spsm
		real*8, dimension(:)  :: p
		integer, dimension(:) :: s
		integer :: dist

		integer :: i, sgn_m
		integer, dimension(size(s)) :: s_c ! cyclic shift of s by $dist

		s_c = CSHIFT(s, dist)

		spsm = 0.0d0
		do i=1, size(s) ! Loop over all pairs and check if we can apply spsp
			if((s(i) == -1) .and. (s_c(i) == 1)) then
				sgn_m = (-1)**((1-mod(i,2))+(1-mod((i+dist),2)))
				spsm = spsm + sgn_m*exp(-2.0d0*p(1)*(s(i)*sum(JF(i,:)*s)&
					&+s_c(i)*sum(JF(i,:)*s_c)-2.0d0*JF(1,dist+1)*s(i)*s_c(i)))
			endif
		enddo
		! Take avg of all pairs shifted by $dist
		spsm = spsm/dble(size(s))
	end function spsm

	! Computes the expectation value of <s|O_i O_j|gf>/<s|gf>
	function crossCor(p, s) result (cc)
		implicit none

		real*8, dimension(:) :: p
		integer, dimension(:) :: s
		real*8, dimension(size(p),size(p)) :: cc

		! Implementation of individual operators corresponding to derivatives
		! wrt to variational parameters

		! In case of  gf = sgn_m*exp( \alpha/2 \sum_{i!=j} v_ij*S^z_i*S^z_j )
		! the O_\alpha is simply 1/2 \sum_{i!=j} v_ij*S^z_i*S^z_j
		cc(1,1) = (0.5d0*sum(spread(s,dim=2,ncopies=size(s))&
			&*spread(s,dim=1,ncopies=size(s))*JF))**2.0d0
	end function crossCor

	! Computes the expectation value of <s|O_i|gf>/<s|gf>
	function dparam(p,s) result (dp)
		implicit none

		real*8, dimension(:) :: p
		integer, dimension(:) :: s
		real*8, dimension(size(p)) :: dp

		dp(1) = 0.5d0*sum(spread(s,dim=2,ncopies=size(s))&
			&*spread(s,dim=1,ncopies=size(s))*JF)
	end function dparam

	subroutine move(p, s, rI2, rF, acc)
		implicit none

		real*8, dimension(:) :: p
		integer, dimension(:) :: s
		integer :: acc                  ! #number of accepted moves
		integer, dimension(2) :: rI2    ! Random indices for spin-flip pair
		real*8 :: rF                    ! Random float for acceptance

		integer :: i1, i2
		real*8 :: w

		! I have to allow proposal of moves moving outside of singlet sector
		! and count them as regular moves, even though they are rejected,
		! since the variational function implicitly contains projector on 
		! singlet subspace

		i1 = rI2(1)+1 ! Random number from 1..size(s)
		i2 = mod(i1+1,size(s)+1)+(i1+1)/(size(s)+1)
		!i2 = rI2(2)+1

		! if((s(i1)+s(i2))==0) then ! we stay in singlet sector
		! 	!s(i1)=(-1)*s(i1)
		! 	!s(i2)=(-1)*s(i2)
		! 	! compute weight (<s_new|gf> / <s_old|gf>)^2
		! 	w = exp(-4.0d0*p(1)*(s(i1)*sum(s*JF(i1,:)) + &
		! 		&s(i2)*sum(s*JF(i2,:)) - 2.0d0*JF(i1,i2)*s(i1)*s(i2)))
		! 	if(w > rF) then    ! accepted
		!  		s(i1)=(-1)*s(i1)
		! 		s(i2)=(-1)*s(i2)
		!  		acc=acc+1
		!  	endif
		! endif
		if((s(i1)+s(i2))==0) then ! we stay in singlet sector
			!s(i1)=(-1)*s(i1)
			!s(i2)=(-1)*s(i2)
			! compute weight (<s_new|gf> / <s_old|gf>)^2
			w = exp(-4.0d0*p(1)*(s(i1)*sum(s*JF(i1,:)) + &
				&s(i2)*sum(s*JF(i2,:)) - 2.0d0*JF(i1,i2)*s(i1)*s(i2)))
			if(w > rF) then    ! accepted
				s(i1)=(-1)*s(i1)
				s(i2)=(-1)*s(i2)
				acc=acc+1
			endif
		endif
	end subroutine move

	subroutine init(L,s)
		implicit none

		integer :: L
		integer, dimension(L) :: s  ! s - spin conf given in s_z basis

		integer :: i,j 

		! Init by fixed AFM order
		do i=1, L, 2
			s(i)   = -1
			s(i+1) =  1
		enddo

		! Init by FM order - algorithm is confined to singlet sector
		! s(:) = 1

		allocate(JF(L,L))
		do i=1, L
		do j=1, L
			JF(i,j) = log((2.0d0*sin(pi*abs(i-j)/dble(L)))**2.0d0)
		enddo
			JF(i,i) = 0.0d0
		enddo

		allocate(auxI(L))
		do i=1, L
			auxI(i)=i
		enddo
	end subroutine

end module 

program vmc

	use heisenberg
	use mersenne_twister

	implicit none

	integer, parameter :: optUnit = 99
	integer, parameter :: valUnit = 100
	integer, parameter :: cfUnit = 101

	! Parameters
	integer :: nParams
	real*8, dimension(20) :: pIn
	real*8 :: opt_eps
	integer :: L
	integer :: nEQ, nPROD, nOPT
	NAMELIST /MODEL/ nParams, pIn, L, nEQ, nPROD, nOPT, opt_eps
	integer :: maxX
	NAMELIST /CORRELATION/ maxX
	character (len=20) :: outf, corrf
	NAMELIST /OUTPUTS/ outf, corrf

	! Workaround for allocatable arrays as input in namelist
	real*8, allocatable, dimension(:) :: p  ! p - array of params
	integer, allocatable, dimension(:) :: s ! s - spin conf given 
											!     in s_z basis
	! Averages
	real*8 :: el_avg, el_var
	real*8, allocatable, dimension(:)   :: dp_avg ! AVG of derivative ops
	real*8, allocatable, dimension(:)   :: f_avg  ! AVG of "param" forces 
	real*8, allocatable, dimension(:,:) :: metric ! Connected AVG of <O_i O_j>

	! Accumulators
	real*8, allocatable, dimension(:)     :: el_vals
	real*8, allocatable, dimension(:,:)   :: dp_vals
	real*8, allocatable, dimension(:,:,:) :: cc_vals

	! Correlation functions - averages
	real*8, allocatable, dimension(:) :: szsz_avg, spsp_avg, spsm_avg
	real*8, allocatable, dimension(:,:) :: szsz_vals, spsp_vals, spsm_vals



	integer :: seed
	
	integer :: i,j,x
	integer, dimension(2) :: rI2
	real*8 :: rF
	integer :: acc

	character (len=20) :: dbgFmtL, dbgFmtP, dbgFmtCorr

	integer :: numArg
	character (len=20) :: argUuid, argSeed

	! Read in args
	numArg = IARGC()
	if(numArg .ne. 2) then
    	! Unique id for output & input files
    	write(*,'("Invalid Number of Arguments")')
    	call EXIT(0)
	endif
	call GETARG(1, argUuid)
	call GETARG(2, argSeed)
	read(argSeed, *) seed

	read(*, NML=MODEL )
	write(*,'("L= ",I10)') L
	write(dbgFmtL,'("(",I3,"I3)")') L
	write(*,'("nParams= ",I10)') nParams
	if( nParams .gt. 20 ) then
		write(*,'("Number of parameters exceeds 20")')
		call EXIT(0)
	endif
	write(dbgFmtP,'("(",I3,"1f10.5)")') nParams
	write(*,dbgFmtP) pIn(1:nParams)
	write(*,'("nEQ: ",I10)') nEQ
	write(*,'("nPROD: ",I10)') nPROD
	write(*,'("nOPT: ",I10)') nOPT
	write(*,'("opt_eps: ",1E10.5)') opt_eps
	read(*, NML=CORRELATION )
	read(*, NML=OUTPUTS )

	allocate(p(nParams))
	p = pIn(1:nParams)
	write(*,'("p: ")',advance='no')
	write(*,dbgFmtP) p
	allocate(s(L))
	allocate(el_vals(nPROD))
	allocate(dp_avg(nParams))
	allocate(f_avg(nParams))
	allocate(dp_vals(nParams,nPROD))
	allocate(metric(nParams,nParams))
	allocate(cc_vals(nParams,nParams,nPROD))

	! Correlation functions - we compute correlation only for fixed value
	! of variational parameters, eg. when no further optimization is done  
	if(nOPT.eq.1) then
		write(*,'("Max site-site dist considered for corr. f.: ",I3)') maxX
		allocate(szsz_avg(maxX), spsp_avg(maxX), spsm_avg(maxX))
		allocate(szsz_vals(maxX,nPROD), spsp_vals(maxX,nPROD),&
			& spsm_vals(maxX,nPROD))
	endif

	call random_setseed(seed)
	call init(L,s)	! Init by any random / fixed conf into M=0 subspace
					! Init Mersenne Twister
	write(*,dbgFmtL) s

	open(unit=optUnit, file=trim(argUuid)//"-opt.dat", form="formatted")
	write(optUnit,'("OPT STEP",2X,"el_avg",14X,"p",19X,"f_avg")')	
	write(optUnit,*)

do j=1, nOPT
	acc=0
	do i=1, nEQ*L
		call random_index(L,rI2)
		call random_number(rF)
		call move(p, s, rI2, rF,acc)
		! Take expectation values
		!if(mod(i,L) == 0) then
			!write(*,'("STEP: ",I10," e_l= ",1f20.10," acc= ",1f20.10)') &
			!	&i, el(p,s), dble(acc)/dble(i)
		!endif
	enddo
	write(*,'("nOPT Step: ",I4," TOTAL EQ acceptance= ",1f20.10)') &
		& j, dble(acc)/dble(i)

	acc=0
	el_vals=0.0d0
	dp_vals=0.0d0
	cc_vals=0.0d0
	if(nOPT.eq.1) then
		szsz_vals=0.0d0
		spsp_vals=0.0d0
		spsm_vals=0.0d0
	endif
	do i=1, nPROD*L
		call random_index(L,rI2)
		call random_number(rF)
		call move(p, s, rI2, rF,acc)
		if(mod(i,L) == 0) then
			el_vals(i/L)     = el(p,s)
			dp_vals(:,i/L)   = dparam(p,s)
			cc_vals(:,:,i/L) = crossCor(p,s)
			if(nOPT.eq.1) then
				do x=1, maxX
					szsz_vals(x,i/L) = szsz(p,s,x)
					spsp_vals(x,i/L) = spsp(p,s,x)
					spsm_vals(x,i/L) = spsm(p,s,x)
				enddo
			endif
			!write(*,'("STEP: ",I10," e_l= ",1f20.10," acc= ",1f20.10)') &
			!	&i, el(p,s), dble(acc)/dble(i)
		endif
	enddo

	! Compute statistics
	if(nOPT.eq.1) then
		! The variational parameters are optimized, now we are interested
		! in estimation of errors of observables hence only single nOPT cycle
		
		! Write out samples of loc energy
		open(unit=valUnit, file=trim(argUuid)//"-"//trim(outf),&
			&form="formatted")
		write(valUnit,'("PROD STEP",11X,"e_loc")')
		write(valUnit,*)
		do i=1, nPROD
			write(valUnit,'(I10,10X,1f20.10)') i, el_vals(i)
		enddo
		close(valUnit)
		! End write out samples of loc energy

		if(maxX > 0) then
			! Write out samples of corr fc's
			open(unit=cfUnit, file=trim(argUuid)//"-"//trim(corrf),&
				&form="formatted")
			write(cfUnit,'("PROD STEP",3X,"szsz [",I3," columns]",2X,&
				&"spsm [",I3," columns]")') maxX, maxX
			write(cfUnit,*)
			write(dbgFmtCorr,'("(I10,2X,",I3,"f10.5)")') 2*maxX
			do i=1, nPROD
				! Specific for maxX = 40
				write(cfUnit,dbgFmtCorr) i, szsz_vals(:,i),&
					&spsm_vals(:,i)
			enddo
			close(cfUnit)
			! End write out samples of corr fc's

			! Compute avg of correlation functions
			szsz_avg = sum(szsz_vals,2)/dble(nPROD)
			spsp_avg = sum(spsp_vals,2)/dble(nPROD)
			spsm_avg = sum(spsm_vals,2)/dble(nPROD)

			write(*,*)
			write(*,'("X",4X,"s^z_0s^z_X",5X,"s^p_0s^p_X",5X,&
				&"s^p_0s^m_X")')
			do x=1, maxX
				write(*,'(I5,3f15.10)') x, szsz_avg(x), spsp_avg(x),&
				& spsm_avg(x)
			enddo
			write(*,*)
		endif
	endif

	el_avg = sum(el_vals)/dble(nPROD)
	el_var = DOT_PRODUCT(el_vals-el_avg,el_vals-el_avg)/dble(nPROD)
	write(*,'("p(1)= ",1f15.10," el_avg= ",1f15.10," el_var= ",1f15.10)')&
		& p(1), el_avg, el_var
	if(nOPT.gt.1) then
		! Compute expectation values of derivative operators and cross-
		! correlation matrix necessary to update the values of variational 
		! params
		dp_avg = sum(dp_vals,2)/dble(nPROD) ! Sum over second dimensions
		                                    ! e.g. samples
		metric = sum(cc_vals,3)/dble(nPROD) - spread(dp_avg,dim=2,&
			&ncopies=size(p))*spread(dp_avg,dim=1,ncopies=size(p))
		f_avg = 0.0d0
		do i=1, nPROD
			f_avg = f_avg + dp_vals(:,i)*el_vals(i)
		enddo
		f_avg = dp_avg*el_avg - f_avg/dble(nPROD)
		write(*,'("f_avg= ",1f15.10," metric= ",1f15.10)') f_avg(1),&
			& metric(1,1)
		write(optUnit,'(I10,3f20.10)') j, el_avg, p(1), f_avg(1)
		! Compute new values for params p
		p(1) = p(1) + f_avg(1)/metric(1,1)*min(1.0d0, abs(f_avg(1)*opt_eps))
	endif
enddo
close(optUnit)

end program vmc
