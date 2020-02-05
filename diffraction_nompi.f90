!!!  HISTORY: v7: - to preserve memory removed psiA, psiO, psiC, psiD variables in favour of pfftin,pfftout
!!!                 and deallocate unnecessary arrays 
!!!               - removed allocatable from definition of input in write_it_2_fits => able to invoke write_it_2_fits(real(array,kind=4)...)
!!!           v8: - constante phase factors like exp(pi*i*r^2/lam/z0) -> C_FLOAT_COMPLEX
module parameters
  double precision :: z0, f, z1, z1_0, LA2, rA, rIO, rC, lam, d, l, Ap_size, tilty, dalpha
  integer          :: JJ
  character(len=10)  :: IOstr
  character(len=100) :: ext
  character(len=38)  :: filename_psiA
  character(len=35)   :: filename_scatter
  DOUBLE PRECISION,PARAMETER            :: dpi=3.141592653589793238462643383279502884197169399375105d0
  DOUBLE PRECISION,PARAMETER            :: RSun=16.00/60.0/180.0*dpi   ! max=16.258 arcmin, min=
contains
  subroutine initialization
    z0  =144348.0d0          ! ISD distance
    f   =330.341d0           ! focal length
    z1  =331.143d0           ! PO-IO distance
    z1_0=331.143d0           ! expected PI-IO distance by design. Determines d, l etc
    d   =331.143             ! z1         --- distance plane O'- plane C; defines f=z1/2 of the O2
    l   =2*331.143-f         ! 2*z1-f     --- distance plane C - plane D; defines f=z1/2 of the O3 
    LA2 =35.0                ! size of the aperture plane 2*LA2 in mm
    JJ  =4096                ! array size in pixels in each plane
    rA  =25.0
    rIO =1.748
    IOstr = "mmH"
    rC  =24.25               ! usually 0.97*rA; 25.0->24.25; 24.9->24.15
    lam =5.5d-4              ! lambda in mm
    Ap_size = dpi*(rA/(LA2*2)*JJ)**2         ! number of pixels in the Aperture
    tilty =  0.0/3600.0*dpi/180.0            ! arcsec - tilt of the coronagraph around Oy, converted to radians
    dalpha = 0.0/3600.0*dpi/180.0            ! arcsec - shift of the Sun in +x, converted to radians
    ext = ".cart8"                           ! _hxyABC ; .z1p60mkm
    filename_psiA    = "psi_onaxis.z144348.262144pts.A900.fits"
  end subroutine initialization 
end module

program diffraction
  use parameters
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
   
  complex(C_FLOAT_COMPLEX), dimension(:,:), allocatable  :: exp_pi_i_r2_lamz0, exp_pi_i_r2_lam_ddel, exp_C
  complex(C_DOUBLE_COMPLEX), pointer                     :: pfftin(:,:), pfftout(:,:)
  type(C_PTR) :: p1, p2   ! type of fftw
  type(C_PTR) :: plan, plan1 !type of fftw
  real,             dimension(:,:), allocatable  :: x_coord, y_coord, Ar_coord, rr_coord, xfreq, yfreq, Or_coord
  double precision, dimension(:,:), allocatable  :: psiA_00
  real,    dimension(:,:), allocatable  :: ApA, ApO, ApC
  real,    dimension(:,:), allocatable  :: IO_r, IC_r, ID_r
  real,    dimension(:,:), allocatable  :: temp
  real(kind=8), dimension(:,:,:), allocatable  :: temp83
  real,    dimension(:),   allocatable :: rho_arr, dS_arr, alpha_arr, beta_arr                       ! 1D arrays with polar coordinates, area and cartesian coordinates  
  integer                 :: n,i, j, k, x, y, Rpts, left, right
  integer                 :: NN, MM, NUM_pts, pts_proc, NMmin, NMmax      ! NN, MM - number of points; NUM_pts - expeted/real number of points; pts_proc - points/process, NMmin-NMmax - proints from the whole set for the current process
  integer                 :: sampl_type                                   ! 1 - polar, 2 - cartesian
  real                    :: dx, alpha, beta, rho, t
  character(len=128)      :: filenameIO, filenameID, temps
  real                    :: t_start, t_finish, t1_start , t_bruno_start, t_bruno_finish
  integer                 :: imax, imin                                              ! imin, imax - indexes for interpolation of psiA_00 to psiA
  integer                 :: ierr

  call initialization

    
  sampl_type=1                                                            ! 1 - polar, 2 - cartesian
  NN=50                                                                   ! spatial sampling: cartes - number of points (aim), polar - sampling in rho
  MM=2 
  NUM_pts=NN*MM                                                           ! here to initialize arrays; after - real number of points

  allocate(x_coord(JJ,JJ),y_coord(JJ,JJ),Ar_coord(JJ,JJ),rr_coord(JJ,JJ))
  allocate(xfreq(JJ,JJ),yfreq(JJ,JJ),Or_coord(JJ,JJ))
  allocate(IO_r(JJ,JJ),IC_r(JJ,JJ),ID_r(JJ,JJ))
  allocate(ApA(JJ,JJ),ApO(JJ,JJ),ApC(JJ,JJ))
  allocate(temp(JJ,JJ),temp83(JJ,JJ,2))
  allocate(exp_pi_i_r2_lamz0(JJ,JJ),exp_pi_i_r2_lam_ddel(JJ,JJ),exp_C(JJ,JJ))
  allocate(rho_arr(NUM_pts),dS_arr(NUM_pts),alpha_arr(NUM_pts),beta_arr(NUM_pts))

  !i=fftw_init_threads()  
  !write(*,'(A, I0, A)') "fftw_init_threads returns ", i, "  (0 = errors)"
  p1=fftw_alloc_complex(int(JJ*JJ, C_SIZE_T))                                                     !
  p2=fftw_alloc_complex(int(JJ*JJ, C_SIZE_T))                                                     ! FFTW-allocated arrays
  call c_f_pointer(p1, pfftin, [JJ,JJ])                                                          ! 
  call c_f_pointer(p2, pfftout, [JJ,JJ])                                                         !
  j = 0
  !call fftw_plan_with_nthreads(j)
  write(*,'(A,I0,A)') "FFTW uses ", j, " threads"
  ierr=fftw_import_wisdom_from_filename('diff_cartes.fftw_wisdom.FFTW_PATIENT.wisdom')
  write(*,'(A, I0)') "Loaded wisdom: 'diff_cartes.fft_wisdom.FFTW_PATIENT.widsom'; ierr=", ierr
  call cpu_time(t_start)
  write(*,'(A)',advance='no') "Calculating FFTW plan .... "
  plan1 =fftw_plan_dft_2d(JJ,JJ,pfftin,pfftout,FFTW_FORWARD,FFTW_MEASURE)                        ! FFTW-allocated arrays; wisdoms:  FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE, FFTW_WISDOM_ONLY
  call cpu_time(t_finish)
  write(*,'(A, F6.1, A)') "   took ", t_finish-t_start, " sec"
  if ( ((t_finish-t_start) > 100.0) ) ierr=fftw_export_wisdom_to_filename('diff_cartes.fftw_wisdom.FFTW_MEASURE.wisdom')

  write(*,*) " "
  write(*,'(a, 3(a,f0.3, a))') "ASPIICS parameters: ", "z0=", z0, " mm;", "   f=", f, " mm;", "  z1=", z1, " mm;"
  write(*,'(A, A,F8.1,A,F7.3,A,F6.2,A)')    "                    ", "Ap=", rA, "   mm; rIO=", rIO, " mm;  rC=", rC, "  mm;"
  write(*,'(A,I0,A)')          "Size of arrays A, O', C, D --- JJ=", JJ, " points"
  write(*,'(A,F4.1,A)')        "Linear size of Aperture plane 2*LA2=", LA2*2, " mm"
  write(*,'(A,F6.4,A)')        " --> spactial scale in B,D  planes ", JJ/(JJ-1)/2.0/LA2*lam*f*1000, " mkm/pix"
  write(*,'(A,F6.4,A)')        " --> spactial scale in O'   plane  ", JJ/(JJ-1)/2.0/LA2*lam*z1*1000, " mkm/pix"
  write(*,*) " "

  call read_psi_fits(filename=filename_psiA,output=psiA_00)                                       ! PsiA_00 (1:Rpts,1:3)
  Rpts=size(psiA_00,1)
  write(*,'(A, I0, A, I0, A)', advance='no') "Psi_A input array ", Rpts, "x", size(psiA_00,2), " points"
  !write(*,'(A,F6.2,A,F6.2,A,F9.6,A)',advance='no') " from ", psiA_00(1,1), " to ", psiA_00(Rpts,1), "mm, dx=", (psiA_00(Rpts,1)-psiA_00(1,1))/(Rpts-1)*1000, " mkm"
  write(*,'(A, F9.6, A)') " (", (psiA_00(2,1)-psiA_00(1,1))*1000, " mkm)"
  write(*,*) " "

  dx=2*LA2/(JJ-1)                                                                                 ! dx in Aperture plane
  do i=1,JJ
    x_coord(i,:)=-LA2+dx*(i-1)
    y_coord(:,i)=-LA2+dx*(i-1)
    xfreq(i,:) = (-real(JJ)/2 + (i-1)*JJ/(JJ-1))/2.0/LA2
    yfreq(:,i) = (-real(JJ)/2 + (i-1)*JJ/(JJ-1))/2.0/LA2
  end do
  Ar_coord = sqrt(x_coord**2 + y_coord**2)
  Or_coord = sqrt( (xfreq*lam*z1)**2 + (yfreq*lam*z1)**2 )
  ApA = Ar_coord
  where (ApA .le. rA)
    ApA = 1.0
  elsewhere
    ApA = 0.0
  endwhere
  ApO = 1.0
  where ((Or_coord > 0.489d0) .and. (Or_coord <= rIO)) ApO = Or_coord*0.0
  ApC = Ar_coord  
  where (ApC .le. rC)
    ApC = 1.0
  elsewhere
    ApC = 0.0
  endwhere 
  !exp_pi_i_r2_lamz0=exp(-dpi*cmplx(0,1)*Ar_coord**2.0/lam/z0)				! the case z1=z1_0, i.e. IO is in O' 
  !exp_pi_i_r2_lam_ddel=exp(dpi*cmplx(0,1)*Ar_coord**2.0/lam/(z1*(2*z1-f)/(z1-f)))      !
  exp_pi_i_r2_lamz0   =exp(dpi*cmplx(0,1)*Ar_coord**2.0/lam*(f-z1)/z1/f)                    ! IO is out of O'; using z1 for O'
  exp_pi_i_r2_lam_ddel=exp(dpi*cmplx(0,1)*Ar_coord**2.0/lam*(z1_0-f)/z1_0/(2*z1_0-f))       ! and z1_0 for C, D planes ---- but here should be (z1-f)/z1/f
  exp_C               =exp(dpi*cmplx(0,1)*Or_coord**2.0/lam*(z1_0-z1)/z1/z1_0)              ! goes into C phase; =1 for z1=z1_0 
 
  call cpu_time(t_start)
  t1_start=t_start 
  if (sampl_type==1) then
    call solar_sampling_polar(NN,MM,NUM_pts,alpha_arr,beta_arr,dS_arr,rho_arr)        ! NN - sampling in rho, MM - sampling in phi
    write(*,'(A,I0,A,I0,A,I0,A)') "Using polar sampling; NN=", NN, " (in rho); MM=", MM, " (in phi). Total: ", NUM_pts, " points"
  else
    call solar_sampling_cartes(NN,NUM_pts,alpha_arr,beta_arr,dS_arr,rho_arr)          ! NN - aimed numer of points, NUM_pts - obtained number
    write(*,'(A,I0,A,I0,A,F9.2,A)') "Using cartesian sampling. Aim: NN=", NN, " points. Obtained NUM_pts=", NUM_pts, ", i.e.", &
                                     RSun/sqrt(dS_arr(1)), " points in rho"
  end if  
  alpha_arr=alpha_arr + dalpha                    ! and here we can shift the Sun
  
  deallocate(Ar_coord,Or_coord,xfreq,yfreq)       ! they are not used in the next part; to preserve memory

  !  *****************************************************PROBLEEEEMMMMMMMM HEEREEEE*********************************************************
  !  ****************** main loop with FFTs over all solar points corresponding to this process ******************************************
  do j=1,NUM_pts
    rho=rho_arr(j)
    if (j==1) then
      write(*,*) " "
      write(*,'(A,I0,A,F7.5,A,F4.2)') "j=", j, "  rho=", rho*1000.0, "/4.6542 mrad;  B_Sun=", BSun(rho)
    end if
	
	
	
    alpha=alpha_arr(j) 
    beta=beta_arr(j)   
    if (j==1) write(*,*) "... checking typical performance for j=1 ..."
    if (j==1) write(*,'(A)',advance='no') "         Interpolating psiA to the Aperture plane "
    rr_coord = sqrt( (x_coord + alpha*z0)**2.0 + (y_coord + beta*z0)**2.0 )                      ! For the inclined beam (rho,phi), Aperture has new coordinates
    if (j==1) write(*,'(A, F6.2, A, F6.2, A)', advance='no') "(", minval(rr_coord), ":", maxval(rr_coord), " mm) ..."
        
		
		
		
    !!! **********************  version7 - psiA(*,2:3) represents complex, interpolating via real ******************
    call r8vec_bracket(Rpts, psiA_00(:,1), real(minval(rr_coord),kind=8), imin, i)               ! Interpolating psiA_00 to rr_coord  
    call r8vec_bracket(Rpts, psiA_00(:,1), real(maxval(rr_coord),kind=8), i, imax)               ! omitting unused points in psiA_00  by finding apropriate imin and imax  - r8vec_bracket ( n, x, xval, left, right )
    !if (j==1) write(*,'(A, I0, A, I0, A, F5.2, A, F5.2, A)', advance='yes') "       Searching in psiA_00 inside indexes [", imin, ":", imax, "] with coords (", psiA_00(imin,1), ":", psiA_00(imax,1), ")"
    call interp_linear(2, imax-imin+1, psiA_00(imin:imax,1), psiA_00(imin:imax,2:3), JJ**2, real(rr_coord,kind=8), temp83)
    pfftin = cmplx(temp83(:,:,1),temp83(:,:,2),kind=8)                                           ! **required**  temp83(JJ,JJ,2)) ;;   psiA -> pfftin        
    !!! *************************************************************************************************
    
    call cpu_time(t_finish)
    if (j==1) write(*,'(A, F6.1, A)') "   took ", t_finish-t1_start, " sec" 
	
	
		! ****************   TIME START ********************************************
	
    t1_start=t_finish
    if (j==1) write(*,'(A)',advance='no') "         Performing FFTW on psiO, psiC, psiD                         ..."
	
	
    !  ****************** calculating psiA --- A plane (psiA -> pfftin) ******************************************
    
    call cpu_time(t_bruno_start)
		!$acc parallel loop    
		do i = 1, JJ  ! 0.4s in multicore aceleration of 2.5X
		  !$acc loop 
		  do n = 1, JJ
		   pfftin(n,i)= pfftin(n,i)*exp(-2*dpi*cmplx(0,1)*(alpha*x_coord(n,i)+beta*y_coord(n,i))/lam) * &
						   exp(-dpi*cmplx(0,1)*(alpha**2+beta**2)*z0/lam)
			end do
		end do
		if (j==1) write(*,*) pfftin(1,1)
    call cpu_time(t_bruno_finish)
    if (j==1) write(*,'(A, F6.3, A)') "   took:Bruno 1", t_bruno_finish-t_bruno_start, " sec"
	
	  
    
    call cpu_time(t_bruno_start)
		!$acc parallel loop    
		do i = 1, JJ  ! 0.1s in multicore aceleration of 3X
			!$acc loop 
			do n = 1, JJ
				pfftin(n,i)= pfftin(n,i)*exp(-2*dpi*cmplx(0,1)*tilty*x_coord(n,i)/lam)
			end do
		end do
    call cpu_time(t_bruno_finish)
    if (j==1) write(*,'(A, F6.3, A)') "   took:Bruno 2", t_bruno_finish-t_bruno_start, " sec"
    !  ****************** calculating psiO' --- O' plane (psiO -> pfftin) ******************************************
    
	call cpu_time(t_bruno_start)
		  !$acc parallel loop    
		  do i = 1, JJ  ! 0.1s in multicore aceleration of 3X
			  !$acc loop 
			  do n = 1, JJ
				pfftin(n,i) = pfftin(n,i)*ApA(n,i)*exp_pi_i_r2_lamz0(n,i)   !0.1s
			  end do
		  end do
    call cpu_time(t_bruno_finish)
    if (j==1) write(*,'(A, F6.3, A)') "   took:Bruno 3", t_bruno_finish-t_bruno_start, " sec"
    
    
	  !pfftin = pfftin*ApA*exp_pi_i_r2_lamz0   !0.1s                                                ! *exp(-!dpi*i*rr_coor^2/lamz0)  - to propagate through a PO
    
	
	call cpu_time(t_bruno_start)
		call fftw_execute_dft(plan1, pfftin, pfftout) !0.2s
	call cpu_time(t_bruno_finish)
	if (j==1) write(*,'(A, F6.3, A)') "   took:Bruno 4", t_bruno_finish-t_bruno_start, " sec"
	
	!pfftin = pfftout/real(JJ)   !0.2s   
	call cpu_time(t_bruno_start)
		!$acc parallel loop    
		do i = 1, JJ  
			!$acc loop 
			do n = 1, JJ
				pfftin(n,i) = pfftout(n,i)/real(JJ)  
			end do
		end do
	call cpu_time(t_bruno_finish)
	if (j==1) write(*,'(A, F6.3, A)') "   took:Bruno 5", t_bruno_finish-t_bruno_start, " sec"
   	
	call cpu_time(t_bruno_start)
		pfftin = cshift(pfftin, JJ/2, 1) !0.2s shiftar nas colunas
	call cpu_time(t_bruno_finish)
	if (j==1) write(*,'(A, F6.3, A)') "   took:Bruno 6", t_bruno_finish-t_bruno_start, " sec"
	
	
	call cpu_time(t_bruno_start)
		pfftin = cshift(pfftin, JJ/2, 2) !1.0s shiftar nas rows, Is it possible to parallelize a shift?
	call cpu_time(t_bruno_finish)
	if (j==1) write(*,'(A, F6.3, A)') "   took:Bruno 7", t_bruno_finish-t_bruno_start, " sec"
    
    
    call cpu_time(t_bruno_start)
    !$acc parallel loop    
    do i = 1, JJ  ! *5 acelleration
      !$acc loop 
      do n = 1, JJ !colloumn
            IO_r(n,i) = IO_r(n,i) + real(abs(pfftin(n,i))**2.0*BSun(rho)*dS_arr(j),kind=4) !0.1s 
  	  end do
  	end do
    call cpu_time(t_bruno_finish)
    if (j==1) write(*,'(A, F6.3, A)') "   took:Bruno 8", t_bruno_finish-t_bruno_start, " sec"
	  
	
    !  ****************** calculating psiC --- C plane (psiC -> pfftin) ******************************************
	!pfftin = pfftin*ApO
	call cpu_time(t_bruno_start)
		!$acc parallel loop    
		do i = 1, JJ  
			!$acc loop 
			do n = 1, JJ
				pfftin(n,i) = pfftin(n,i)/ApO(n,i)  
			end do
		end do
	call cpu_time(t_bruno_finish)
	if (j==1) write(*,'(A, F6.3, A)') "   took:Bruno 9", t_bruno_finish-t_bruno_start, " sec"
	
	
    if ( abs(z1-z1_0) > 0.001 ) pfftin = pfftin*exp_C !0s
		
    
    call cpu_time(t_bruno_start)
    	call fftw_execute_dft(plan1, pfftin, pfftout) !0.2s
    call cpu_time(t_bruno_finish)
    if (j==1) write(*,'(A, F6.3, A)') "   took:Bruno 10", t_bruno_finish-t_bruno_start, " sec"
	
	!pfftin = pfftout/real(JJ) !0.2s
	call cpu_time(t_bruno_start)
		!$acc parallel loop    
		do i = 1, JJ  
			!$acc loop 
			do n = 1, JJ
				pfftin(n,i) = pfftout(n,i)/real(JJ)  
			end do
		end do
	call cpu_time(t_bruno_finish)
	if (j==1) write(*,'(A, F6.3, A)') "   took:Bruno 11", t_bruno_finish-t_bruno_start, " sec"
    
    !IC_r = IC_r + real(abs(pfftin)**2.0*BSun(rho)*dS_arr(j),kind=4)   !0.4s
    call cpu_time(t_bruno_start)
		!$acc parallel loop    
		do i = 1, JJ  
			!$acc loop 
			do n = 1, JJ
				IC_r(n,i) = IC_r(n,i)+ real(abs(pfftin(n,i))**2.0*BSun(rho)*dS_arr(j),kind=4)
			end do
		end do
    call cpu_time(t_bruno_finish)
    if (j==1) write(*,'(A, F6.3, A)') "   took:Bruno 12", t_bruno_finish-t_bruno_start, " sec"
	
	
    !  ****************** calculating psiD --- D plane (psiD -> pfftin) ******************************************
	!pfftin = pfftin*ApC*exp_pi_i_r2_lam_ddel !0.1s
    call cpu_time(t_bruno_start)
		!$acc parallel loop    
		do i = 1, JJ  
			!$acc loop 
			do n = 1, JJ
				pfftin(n,i) = pfftin(n,i)*ApC(n,i)*exp_pi_i_r2_lam_ddel(n,i)
			end do
		end do
	call cpu_time(t_bruno_finish)
    if (j==1) write(*,'(A, F6.3, A)') "   took:Bruno 13", t_bruno_finish-t_bruno_start, " sec"
	
    call cpu_time(t_bruno_start)
    	call fftw_execute_dft(plan1, pfftin, pfftout) !0.2s
    call cpu_time(t_bruno_finish)
    if (j==1) write(*,'(A, F6.3, A)') "   took:Bruno 14", t_bruno_finish-t_bruno_start, " sec"
    
	!pfftin = pfftout/real(JJ) !0.2s
	call cpu_time(t_bruno_start)
		!$acc parallel loop    
		do i = 1, JJ  
			!$acc loop 
			do n = 1, JJ
				pfftin(n,i) = pfftout(n,i)/real(JJ)
			end do
		end do
	call cpu_time(t_bruno_finish)
	if (j==1) write(*,'(A, F6.3, A)') "   took:Bruno 15", t_bruno_finish-t_bruno_start, " sec"

  
    call cpu_time(t_bruno_start)
		!$acc parallel loop    
		do i = 1, JJ  ! in multicore aceleration of 5X , rows
		  !$acc loop 
		  do n = 1, JJ !colloumn
				ID_r(n,i) = ID_r(n,i) + real(abs(pfftin(n,i))**2.0*BSun(rho)*dS_arr(j),kind=4) !0.1s 
		  end do
		end do
    call cpu_time(t_bruno_finish)
    if (j==1) write(*,'(A, F6.3, A)') "   took:Bruno 16", t_bruno_finish-t_bruno_start, " sec"
  
	
	!!!! Time total 1+0.3+0.1+0.2+0.2+0.2+1+0.5+0.2+0.2+0.4+0.1+0.2+0.2+0.5 = 5.3s
	
	! ****************   TIME FINISH ********************************************
    call cpu_time(t_finish)
    if (j==1) write(*,'(A, F6.3, A)') "   Final ", t_finish-t1_start, " sec"
    if (j==1) write(*,'(A)',advance='no') " j="
    write(*,'(I0,A)',advance='no') j, "  "
    t1_start=t_finish 
	
	
	
  end do

  IO_r = IO_r/Ap_size
  IC_r = IC_r/Ap_size
  ID_r = ID_r/Ap_size
  
  call cpu_time(t_finish)
  write(*,'(A, F6.1, A)') "Interpolation and FFT took ", t_finish-t_start, " sec" 
  t1_start=t_finish

  call fftw_destroy_plan(plan1)  

  !call write_it_2_fits(real(IC_r,kind=4),filename='IC_r.fits')                    !aimag(psiO)
  temps=trim(generate_filebasename(NN,MM,NUM_pts,sampl_type))
  filenameIO = 'IO_r'//temps
  call write_it_2_fits(real(IO_r,kind=4),trim(filenameIO),NN,MM,NUM_pts,sampl_type)
  filenameID = 'ID_r'//temps
  call write_it_2_fits(real(ID_r,kind=4),trim(filenameID),NN,MM,NUM_pts,sampl_type)
 
  call cpu_time(t_finish)
  t1_start=t_finish
  write(*,'(A, F9.3, A)') "Interpolation, FFT and saving took ", t_finish-t_start, " sec"
 
  deallocate(psiA_00)
  deallocate(x_coord, y_coord, rr_coord)
  deallocate(IO_r,IC_r,ID_r)
  deallocate(ApA,ApO,ApC)
  deallocate(temp,temp83) 
  deallocate(exp_pi_i_r2_lamz0,exp_pi_i_r2_lam_ddel,exp_C)
  deallocate(rho_arr,dS_arr,alpha_arr,beta_arr)
  call fftw_free(p1)
  call fftw_free(p2)
  !call fftw_cleanup_threads()

contains 
  subroutine solar_sampling_polar(NN,MM,NUM_pts,alpha,beta,dS,rho)
    implicit none
    integer, intent(in)         :: NN,MM,NUM_pts
    real, dimension(NUM_pts), intent(out) :: alpha,beta,dS,rho
    real, dimension(NN,MM)                :: rho_arr2, phi_arr2, dS_arr2
    real, dimension(NUM_pts)              :: rho_arr, phi_arr
    real          :: drho, dphi
    integer       :: i

    drho=RSun/NN
    dphi=2*dpi/MM
    do i=1,NN
      rho_arr2(i,:)=drho*(i-1+0.5)
      dS_arr2(i,:) =drho*(i-1+0.5)*drho*dphi
    end do
    do i=1,MM
      phi_arr2(:,i)=dphi*(i-1)
    end do
    rho_arr=reshape(rho_arr2,(/NN*MM/))
    phi_arr=reshape(phi_arr2,(/NN*MM/))
    dS=reshape(dS_arr2,(/NN*MM/))
    alpha=rho_arr*cos(phi_arr)
    beta=rho_arr*sin(phi_arr)
    rho=rho_arr
  end subroutine solar_sampling_polar

  subroutine solar_sampling_cartes(NN,MM,alpha,beta,dS,rho)
    implicit none
    integer, intent(in)        :: NN       ! expected number of points
    integer, intent(inout)     :: MM       ! array size/real number of points
    real, dimension(MM), intent(out)   :: alpha,beta,dS,rho
    real        :: dx
    integer     :: i, j                

    dx=RSun * sqrt(dpi/NN)
    dS = dx**2
    alpha(1)=0.0
    beta(1)=0.0
    rho(1)=0.0
    MM=1
    i=1
    cyclex: do
      j=0
      cycley: do
        alpha(MM+1:MM+4)=(/ i*dx, -j*dx, -i*dx, j*dx /)
        beta(MM+1:MM+4) =(/ j*dx,  i*dx, -j*dx,-i*dx /)
        rho(MM+1:MM+4)  =sqrt((i*dx)**2+(j*dx)**2) 
        MM = MM+4
        j = j+1
        if ( sqrt((i*dx)**2 + (j*dx)**2) > RSun ) exit cycley
      end do cycley
      i = i+1
      if ( i*dx > RSun ) exit cyclex
    end do cyclex
  end subroutine solar_sampling_cartes

  function BSun(rho)
    real, intent(in) :: rho
    real      :: Bsun
    if ( rho < RSun ) then
      BSun=1.0-0.762*(1-sqrt(1-(rho/RSun)**2))-0.232*(1-(rho/RSun)**2)*log10(sqrt(1-(rho/RSun)**2)) !  BSun=(rho_axis lt RSun)  
    else
      BSun=0.0
    end if
  end function
 
  function generate_filebasename(NN,MM,NUM_pts,sampl_type)
    use parameters
    character(len=128)   :: generate_filebasename
    character(len=128)   :: tmp1
    integer,intent(in)   :: NN,MM,NUM_pts,sampl_type
    character(len=64)    :: Tstr,Sstr,NNstr

    if ( tilty >= 0.0 ) then
      write(Tstr,'(A,I2.2)') ".T", int(180.0/dpi*3600.0*tilty)    !tilty = 0.0/3600.0*dpi/180.0  
    else 
      write(Tstr,'(A,I3.2)') ".T", int(180.0/dpi*3600.0*tilty) 
    endif
    if ( dalpha >= 0.0 ) then
      write(Sstr,'(A,I2.2)') "_S", int(180.0/dpi*3600.0*dalpha)   !dalpha = 0.0/3600.0*dpi/180.0        
    else   
      write(Sstr,'(A,I3.2)') "_S", int(180.0/dpi*3600.0*dalpha) 
    endif
    if ( sampl_type==1 ) then
      write(NNstr,'(A,I0,A,I0)') ".rho",NN,".phi",MM
    else
      write(NNstr,'(A,I0)') ".NN", NN
    endif
    !        (A,   I6,        A,     A,      I04,          A,  
    !         A,     I4,         A,     I2        A,    I0,  A,            A            A,         A)
    write(tmp1,'(A,I6,A,A,I04,A,A,I4,A,I2,A,I0,A,A,A,A)') & 
            ".z", int(z0), trim(NNstr), ".IO", int(rIO*1000), trim(IOstr), & 
            ".C", int(rC*100), ".LA2", int(LA2), ".JJ", JJ, trim(Tstr), trim(Sstr), trim(ext), ".fits"

    generate_filebasename = tmp1 
  end function
 
  subroutine write_it_2_fits(input,filename,NN,MM,NUM_pts,sampl_type)
    use parameters
    real, dimension(:,:), intent(in) :: input
    integer, dimension(2)            :: naxes
    integer                      :: iunit, naxis, status, bitpix, nelements, group, fpixel, blocksize
    logical                      :: simple, extend
    !character(len=24)           :: filename
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: NN, MM, NUM_pts, sampl_type
    character(len=128)           :: tmp
    double precision             :: scaleB, scaleO, dtemp

    !filename = 'my_test_fits2.fits'
    status   = 0
    naxis    = 2
    naxes(1) = JJ
    naxes(2) = JJ
    bitpix   = -32        ! 32-bit floating point
    simple   = .true.
    extend   = .false.
    group    = 1          ! some deprecated parameter
    fpixel   = 1
    nelements = JJ*JJ
    scaleB = JJ/(JJ-1)/2.0/LA2*lam*f*1000
    scaleO= JJ/(JJ-1)/2.0/LA2*lam*z1*1000

    call ftopen(iunit, filename, 1, blocksize, status)
    if (status .eq. 0)then
       call ftdelt(iunit, status)
       write(*,'(A,A)') "%WRITE_IT_2_FITS: Deleting old ", filename
    else
       call ftcmsg
    end if

    status = 0
    write(*,'(A,A)') "%WRITE_IT_2_FITS: Saving ", filename
    call ftgiou(iunit, status)
    if (status .gt. 0) write(*,'(A,I6)') "  Could not get iunit", status
    call ftinit(iunit, filename, 1, status)
    if (status .gt. 0) write(*,'(A,I6)') "  Could not init fits", status
    call ftphpr(iunit, simple, bitpix, naxis, naxes, 0, 1, extend, status)
    if (status .gt. 0) write(*,'(A)') "  Could not write fits header"
    call ftppre(iunit, group, fpixel, nelements, input, status )
    if (status .gt. 0) write(*,'(A)') "  Could not write fits file"
    call ftpkyg(iunit, 'Z0', z0, 1, "[mm] ISD distance", status)
    call ftpkyg(iunit, 'F',  f,  3, "[mm] Primary objective focal length", status)
    call ftpkyg(iunit, 'Z1', z1, 3, "[mm] PO-IO distance", status)
    call ftpkyd(iunit, 'LAMBDA', lam, 1, "[mm] Lambda", status)
    call ftpkyg(iunit, 'R_A',   rA, 3, "[mm] Entrance aperture radius", status)
    call ftpkyg(iunit, 'R_IO', rIO, 3, "[mm] Internal occulter radius", status)
    call ftpkyg(iunit, 'R_C',   rC, 3, "[mm] Lyot stop radius", status)
    call ftpkyg(iunit, 'D',   d, 3, "[mm] distance plane O' - plane C; f_IO=z1/2", status)
    call ftpkyg(iunit, 'L',   l, 3, "[mm] ??distance plane C  - plane D; ??f_O3=z1/2??", status)
    call ftpkyg(iunit, 'TILT_Y', tilty*180.0/dpi*3600.0, 1, "[arcsec] Tilt of the coronograph WRT Oy", status)
    call ftpkyg(iunit, 'SHIFT_S', dalpha*180.0/dpi*3600.0, 1, "[arcsec] Shift of the Sun along Ox", status)
    call ftpkyj(iunit, 'JJ', JJ, "[pix] Linear size of arrays", status)
    call ftpkyg(iunit, 'LA2', LA2, 1, "[mm] Half-size of entrance array", status)
    call ftpkyg(iunit, 'PSCALE_B', scaleB, 4, "[mkm/pix] Plate scale in B, D planes", status)
    call ftpkyg(iunit, 'PSCALE_O', scaleO, 4, "[mkm/pix] Plate scale in O' plane", status)
    if (sampl_type==1) then
      call ftpkyj(iunit,'NN',NN, '[pts] Polar sampling, points in rho', status)
      call ftpkyj(iunit,'MM',MM, '[pts] Polar sampling, points in phi', status)
    else
      call ftpkyj(iunit, 'NN_aim', NN, '[pts] Cartesian sampling, number of points - as intended', status)
      call ftpkyj(iunit, 'NN_real', NUM_pts, '[pts] Cartesian sampling, number of points', status)
    end if
    if (status .gt. 0) write(*,'(A)') "  Could not write some keyword"
    call ftclos(iunit, status)
    call ftfiou(iunit, status)
  end subroutine write_it_2_fits
 
  subroutine read_fits(filename,output)
    real, dimension(:,:), allocatable, intent(out) :: output
    real, dimension(:,:), allocatable :: temp
    character(len=*), intent(in)      :: filename
    integer, dimension(2)             :: naxes
    integer   :: iunit, naxis, status, nelements, group, readwrite, blocksize, nfound, firstpix, nullval, nbuffer, j
    logical   :: anynull, verbose
    real, dimension(:), allocatable :: buffer

    verbose=.false.
    if (verbose) write(*,'(A,A)') "%READ_FITS: reading ", filename
    status=0
    readwrite=0
    call ftgiou(iunit, status)
    call ftopen(iunit, filename, readwrite, blocksize, status)
    if (status .gt. 0) write(*,'(A,I6)') "Coulnd not open file for reading ", status
    call ftgknj(iunit, 'NAXIS', 1, 2, naxes, nfound, status)
    if (nfound .ne. 2) write(*,'(A,I6)') "Failed to read NAXIS keyword ", status
    if (verbose) write(*,'(A,I0,A,I1,A)') "  Reading file by ", naxes(1), "x", naxes(2), " pixels"
    nelements=naxes(1)*naxes(2)
    group=1
    firstpix=1
    nullval=-999
    nbuffer=naxes(1)
    allocate(buffer(nbuffer))
    allocate(temp(nbuffer,naxes(2)))
    allocate(output(nbuffer,naxes(2)))
    if (verbose) write(*,'(A, I0, A)') "  Allocating buffer with ", nbuffer, " elements"
    do j=1,naxes(2)
      if (verbose) write(*,'(A,I2, A, I0, A, I0)') "  j=", j, ";    elements to read: ", nbuffer, ";   firstpix: ", firstpix
      call ftgpve(iunit, group, firstpix, nbuffer, nullval, buffer, anynull, status)
      nelements=nelements-nbuffer
      firstpix=firstpix+nbuffer
      temp(:,j)=buffer
    end do 
    output=temp
    call ftclos(iunit,status)
    call ftfiou(iunit,status)
    deallocate(temp)
    deallocate(buffer)
    if (verbose) write(*,'(A)') "%READ_FITS: exiting "
  end subroutine read_fits
 
  subroutine read_psi_fits(filename,output)
    double precision, dimension(:,:), allocatable, intent(out) :: output
    real, dimension(:,:), allocatable :: temp
    character(len=*), intent(in)      :: filename
    integer, dimension(2)             :: naxes
    integer   :: iunit, naxis, status, nelements, group, readwrite, blocksize, nfound, firstpix, nullval, nbuffer, j
    logical   :: anynull, verbose
    real, dimension(:), allocatable :: buffer

    verbose=.false.
    if (verbose) write(*,'(A,A)') "%READ_PSI_FITS: reading ", filename
    status=0
    readwrite=0
    call ftgiou(iunit, status)
    call ftopen(iunit, filename, readwrite, blocksize, status)
    if (status .gt. 0) write(*,'(A,I6)') "Coulnd not open file for reading ", status
    call ftgknj(iunit, 'NAXIS', 1, 2, naxes, nfound, status)
    if (nfound .ne. 2) write(*,'(A,I6)') "Failed to read NAXIS keyword ", status
    if (verbose) write(*,'(A,I0,A,I1,A)') "  Reading file by ", naxes(1), "x", naxes(2), " pixels"
    nelements=naxes(1)*naxes(2)
    group=1
    firstpix=1
    nullval=-999
    nbuffer=naxes(1)
    allocate(buffer(nbuffer))
    allocate(temp(nbuffer,3))
    allocate(output(nbuffer,3))
    if (verbose) write(*,'(A, I0, A)') "  Allocating buffer with ", nbuffer, " elements"
    do j=1,3
      if (verbose) write(*,'(A,I2, A, I0, A, I0)') "  j=", j, ";    elements to read: ", nbuffer, ";   firstpix: ", firstpix
      call ftgpve(iunit, group, firstpix, nbuffer, nullval, buffer, anynull, status)
      nelements=nelements-nbuffer
      firstpix=firstpix+nbuffer
      temp(:,j)=buffer
    end do 
    output=temp
    call ftclos(iunit,status)
    call ftfiou(iunit,status)
    deallocate(temp)
    deallocate(buffer)
    if (verbose) write(*,'(A)') "%READ_PSI_FITS: exiting "
  end subroutine read_psi_fits

  function r8vec_ascends_strictly ( n, x )
    implicit none

    integer(kind=4)  :: n, i
    logical          :: r8vec_ascends_strictly
    real(kind=8)     :: x(n)

    do i = 1, n - 1
      if ( x(i+1) <= x(i) ) then
        r8vec_ascends_strictly = .false.
        return
      end if
    end do
    r8vec_ascends_strictly = .true.
    return
  end function 
  
  subroutine r8vec_bracket ( n, x, xval, left, right )
    implicit none

    integer(kind=4)   :: n, i, left, right
    real ( kind = 8 ) :: x(n)
    real ( kind = 8 ) :: xval
    integer(kind=4)   :: l,r,m
    logical           :: enter_do

    l = 1
    r = n
    enter_do = .true.
    !!!!! ********* I really don't understand why this optimization hangs the whole procedure**********
    !!if (xval >= maxval(x)) then 
    !!  l=n-1 !left = n-1    
    !!  r=n   !right= n
    !!  enter_do = .false.
    !!  !write(*,*) "%r8vec_bracket: xval>=maxval(x); xval=", xval, "  maxval(x)=", maxval(x)
    !!  !return
    !!end if
    !!if (xval <= minval(x)) then
    !!  l=1   !left = 1
    !!  r=2   !right = 2
    !!  enter_do = .false.
    !!  !write(*,*) "%r8vec_bracket: xval<=maxval(x); xval=", xval, "  minval(x)=", minval(x)
    !! !return
    !!end if
    do while (r-l > 1)
      !if ( .not. enter_do ) write(*,*) "r8vec_bracket: We should never enter here: l,r = ", l, r
      m = int( (l+r)/2 )
      if ( x(m) > xval ) then
        r = m
      else
        l = m
      end if
    end do 
    left = l
    right = r
    return
  end subroutine r8vec_bracket


  !!! ************* this version of interpolation via real and 2-dimensional->complex works faster then others *****************
  subroutine interp_linear ( m, data_num, x_data, f_data, interp_num, x_interp, f_interp )
    implicit none
    integer ( kind = 4 ) :: data_num, interp_num, m
    integer ( kind = 4 ) :: left, right, interp
    real ( kind = 8 )    :: t
    real ( kind = 8 ), dimension(data_num,m)   :: f_data
    real ( kind = 8 ), dimension(interp_num,m) :: f_interp
    real ( kind = 8 ), dimension(data_num)   :: x_data
    real ( kind = 8 ), dimension(interp_num) :: x_interp

    !if ( .not. r8vec_ascends_strictly ( data_num, x_data ) ) then
    !  write ( *, '(a)' ) ' '
    !  write ( *, '(a)' ) 'INTERP_LINEAR - Fatal error!'
    !  write ( *, '(a)' ) '  Independent variable array T_DATA is not strictly increasing.'
    !  stop 1
    !end if

    do interp = 1, interp_num
      t = x_interp(interp)
      call r8vec_bracket ( data_num, x_data, t, left, right )
      f_interp(interp,1:m) = &
        ( ( x_data(right) - t                ) * f_data(left,1:m)   &
        + (                 t - x_data(left) ) * f_data(right,1:m) ) &
        / ( x_data(right)     - x_data(left) )
    end do
    return
  end subroutine interp_linear

end program diffraction
