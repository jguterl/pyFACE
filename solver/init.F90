module init
    use Header
    use functions
    use cap_functions
    use InventoryHeader
    use io
    use error
    use compute


    implicit none

    integer ::unit_Tramp=200
!  real::K0abs_l(nspc)
contains
 subroutine test
   write(*,*) 'hello test init'
 end subroutine test


    subroutine initialize
        !     ******************************************************************
        !     * initialization of arrays and parameters                        *
        !     *                                                                *
        !     * Author: Roman D. Smirnov                                       *
        !     * E-mail: rsmirnov@ucsd.edu; rosmirnov@yahoo.com                 *
        !     *                                                                *
        !     ******************************************************************

        !      initialization of species parameters

        if (verbose_init) call print_milestone('allocation done')
        call init_misc
        if (verbose_init)call print_milestone('init misc done')
        call init_time
        if (verbose_init) call print_milestone('init time done')
        ! !     call init_seed
         call get_ero_flx
        ! call init_temp
         if (verbose_init) call print_milestone('init temp done')
         call init_options
         if (verbose_init) call print_milestone('init options done')
         call compute_inflx
         if (verbose_init) call print_milestone('init compute inflx done')
         call init_volume_species
        if (verbose_init) call print_milestone('init volume species done')
         call init_source
         if (verbose_init) call print_milestone('init source done')
         call init_boundary
         if (verbose_init) call print_milestone('init boundary done')
         call init_reactions
        if (verbose_init) call print_milestone('init reaction done')

         call compute_onthefly_inventory ! init onthefly inventory
         
         if (verbose_debug) call print_milestone('initialization done')
         is_initialized = .TRUE.
         call print_milestone('initialization done')
    end subroutine initialize

    subroutine init_misc()
        integer k


        ! some material constants


        rhocp=rho*cp

        ! tmp=rand(seed)

        ! counter for filename numbers
        !sfln_voldata=0
        !sfln_srfdata=0
        !sfln_heatdata=0

        iter_solver =0

        ! number fo equations to solve
         if(compute_spc) then
        if (solve_heat_eq) then
            neq=nspc*(ngrd+3)+ngrd+1
        else

            neq=nspc*(ngrd+3)
        endif
        else
        if (solve_heat_eq) then
            neq=ngrd+1
        else
            neq=0
            call face_error('nspc=0 and no heat eq solving:need at least one equation to solve...')
        endif
        endif
        do k=1,nspc
            trace_flux_sum_inflx(k)=0.0
            trace_flux_sum_qflx(k)=0.0
            trace_flux_sum_Gdes_l(k)=0.0
            trace_flux_sum_Gdes_r(k)=0.0
            trace_flux_sum_Q_l(k)=0.0
            trace_flux_sum_Q_r(k)=0.0
            trace_flux_min_Gdes_l(k)=1.d99
            trace_flux_min_Gdes_r(k)=1.d99
            trace_flux_max_Gdes_l(k)=0.0
            trace_flux_max_Gdes_r(k)=0.0
            trace_flux_sig_Gdes_l(k)=0.0
            trace_flux_sig_Gdes_r(k)=0.0
        enddo

        do k=1,nspc
            onthefly_int_dens(k)=0.0
            onthefly_int_dsrf(k)=0.0
            onthefly_net_int_dens(k)=0.0
            onthefly_net_int_dsrf(k)=0.0
            onthefly_int_des(k)=0.0
            onthefly_int_dens(k)=0.0
            onthefly_int_src(k)=0.0
        enddo

        !call open_timedata_files

    end subroutine init_misc


    !      subroutine init_seed()
    !      integer i
    !      character*10 cdate, ctime, czone
    !      call date_and_time (cdate, ctime, czone, val)
    !      seed=0
    !
    !      do i=1,8
    !       seed=seed+abs(val(i))
    !      enddo
    !      end subroutine init_seed

    subroutine init_reactions()
        integer i,j,k,l,m

        !     ------------------------------------------------------------------
        !      initialization of reaction constatnts
        !     ------------------------------------------------------------------
        call init_kbin0
        call init_nuth0
        do k=1,nspc
            do l=1,nspc
                do m=1,nspc
                    ebin (k,l,m)=ebinar(k,l,m)
                    do i=1,ndt
                        do j=0,ngrd
                            kbin(i,j,k,l,m)=kbin0(k,l,m)*exp(-ee*ebin(k,l,m)/(kb*temp(i,j)))
                        enddo
                    enddo
                enddo
                eth  (k,l)=etherm(k,l)
                do i=1,ndt
                    do j=0,ngrd
                        nuth(i,j,k,l)=nuth0(k,l)*exp(-ee*eth(k,l)/(kb*temp(i,j)))
                    enddo
                enddo
            enddo
        enddo
        do i=1,ndt
            do j=0,ngrd
                do k=1,nspc
                    rct(i,j,k)=0.0
                    do l=1,nspc
                        rct (i,j,k)=rct(i,j,k)+nuth(i,j,k,l)*dens(i,j,l)*ctherm(i,j,k,l)
                        do m=1,l
                            rct(i,j,k)=rct(i,j,k)+kbin(i,j,k,l,m)*dens(i,j,l)*dens(i,j,m)*cbinar(i,j,k,l,m)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        if (verbose_init) write(iout,*) "Initialization reaction terms: DONE"
    end subroutine init_reactions

    subroutine init_time()

        !      time step initialization
        dt_face=dt0_face
        dt_face_old=dt_face
        !      dt=cdt*ttm
        time_savevol=start_time
        time_savetime=start_time
        time=start_time
        solver_step_count=0
        if (verbose_init) write(iout,*) "Initialization time parameters: DONE"
    end subroutine init_time


    subroutine get_ero_flx
        integer:: j,ngrd2,k
        real::dx0
        real:: a,b,sa,sm

        !      initialization of grid arrays
        x(0)=0.0
        if (grid_type.eq."U") then
            dx0=length/ngrd
            do j=1,ngrd
                x  (j  )=j*dx0
                dx (j-1)=x(j)-x(j-1)
            enddo
            dx (ngrd)=dx(ngrd-1)
        elseif (grid_type.eq."S") then
            if (grid_gen_mode.eq."alpha") then

                ngrd2=ngrd/2
                if (mod(ngrd,2) .ne. 0) then
                    dx0=length*(alpha-1.d0)/((1.d0+alpha)*alpha**ngrd2-2.d0)
                    x (0)=0.0
                    dx(0)=dx0
                    do j=1,ngrd2
                        x (j)=x(j-1)+dx(j-1)
                        dx(j)=dx(j-1)*alpha
                    enddo
                    do j=ngrd2+1,ngrd
                        x (j)=x(j-1)+dx(j-1)
                        dx(j)=dx(j-1)/alpha
                    enddo
                else
                    dx0=0.5d0*length*(alpha-1.d0)/(alpha**ngrd2-1.d0)
                    x (0)=0.0
                    dx(0)=dx0
                    do j=1,ngrd2-1
                        x (j)=x(j-1)+dx(j-1)
                        dx(j)=dx(j-1)*alpha
                    enddo
                    x (ngrd2)=x (ngrd2-1)+dx(ngrd2-1)
                    dx(ngrd2)=dx(ngrd2-1)
                    do j=ngrd2+1,ngrd
                        x (j)=x(j-1)+dx(j-1)
                        dx(j)=dx(j-1)/alpha
                    enddo
                endif
                dx (ngrd)=dx(ngrd-1)
                if (x(ngrd).ne.length) then
                    call face_error("x(ngrd).ne.length",x(ngrd))
                endif

            else
                call face_error('Grid generation not implemented for this type of grid',grid_type)
            endif
        elseif (grid_type.eq."A") then
            if (grid_gen_mode.eq."seed") then
                if (grid_dx0.ge.length) then
                    call face_error("dx0 cannot be >= than total grid length in this grid mode")
                endif

                a=1.0000001d0;
                b=100d0;
                do while (abs(b-a)>1.d-10)
                    sa=sum_geo(a,ngrd)-length/grid_dx0
                    sm=sum_geo((b+a)/2d0,ngrd)-length/grid_dx0
                    if (sa*sm<0d0) then
                        b=(a+b)/2d0
                    else
                        a=(a+b)/2d0
                    endif
                    if (verbose_init) write(iout,*) "a=",a,"b=",b
                enddo
                alpha=a
                if (alpha.gt.99d0) then
                    alpha=1d0
                endif
                dx0=grid_dx0
                if (verbose_init) write(iout,*) "seed mode: alpha=",alpha," dx0=",dx0
            elseif (grid_gen_mode.eq."alpha") then
                dx0=length*(1.d0-alpha)/(1-alpha**ngrd)
                if (verbose_init) write(iout,*) "alpha mode: alpha=",alpha," dx0=",dx0
            else
                call face_error("unknown grid gen mode")
            endif

            x (0)=0.0
            dx(0)=dx0
            do j=1,ngrd
                x (j)=x(j-1)+dx(j-1)
                dx(j)=dx(j-1)*alpha
            enddo
            dx (ngrd)=dx(ngrd-1)

            if (abs(x(ngrd)-length)>1d-7) then

                call face_error("x(ngrd).ne.length",x(ngrd))

            endif
        else
            call face_error('Unknown type of grid',grid_type)
        endif
        ! find index of depth
        do k=1,nspc
            j_implantation_depth(k)=ngrd
            do j=0,ngrd
                if (x(j).gt.implantation_depth(k)+x(0)) then
                    j_implantation_depth(k)=j-1
                    exit
                endif
            enddo
            if (verbose_init) write(iout,*) " j_implantation_depth(k)=",j_implantation_depth(k)," k=",k
        enddo

        do k=1,nspc
            j_diagnostic_depth(k)=ngrd
            do j=0,ngrd
                if (x(j).gt.diagnostic_depth(k)+x(0)) then
                    j_diagnostic_depth(k)=j-1
                    exit;
                endif
            enddo
            if (verbose_init) write(iout,*) " j_diagnostic_depth(k)=",j_diagnostic_depth(k)," k=",k
        enddo

        call write_grid
        if (verbose_init) write(iout,*) " -- Initialization x grid completed"
    end subroutine get_ero_flx


    subroutine init_src_profile

        real::s
        integer:: j,k

        do k=1,nspc
            do j=0,ngrd
                src_profile(j,k)=0d0

                if (implantation_model(k).eq.'G') then

                    if (implantation_width(k).gt.0d0) then
                        src_profile(j,k)=exp(-0.5d0*abs((x(j)-implantation_depth(k))/implantation_width(k))**2.d0)
                    else
                        call face_error("implantation_width(k)=0 with G implantation model: k=",k)
                    endif
                elseif  (implantation_model(k).eq.'S') then

                    if ((j_implantation_depth(k).gt.0).and.j.le.j_implantation_depth(k)) then
                        src_profile(j,k)=1d0
                    endif


                elseif  (implantation_model(k).eq.'E') then
                    if (implantation_width(k).gt.0d0) then
                        src_profile(j,k)=1.d0-erf((x(j)-implantation_depth(k))/(sqrt2*implantation_width(k)))
                    else
                        call face_error("implantation_width(k)=0 with E implantation model: k=",k)
                    endif


                else

                    call face_error("Unknown implantation model : ",implantation_model(k),"; k=",k)
                endif

            enddo
            s=integrale_src_profile(k)
            if (s.le.0d0) then
                s=1d0
            endif

            do j=0,ngrd
                src_profile(j,k)=src_profile(j,k)/s
            enddo

        enddo

    end subroutine init_src_profile

! subroutine init_gamma_pulse
!   integer k
!   do k=1,nspc
! 
!     if(Gamma_in_pulse(k).eq."S") then
!         Gamma_in_pulse_int(k)=Gamma_in_pulse_S
! 
!     elseif (Gamma_in_pulse(k).eq."N") then
!         Gamma_in_pulse_int(k)=Gamma_in_pulse_N
! 
!     elseif (Gamma_in_pulse(k).eq."B") then
!         Gamma_in_pulse_int(k)=Gamma_in_pulse_B
! 
!     elseif (Gamma_in_pulse(k).eq."R") then
!         Gamma_in_pulse_int(k)=Gamma_in_pulse_R
! 
!     else
!         call face_error("Unknown Gamma_in_pulse::",Gamma_in_pulse(k),"at k=",k)
!     endif
! enddo
! end subroutine init_gamma_pulse

    subroutine init_source

        integer::i,j,k,l
        !call init_gamma_pulse
        call init_src_profile

        call compute_inflx
        !      Initialization of sources
        do i=1,ndt
            do j=0,ngrd
                do k=1,nspc
                    srs (i,j,k)=inflx(k)*src_profile(j,k)
                    src (i,j,k)=srs   (i,j,k)*csours(i,j,k)
                    do l=1,nspc
                        srb (i,j,k,l)=srcbin(  j,k,l)
                        src (i,j,k  )=src   (i,j,k  )+srb(i,j,k,l)*dens(i,j,l)*csrbin(i,j,k,l)
                    enddo
                enddo
            enddo
        enddo
        do i=1,ndt
            do k=1,nspc
                do j=1,ngrd-1
                    jout(i,k)=jout(i,k)+srs(i,j,k)*(1.d0-csours(i,j,k))*0.5d0*(dx(j-1)+dx(j))
                    do l=1,nspc
                        jout(i,k)=jout(i,k)+srb(i,j,k,l)*dens(i,j,l)*(1.d0-csrbin(i,j,k,l))*0.5d0*(dx(j-1)+dx(j))
                    enddo
                enddo
                jout(i,k)=jout(i,k)+srs(i,0   ,k)*(1.d0-csours(i,0   ,k))*0.5d0*dx(0   )
                jout(i,k)=jout(i,k)+srs(i,ngrd,k)*(1.d0-csours(i,ngrd,k))*0.5d0*dx(ngrd)
                do l=1,nspc
                    jout(i,k)=jout(i,k)+srb(i,0   ,k,l)*dens(i,0   ,l)*(1.d0-csrbin(i,0   ,k,l))*0.5d0*dx(0   )
                    jout(i,k)=jout(i,k)+srb(i,ngrd,k,l)*dens(i,ngrd,l)*(1.d0-csrbin(i,ngrd,k,l))*0.5d0*dx(ngrd)
                enddo
            enddo
        enddo
        if (verbose_init) write(iout,*) " -- Initialization source terms completed"
    end subroutine

    subroutine init_boundary
        integer ::i,k
        real::tmp, Edes_lc,Edes_rc

        !      surface parameters
        if (verbose_init) write(iout,*) "Initialization boundary variables"




        do k=1,nspc



            if (mass(k)*gas_temp(k) .ne. 0.0) then
                j0(k)=j0(k)+gas_pressure(k)/sqrt(twopi*mass(k)*ee*gas_temp(k))
            endif


            if (active_cap_surface) then
                tmp=2d0*(dsrfl0(k)/dsrfm(k)*2d0-1.d0)
                Edes_lc=Edes_lsat(k)+(Edes_l(k)-Edes_lsat(k))*0.5d0*(1.d0-erf(tmp))
                tmp=2d0*(dsrfr0(k)/dsrfm(k)*2d0-1.d0)
                Edes_rc=Edes_rsat(k)+(Edes_r(k)-Edes_rsat(k))*0.5d0*(1.d0-erf(tmp))
            else
                Edes_lc=Edes_l(k)
                Edes_rc=Edes_r(k)
            endif
            ! left
            if (left_surface_model_int(k).eq.surf_model_S) then
                K0abs_l(k)=1.d0
                K0des_l(k)=nu(k)*lambda**(2*order_desorption_left(k)-2)*csrf

                K0b_l(k)=nu(k)*clng
                K0ads_l(k)=nu(k)*lambda*clng
                Kabs_l(k)=j0(k)*K0abs_l(k)*exp(-  ee*Eabs_l(k) /(kb*temp(ndt,0)))
                Kdes_l(k)=2.d0*K0des_l(k)*exp(-  ee*Edes_lc /(kb*temp(ndt,0)))
                Kb_l(k)=        K0b_l(k)  *exp(-  ee*Eb_l(k)   /(kb*temp(ndt,0)))
                Kads_l(k)=      K0ads_l(k)*exp(-  ee*Eads_l(k) /(kb*temp(ndt,0)))
            elseif (left_surface_model_int(k).eq.surf_model_B) then

                K0abs_l(k)=0d0
                K0des_l(k)=nu(k)*lambda**(3*order_desorption_left(k)-2)*csrf
                if (verbose_surface) then
                    write(iout,*) 'init: K0des_l(k)',K0des_l(k),'nu(k)=',nu(k),'csrf=',csrf,' lambda=',lambda
                endif
                K0b_l(k)=0d0
                K0ads_l(k)=0d0

                Kabs_l(k)=min_rate_surface
                Kdes_l(k)=2.d0*K0des_l(k)*exp(-  ee*Edes_lc /(kb*temp(ndt,0)))
                if (verbose_surface) then
                    write(iout,*) 'init: Kdes_l(k)',Kdes_l(k)
                endif
                Kb_l(k)=0d0
                Kads_l(k)=0d0
            elseif (left_surface_model_int(k).eq.surf_model_N) then
                K0abs_l(k)=0d0
                K0des_l(k)=0d0
                K0b_l(k)=0d0
                K0ads_l(k)=0d0
                Kabs_l(k)=0d0
                Kdes_l(k)=0d0
                Kb_l(k)=0d0
                Kads_l(k)=0d0
            else
                call face_error("Init:Unknown left surface model:",left_surface_model_int(k))
            endif

            !right
            if (right_surface_model_int(k).eq.surf_model_S) then

                K0abs_r(k)=1.d0
                K0des_r(k)=nu(k)*lambda**(2*order_desorption_right(k)-2)*csrf
                K0b_r(k)=nu(k)*clng
                K0ads_r(k)=nu(k)*lambda*clng
                Kabs_r(k)=j0(k)*K0abs_r(k)*exp(-  ee*Eabs_r(k) /(kb*temp(ndt,0)))
                Kdes_r(k)=2.d0*K0des_r(k)*exp(-  ee*Edes_rc /(kb*temp(ndt,0)))
                Kb_r(k)=        K0b_r(k)  *exp(-  ee*Eb_r(k)   /(kb*temp(ndt,0)))
                Kads_r(k)=      K0ads_r(k)*exp(-  ee*Eads_r(k) /(kb*temp(ndt,0)))
            elseif (right_surface_model_int(k).eq.surf_model_B) then
                K0abs_r(k)=0d0
                K0des_r(k)=nu(k)*lambda**(3*order_desorption_right(k)-2)*csrf
                K0b_r(k)=0d0
                K0ads_r(k)=0d0

                Kabs_r(k)=0d0
                Kdes_r(k)=2.d0*K0des_r(k)*exp(-  ee*Edes_rc /(kb*temp(ndt,0)))
                Kb_r(k)= 0d0
                Kads_r(k)=0d0
            elseif (right_surface_model_int(k).eq.surf_model_N) then
                Kabs_r(k)=0d0
                Kdes_r(k)=0d0
                Kb_r(k)=0d0
                Kads_r(k)=0d0
                K0abs_r(k)=0d0
                K0des_r(k)=0d0
                K0b_r(k)=0d0
                K0ads_r(k)=0d0
            else
                call face_error("unknown right surface model:",right_surface_model(k))
            endif

            do i=1,ndt
                Gsrf_l(i,k)=0.0
                Gsrf_r(i,k)=0.0

                if (dsrfl0(k) .eq. 0.0) then
                    dsrfl(i,k)=1.d0
                else
                    dsrfl(i,k)=dsrfl0(k)
                endif
                if (dsrfr0(k) .eq. 0.0) then
                    dsrfr(i,k)=1.d0
                else
                    dsrfr(i,k)=dsrfr0(k)
                endif


                ! left
                if ((left_surface_model_int(k).eq.surf_model_S)) then
                    Gabs_l (i,k)=Kabs_l(k)
                    Gdes_l (i,k)=Kdes_l(k) *dsrfl(i,k)**order_desorption_left(k)
                    Gb_l (i,k)  =Kb_l(k)   *dsrfl(i,k)
                    Gads_l (i,k)=Kads_l(k) *dens(i,0   ,k)
                elseif (left_surface_model_int(k).eq.surf_model_N) then
                    Gabs_l(ndt,k)=0d0                           ! Gabsorp=K(gas)
                    Gdes_l (ndt,k)=0d0          ! Gdesorp=K*ns^2
                    Gb_l (ndt,k)  =dsrfl(ndt,k)               ! Gbulk  =K*ns
                    Gads_l (ndt,k)=0d0
                elseif (left_surface_model_int(k).eq.surf_model_B) then
                    Gabs_l (i,k)=0d0
                    Gdes_l (i,k)=Kdes_l(k) *dens(i,0   ,k)**order_desorption_left(k)
                    Gb_l (i,k)  =dsrfl(ndt,k)
                    Gads_l (i,k)=0d0
                endif


                ! right
                if ((right_surface_model_int(k).eq.surf_model_S)) then
                    Gabs_r (i,k)=Kabs_r(k)
                    Gdes_r (i,k)=Kdes_r(k) *dsrfr(i,k)**order_desorption_right(k)
                    Gb_r (i,k)  =Kb_r(k)   *dsrfr(i,k)
                    Gads_r (i,k)=Kads_r(k) *dens(i,ngrd   ,k)
                elseif(right_surface_model_int(k).eq.surf_model_N) then
                    Gabs_r (i,k)=0d0
                    Gdes_r (i,k)=0d0
                    Gb_r (i,k)  =dsrfr(i,k)
                    Gads_r (i,k)=0d0
                elseif (right_surface_model_int(k).eq.surf_model_B) then
                    Gabs_r (i,k)=0d0
                    Gdes_r (i,k)=Kdes_r(k) *dens(i,ngrd   ,k)**order_desorption_right(k)
                    Gb_r (i,k)  =dsrfr(i,k)
                    Gads_r (i,k)=0d0
                endif


                    call compute_cap_factor_surface(k,i)


                jout(i,k)=jout(i,k)+Gdes_l(i,k)

            enddo
        enddo
        if (verbose_init) write(iout,*) "Initialize boundary: DONE"
    end subroutine init_boundary

!     subroutine read_Tramp_file
!         integer ios,i
! 
!         open(unit=unit_Tramp, file=trim(framp_string), iostat=ios,action='read',status='old')
!         if ( ios /= 0 ) then
!             write(iout,*) 'Opening of temperature ramp file "', framp_string ,'" : FAIL '
!             stop
!         endif
!         if (verbose_init)  write(iout,*) 'Opening of temperature ramp file "', trim(framp_string) ,'" : DONE '
! 
!         read (unit_Tramp, '(i4)') nramp
!         allocate(rtime(nramp))
!         allocate(rtemp(nramp))
!         do i=1,nramp
!             read (20, *) rtime(i), rtemp(i)
!         enddo
!         close(unit_Tramp)
!         write(iout,*) 'Reading of temperature ramp file "', trim(framp_string) ,'" : DONE '
!     end subroutine read_Tramp_file


!     subroutine init_temp
!         integer::i,j,n
!  if (verbose_init)  write(iout,*) 'Initialization temperature'
! !         if (tramp1 .ne. tramp0) then
! !             dtemp=(temp_final-temp_init)/(tramp1-tramp0)
! !         else
! !             dtemp=0.0
! !         endif
! ! 
! !         if (framp_string .ne. 'none') then
! !             framp=framp_readfile
! !             call read_Tramp_file
! !         else
! !             framp=framp_none
! !         endif
! 
!         !if (.not.solve_heat_eq) then
! !          if (verbose_init)  write(iout,*) 'Not solving heat equation'
! !             if (framp.ne. framp_none) then
! ! 
! ! 
! !                 do j=0,ngrd
! !                     do i=1,ndt
! !                         if (rtime(1) .gt. time) then
! !                             temp(i,j)=rtemp(1)
! !                         elseif (rtime(nramp) .gt. time) then
! !                             temp(i,j)=rtemp(nramp)
! !                         else
! !                             do n=2,nramp-1
! !                                 if (rtime(n) .gt. time) then
! !                                     temp(i,j)=rtemp(n-1)+(rtemp(n)-rtemp(n-1))*(time-rtime(n-1))/(rtime(n)-rtime(n-1))
! !                                     exit
! !                                 endif
! !                             enddo
! !                         endif
! !                     enddo ! i
! !                 enddo ! j
! ! 
! !             else
! !             if (verbose_init)  write(iout,*) 'Not reading temperature file'
! ! 
! !                 do j=0,ngrd
! !                     do i=1,ndt
! !                         temp(i,j)=temp_init
! !                     enddo
! !                 enddo
! !             endif
! 
!          !  solving heat equation -> initial linear profile of T in bulk
!         !elseif(solve_heat_eq) then
!            ! do j=0,ngrd
!                 !do i=1,ndt
!                  !   temp(i,j)=temp_init+(temp_final-temp_init)*x(j)/length
!                ! enddo
!            ! enddo
!         !endif
! 
!         if (verbose_init)  write(iout,*) 'Initialization of temperature : DONE '
!     end subroutine init_temp

    subroutine init_volume_species()

        integer:: i,j,k
        do k=1,nspc
            if (verbose_init) write(iout,*) 'initial profile of density : ',gprof(k)
            do j=0,ngrd
                do i=1,ndt
                    flx (i,j,k)=0.0
                    ero_flx (i,j,k)=0.0
                    dif_flx (i,j,k)=0.0
                    cdif(i,j,k)=cdif0(k)*exp(-ee*edif(k)/(kb*temp(i,j)))
                    rate_d (i,j,k)=0.0
                    if (gprof(k) .eq. 'S') then
                        if (x(j)>gxmax(k)) then
                            dens(i,j,k)=0
                        else
                            dens(i,j,k)=dens0(k)
                        endif
                    elseif(gprof(k) .eq. 'F') then
                        dens(i,j,k)=dens0(k)
                    elseif(gprof(k) .eq. 'G') then
                        if (gsigm(k)<=0) then
                            call face_error('gsigm(k)<=0 for gaussian profile')
                        endif

                        dens(i,j,k)=dens0(k) *exp(-0.5d0*abs((x(j)-gxmax(k))/gsigm(k))**2.d0)
                      !  if (verbose_init) write(iout,*)"dens0(k)=",dens0(k),"gxmax(k)=",gxmax(k),"gsigm(k)=",gsigm(k)
                    else
                        call face_error('unknow option for n0_profile: gprof(k)=', gprof(k),' k=',k)
                    endif
                    if (dens(i,j,k) .lt. 1.d5) then
                        dens(i,j,k)=1.d5
                    endif
                    if (isnan(dens(i,j,k))) call face_error("dens(i,j,k) is NAN i=",i,"j=",j,"k=",k)
                enddo
            enddo
        enddo
        do k=1,nspc
            do i=1,ndt
                jout(i,k)=0.0
            enddo
        enddo


    end subroutine init_volume_species

    subroutine init_casename(case_name)
        character(*)::case_name
        casename=trim(case_name)
        if (verbose_init) write(iout,*) 'Casename : ', casename
        casename=trim(casename)
    end subroutine init_casename


    subroutine init_path(path)

        character(*)::path
        character(30)::name
        integer ios,unit_testpath

        name='test.path'
        path_folder=trim(path)
        if (verbose_init) write(iout,*) "Files will be saved in the folder :", path_folder
        call system('mkdir -p '//trim(path_folder))
        path_folder=trim(path_folder)//'/'
        dat_folder=trim(path_folder)//trim(casename)//'_dat'
        call system('mkdir -p '//dat_folder)

        call set_unit(unit_testpath)
        open (unit=unit_testpath, file=trim(path_folder)//name,status='replace', form='unformatted', iostat=ios)
        if (ios.ne.0) then
            call face_error('cannot write in the folder : ', trim(path_folder)//name)
        endif
        close(unit_testpath)
         ! default restart filename
        restart_filename=trim(path_folder)//"dsave.rst"
        final_state_file=trim(path_folder)//trim(casename)//".state"
    end subroutine init_path

    real function sum_geo(q,N)
        real,intent(in):: q
        integer,intent(in) :: N
        if (q==1d0) then
            sum_geo=real(N)
        else
            sum_geo=(1-q**real(N))/(1-q)
        endif
    end function sum_geo

    subroutine init_kbin0
        integer :: kk,ll,mm,n
        do kk=1,nspc
            do ll=1,nspc
                do mm=1,nspc
                    kbin0(kk,ll,mm)=0.0
                enddo
            enddo
        enddo
        if (nspc.ge.3) then
            do n=2,nspc-1,2
                kbin0(1  ,n,1)=-nu(n)*lambda**3*cvlm ! we have introduced variable pre-exponential factor read from the input file> It allows flexible rates
                                                     ! with no restriction to activation energy only.
                kbin0(n  ,n,1)=-nu(n)*lambda**3*cvlm
                kbin0(n+1,n,1)=+nu(n)*lambda**3*cvlm
            enddo
        endif
    end subroutine init_kbin0

    subroutine init_nuth0
        integer :: kk,ll,n
        do kk=1,nspc
            do ll=1,nspc
                nuth0(kk,ll)=0.0
            enddo
        enddo
        if (nspc.ge.3) then
            do n=3,nspc,2
                nuth0(1  ,n)=+nu(n)                 ! we have introduce variable pre-exponential factor read from the input file> It allows flexible rates
                                                    ! with no restriction to activation energy only.
                nuth0(n-1,n)=+nu(n)
                nuth0(n  ,n)=-nu(n)
            enddo
        endif
    end subroutine init_nuth0
    
subroutine init_options
integer k,kk 
character*1 ::tmpstr='U'

            if (nspc<1) then
            nspc = 1
            compute_spc=.false. ! we turn off all the computation with R-D equations but we keep nspc=1 to allow allocation of variable (avoid to modify all the output when nspc=0)
              !  write(iout,*) "ERROR: n_species must be at least equal to 1"
              !  stop
            else
            compute_spc=.true.
            endif


        !check order of the numerical solver
            if (order_solver.ne.1.and.order_solver.ne.2.and.order_solver.ne.5) then
                call face_error("order of the solver must be 1 2 or 5. current order : ", order_solver)
            else
            ndt=order_solver+1
            if (verbose_init) write(iout,*) ' order of numerical order: ', order_solver, ' -> ndt =', ndt
            endif


!             if(Q_in_pulse.eq."S") then
!                 Q_in_pulse_int=Q_in_pulse_S
! 
!             elseif (Q_in_pulse.eq."N") then
!                 Q_in_pulse_int=Q_in_pulse_N
! 
!             elseif (Q_in_pulse.eq."B") then
!                 Q_in_pulse_int=Q_in_pulse_B
! 
!             elseif (Q_in_pulse.eq."R") then
!                 Q_in_pulse_int=Q_in_pulse_R
! 
!             else
!                 call face_error("Unknown Q_in_pulse::",Q_in_pulse)
!             endif
! 
! if (verbose_input) write(iout,*) "checking input:",'T_pulse'
!             if (T_pulse.ne."N".AND.Q_in_pulse.ne."N") then
!                 call face_error("cannot run T pulse and Q pulse simultaneously")
!             endif
! 
!         if (verbose_input)  write(iout,*) 'T_pulse:',T_pulse
!             if(T_pulse.eq."S") then
!                 T_pulse_int=T_pulse_S
! 
!             elseif (T_pulse.eq."N") then
!                 T_pulse_int=T_pulse_N
! 
!             elseif (T_pulse.eq."B") then
!                 T_pulse_int=T_pulse_B
! 
!             elseif (T_pulse.eq."R") then
!                 T_pulse_int=T_pulse_R
! 
!             else
!                 call face_error("Unknown T_pulse:",T_pulse)
!             endif


       if (verbose_input) write(iout,*) "checking input:",'left_surface_model_int'
        
        do k=1,nspc
        if (verbose_input) write(iout,*) "k=",k,"; left_surface_model(k)='",trim(left_surface_model(k)),'"'
          if(left_surface_model(k).eq.'S') then
                left_surface_model_int(k)=surf_model_S
            elseif (left_surface_model(k).eq."N") then
                left_surface_model_int(k)=surf_model_N
            elseif (left_surface_model(k).eq."B") then
                left_surface_model_int(k)=surf_model_B
            else
                call face_error("Unknown left_surface_model:'",trim(left_surface_model(k)),"at k=",k)
            endif
        enddo
        
        if (verbose_input) write(iout,*) "checking input:",'right_surface_model'
        do k=1,nspc
           if (right_surface_model(k).eq."S") then
               right_surface_model_int(k)=surf_model_S
           elseif (right_surface_model(k).eq."N") then
               right_surface_model_int(k)=surf_model_N
           elseif (right_surface_model(k).eq."B") then
               right_surface_model_int(k)=surf_model_B
           else
               call face_error("Unknown right_surface_model:" ,right_surface_model(k),"at k=",k)
           endif
       enddo
        
!         if (verbose_input) write(iout,*) "checking input:",'Gamma_pulse'
!           do k=1,nspc
! 
!             if(Gamma_in_pulse(k).eq."S") then
!                 Gamma_in_pulse_int(k)=Gamma_in_pulse_S
! 
!             elseif (Gamma_in_pulse(k).eq."N") then
!                 Gamma_in_pulse_int(k)=Gamma_in_pulse_N
! 
!             elseif (Gamma_in_pulse(k).eq."B") then
!                 Gamma_in_pulse_int(k)=Gamma_in_pulse_B
! 
!             elseif (Gamma_in_pulse(k).eq."R") then
!                 Gamma_in_pulse_int(k)=Gamma_in_pulse_R
! 
!             else
!                 call face_error("Unknown Gamma_in_pulse::",Gamma_in_pulse(k),"at k=",k)
!             endif
!         enddo
        ! check that name of species are unique
        
           do k=1,nspc
             do kk=1,nspc
                 if (k.ne.kk.and.namespc(k).eq.namespc(kk)) then
                call face_error("two species have the same name: k=",k," and k=",kk," namespc=",namespc(k),&
                " and namespc=",namespc(kk))
            endif
       enddo
       enddo









! check that the reduction factor is >0 and lt 1

if (reduction_factor_dt_spc.le.0d0.or.reduction_factor_dt_spc.gt.1d0) then
call face_error("reduction factor dt spc cannot be <=0 and >1: reduction factor dt spc=",reduction_factor_dt_spc)
endif


if (reduction_factor_dt_heat.le.0d0.or.reduction_factor_dt_heat.gt.1d0) then
call face_error("reduction factor dt heat cannot be <=0 and >1: reduction factor dt heat=",reduction_factor_dt_heat)
endif


if ((order_solver.gt.1).AND. (variable_timestep)) then
call face_error("Cannot have variable timestep when order of numerical time scheme >1 (order_solver>1)")
endif


end subroutine init_options
end module init

subroutine initialize_wrapper
  use init
  call initialize()
end subroutine
  subroutine  allocate
  use io
      call gchange('InventoryHeader',0)
      call gchange('SpeciesHeader',0)
      call gchange('BulkHeader',0)
      call gchange('ReactionRatesHeader',0)
      call gchange('GridHeader',0)
      call gchange('ParticleFluxHeader',0)
      call gchange('SurfaceHeader',0)
      call set_default_strings()
      call set_default_arrays()
      call print_milestone('allocation completed')
  end subroutine
  
  subroutine set_default_strings()
  use header 
  integer k 
  if (nspc.gt.0) then
  do k =1, nspc
  left_surface_model(k) = 'S'
  right_surface_model(k) = 'S'
  implantation_model(k) = 'S'
  gprof(k) = 'F'
  enddo
  endif
  end subroutine set_default_strings
   subroutine set_default_arrays()
   use header 

   end subroutine set_default_arrays