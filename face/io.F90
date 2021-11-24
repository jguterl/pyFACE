module IO
    use header
    use error
    use InventoryHeader

    !use modFACE_allocate
    !use modFACE_compute
    implicit none
    save
    integer,allocatable::unit_timedata(:)
    integer::current_ifile=10
contains
! ** time data
    subroutine open_timedata_files()
        integer:: ios,k
        character*256::filename
        ! open time file
        allocate(unit_timedata(nspc))

        do k=1,nspc
            call set_unit(unit_timedata(k))
            write (filename, '(a,a, i2.2, a)') trim(dat_folder),'/time_', k, '.dat'
            open (unit_timedata(k), file=trim(filename),status='replace', iostat=ios)
            if (ios.ne.0) then
                call face_error('Cannot open file ', trim(filename))
            endif
            write (unit_timedata(k), '(24a22)')&
                'time',&
                'tempL',&
                'tempR',&
                'densL',&
                'densR',&
                'NsrfL',&
                'NsrfR',&
                'Gabs_l',&
                'Gabs_r',&
                'Gdes_l',&
                'Gdes_r',&
                'Gb_l',&
                'Gb_r',&
                'Gads_l',&
                'Gads_r',&
                'Qnty',&
                'Src',&
                'Rct',&
                'rate_d',&
                'inflx',&
                'qrad',&
                'Kdes_l',&
                'Kdes_r',&
                'diag_flx'
        enddo
    end subroutine open_timedata_files

   subroutine close_timedata_files()
        integer:: k
        logical:: test
do k=1,nspc
        inquire(unit_timedata(k),opened=test)
        if (test) then
        close(unit_timedata(k))
        endif
          enddo

         if (allocated(unit_timedata)) deallocate(unit_timedata)
    end subroutine close_timedata_files

subroutine save_timedata
    integer::j,k
    real:: qnty,frmn,rctn,diag_flx
    character*256  myfmt1,myfmt2

    write(myfmt1,*) &
        "('+', ' time=', es12.2e3, ' s; T_l=', es12.3,' K; dt=', es12.2, ' s; number of iterations ', i4)"
    write(myfmt2,*) "(es22.11e3,24es22.11e3)"
    if (verbose_step) write (iout, myfmt1) time, temp(ndt,0), dt_face, iter_solver

    do k=1,nspc
        qnty=0.d0
        frmn=0.d0
        rctn=0.d0
        do j=0,ngrd-1
            qnty=qnty+0.5d0*(dens(ndt,j,k)+dens(ndt,j+1,k))*dx(j)
            frmn=frmn+0.5d0*(src (ndt,j,k)+src (ndt,j+1,k))*dx(j)
            rctn=rctn+0.5d0*(rct (ndt,j,k)+rct (ndt,j+1,k))*dx(j)
        enddo
        diag_flx=dif_flx(ndt,j_diagnostic_depth(k),k)

        write (unit_timedata(k),myfmt2)&
            time, &
            temp (ndt,0),&
            temp (ndt,ngrd),&
            dens (ndt,0,   k),&
            dens (ndt,ngrd,k),&
            dsrfl(ndt,     k),&
            dsrfr(ndt,     k),&
            Gabs_l(ndt,k),&
            Gabs_r(ndt,k),&
            Gdes_l(ndt,k),&
            Gdes_r(ndt,k),&
            Gb_l(ndt,k),&
            Gb_r(ndt,k),&
            Gads_l(ndt,k),&
            Gads_r(ndt,k),&
            qnty,&
            frmn,&
            rctn,&
            rate_d  (ndt,0,k),&
            inflx(k),&
            rad,&
            Kdes_l(k),&
            Kdes_r(k),&
            diag_flx
    enddo

end subroutine save_timedata

! ** volume data

    subroutine save_voldata

        integer j, k,ios,unit_voldata

        character*256::  filename,myfmt1,myfmt2,myfmt3

        !     --- saving snapshot of volume distributions ---

        write(myfmt1,*)    "(a, 1pe18.9e2, a)"
        write(myfmt2,*) "(9a22)"
        write(myfmt3,*) "(i22.4, 8es22.12)"
        call set_unit(unit_voldata)
        do k=1,nspc
        if (dump_vol_append) then

         write (filename, '(a,a, i2.2, a)') trim(dat_folder),'/vol', k, '.dat'
         if (first_voldump) then
         open (unit=unit_voldata, file=trim(filename), access='append',status='replace', iostat=ios)
         first_voldump=.false.
         else
         open (unit=unit_voldata, file=trim(filename), access='append',status='unknown', iostat=ios)
         endif

        else
         write (filename, '(a,a, i2.2, a, i4.4, a)') trim(dat_folder),'/vol', k, '_', sfln_voldata, '.dat'
         open (unit=unit_voldata, file=trim(filename), status='replace', iostat=ios)
        endif


            if (ios.eq.0) then

                write (unit_voldata,myfmt1) 'time=', time, ' s'

                write (unit_voldata,myfmt2) &
                    ' icell', &
                    '  x',&
                    ' dens',&
                    'dif_flx',&
                    'ero_flx',&
                    '   src',&
                    '  rct',&
                    '   rate_d',&
                    '  cdif'

                do j=0,ngrd
                    write (unit_voldata, fmt=myfmt3) &
                        j, &
                        x(j), &
                        dens(ndt,j,k), &
                        dif_flx (ndt,j,k), &
                        ero_flx(ndt,j,k), &
                        src (ndt,j,k), &
                        rct (ndt,j,k), &
                        rate_d (ndt,j,k), &
                        cdif(ndt,j,k)
                enddo

                close (unit_voldata)


            else
                write (iout, '(2a)') ' *** error saving ', trim(filename)

            endif

        enddo
        sfln_voldata=sfln_voldata+1
    end subroutine save_voldata


    subroutine save_heatdata()

        integer j, ios,unit_heatdata
        character*256:: filename,myfmt1,myfmt2,myfmt3

        write(myfmt1,*)    "(a, 1pe18.9e2, a)"
        write(myfmt2,*) "(9a22)"
        write(myfmt3,*) "(i22.4, 5es22.12)"
        call set_unit(unit_heatdata)
        if (dump_vol_append) then
         write (filename, '(a,a)') trim(dat_folder),'/heat.dat'
         open (unit=unit_heatdata, file=trim(filename), access='append',status='unknown', iostat=ios)
        else
         write (filename, '(a,a,i4.4, a)') trim(dat_folder),'/heat_', sfln_heatdata, '.dat'
         open (unit_heatdata, file=trim(filename), status='replace',iostat=ios)
        endif


        if (ios.eq.0) then
            write (unit_heatdata, myfmt1) 'time=', time, ' s'

            write (unit_heatdata, myfmt2)&
                '   icell',&
                ' x',&
                '  temp',&
                '  qflx',&
                '  rate_t',&
                ' ero_qflx'

            do j=0,ngrd
                write (unit_heatdata, myfmt3)&
                    j,&
                    x(j),&
                    temp(ndt,j),&
                    qflx(ndt,j),&
                    rate_t (ndt,j),&
                    ero_qflx(ndt,j)
            enddo

            close (unit_heatdata)
        else
            write (iout, '(2a)') ' *** error saving ', trim(filename)
        endif

        sfln_heatdata=sfln_heatdata+1

    end subroutine save_heatdata

    subroutine save_srfdata
        character*256::filename,myfmt1,myfmt2,myfmt3
        integer k,ios,unit_surfdata
        !     --- saving snapshot of surface parameters
        write(myfmt1,*) "(a, 1pe13.4e2, a)"
        write(myfmt2,*) "(12(a13))"
        write(myfmt3,*) "(i13.2, 11(1pe13.4e2))"
        call set_unit(unit_surfdata)

          if (dump_srf_append) then
         write (filename, '(a,a)') trim(dat_folder),'/srf.dat'
         open (unit=unit_surfdata, file=trim(filename), access='append',status='unknown', iostat=ios)
        else
         write (filename, '(a, a, i4.4, a)') trim(dat_folder),'/srf_', sfln_srfdata, '.dat'
         open (unit=unit_surfdata, file=trim(filename), status='replace', iostat=ios)
        endif

        if (ios.eq.0) then

            write (unit_surfdata, myfmt1) 'time=', time, ' s'
            write (unit_surfdata, myfmt2)&
                'spc#',&
                'NsrfL',&
                'Gabs_l',&
                'Gdes_l',&
                'Gb_l',&
                'Gads_l',&
                'NsrfR',&
                'Gabs_r',&
                'Gdes_r',&
                'Gb_r',&
                'Gads_r',&
                'Jout'

            do k=1,nspc
                write (unit_surfdata, myfmt3) &
                    k, &
                    dsrfl(ndt,k),&
                    Gabs_l  (ndt,k),&
                    Gdes_l  (ndt,k),&
                    Gb_l  (ndt,k),&
                    Gads_l  (ndt,k),&
                    dsrfr(ndt,k),&
                    Gabs_r  (ndt,k),&
                    Gdes_r  (ndt,k),&
                    Gb_r  (ndt,k),&
                    Gads_r  (ndt,k),&
                    jout (ndt,k)
            enddo

            close (unit_surfdata)

        else
            write (iout, '(2a)') ' *** error saving ', trim(filename)
        endif

        sfln_srfdata=sfln_srfdata+1

    end subroutine save_srfdata



    subroutine save
        real::tmp

        ! dump time data
!        if (dump_time_dt.gt.0d0) then
!            tmp=2.d0*abs(time-dump_time_dt*nint(time/dump_time_dt))
!        else
!            tmp=1d99
!        endif
        if (dump_time.and.(time_savetime.ge.dump_time_dt.or.time.le.start_time.or.time.ge.end_time)) then
            call save_timedata
            time_savetime=0d0
        endif

        ! dump space data
!        if (dump_space_dt.gt.0d0) then
!            tmp=2.d0*abs(time-dump_space_dt*nint(time/dump_space_dt))
!        else
!            tmp=1d99
!        endif
        if (dump_space.and.(time_savevol.ge.dump_space_dt.or.time.le.start_time.or.time.ge.end_time)) then
            call save_voldata
            !call save_srfdata
            if (solve_heat_eq) then
            call save_heatdata
            endif
            time_savevol=0d0
        endif

        ! dump restart file
        if (dump_restart_dt.gt.0d0) then
            tmp=2.d0*abs(time-dump_restart_dt*nint(time/dump_restart_dt))
        else
            tmp=1d99
        endif
        if (dump_restart.and.((tmp .lt. dt_face.and.time.gt.start_time).or.(time.ge.end_time))) then
            call store_restart(trim(restart_filename))
        endif

    end subroutine save






! store and restore state files
subroutine store_state(filename)
    character(*):: filename
    integer  unit_store,j, k,ios

    call set_unit(unit_store)
    open (unit=unit_store, file=filename,status='replace', form='formatted', iostat=ios)

    if (ios.eq.0) then
        do k=1,nspc
            do j=0,ngrd
                write (unit_store,*) dens(ndt,j,k)
            enddo
        enddo
        do k=1,nspc
            write (unit_store,*) dsrfl(ndt,k)
            write (unit_store,*) dsrfr(ndt,k)
        enddo

        do j=0,ngrd

            write (unit_store,*) temp(ndt,j)

        enddo
        close (unit_store)

        CALL PRINT_MILESTONE('Simulation state stored in '//TRIM(filename))
    else
        call face_error('cannot open state file: ', filename)
    endif

end subroutine store_state

subroutine restore_state(filename)
    character(*):: filename
    integer unit_restore,j, k, ios,i
     call set_unit(unit_restore)
    open(unit_restore, file=trim(filename), iostat=ios,action='read',status='old',form='formatted')

    if ( ios /= 0 ) then
        call face_error('Cannot open state file :', trim(adjustl(filename)))
    endif

    call print_milestone('Restoring state from file: '// trim(filename))


    do k=1,nspc
        do j=0,ngrd
            read (unit_restore,*) dens(ndt,j,k)
        enddo
    enddo

    do i=1,ndt
    dens(i,0:ngrd,1:nspc)=dens(ndt,0:ngrd,1:nspc)
    enddo


    do k=1,nspc
        read (unit_restore,*) dsrfl(ndt,k)
        read (unit_restore,*) dsrfr(ndt,k)
    enddo

    do i=1,ndt
    dsrfl(i,1:nspc)=dsrfl(ndt,1:nspc)
    dsrfr(i,1:nspc)=dsrfr(ndt,1:nspc)
    enddo

    if (restore_state_temp) then
    do j=0,ngrd
        read (unit_restore,*) temp(ndt,j)
    enddo
    do i=1,ndt
    temp(i,0:ngrd)=temp(ndt,0:ngrd)
    enddo
    endif

    close (unit_restore)

end subroutine restore_state


! ***** store and restore restart files *****

subroutine store_restart(filename)
    character(*):: filename
    integer i, j, k,ios,unit_restart
    call set_unit(unit_restart)
    open (unit=unit_restart, file=filename,status='replace', form='unformatted', iostat=ios)
    if (ios.eq.0) then
        do k=1,nspc
            do j=0,ngrd
                do i=1,ndt
                    write (unit_restart) dens(i,j,k), flx (i,j,k), ero_flx(i,j,k),cdif(i,j,k), rct (i,j,k), rate_d (i,j,k)
                enddo
            enddo
        enddo
        do j=0,ngrd
            write (unit_restart) x(j)
        enddo
        do k=1,nspc
            do i=1,ndt
                write (unit_restart) dsrfl(i,k), Gsrf_l(i,k)
                write (unit_restart) dsrfr(i,k), Gsrf_r(i,k)
                write (unit_restart) Gabs_l(i,k), Gdes_l(i,k), Gb_l(i,k), Gads_l(i,k)
                write (unit_restart) Gabs_r(i,k), Gdes_r(i,k), Gb_r(i,k), Gads_r(i,k)
                write (unit_restart) jout(i,k)
            enddo
        enddo
        do j=0,ngrd
            do i=1,ndt
                write (unit_restart) temp(i,j), qflx(i,j), rate_t(i,j), ero_qflx(i,j)
            enddo
        enddo
        write (unit_restart) time
        write (unit_restart) sfln_voldata,sfln_srfdata,sfln_heatdata
        close (unit_restart)
    else
        call face_error ('Cannot read restart file ', filename)
    endif
end subroutine store_restart

subroutine restore_restart(filename)
    character(*):: filename
    character*256 :: str
    integer i, j, k, ios,unit_restart
    call set_unit(unit_restart)
    open (unit=unit_restart, file=filename, form='unformatted',iostat=ios,&
        action='read',status='old')
    if ( ios .ne. 0 ) then
        call face_error('cannot open restart file : ',filename)
    endif
    call print_milestone('Restoring state from  restartfile: '// trim(filename))
    do k=1,nspc
        do j=0,ngrd
            do i=1,ndt
                read (unit_restart) dens(i,j,k), flx (i,j,k), ero_flx(i,j,k),cdif(i,j,k), rct(i,j,k), rate_d(i,j,k)
            enddo
        enddo
    enddo

    do j=0,ngrd
        read (unit_restart) x(j)
    enddo

    do k=1,nspc
        do i=1,ndt
            read (unit_restart) dsrfl(i,k), Gsrf_l(i,k)
            read (unit_restart) dsrfr(i,k), Gsrf_r(i,k)
            read (unit_restart) Gabs_l(i,k), Gdes_l(i,k), Gb_l(i,k), Gads_l(i,k)
            read (unit_restart) Gabs_r(i,k), Gdes_r(i,k), Gb_r(i,k), Gads_r(i,k)
            read (unit_restart) jout(i,k)
        enddo
    enddo

    do j=0,ngrd
        do i=1,ndt
            read (unit_restart) temp(i,j), qflx(i,j), rate_t(i,j), ero_qflx(i,j)
        enddo
    enddo

    read (unit_restart) time
    read (unit_restart) sfln_voldata,sfln_srfdata,sfln_heatdata

    close (unit_restart)

    call compute_inflx()

    write (str,'(a,1pe13.4e2,a)') '  *** simulation restarted with dbl precision file from t=',time, ' s'
    call print_milestone(str)
end subroutine restore_restart

! ***** restore routine: restore from restart file or state file *****
subroutine restore
    character*256 :: restart_file
    character*256 :: state_file
    if (verbose_restore) write(iout,*) 'read_restart_file:',read_restart_file
    if (verbose_restore) write(iout,*) 'read_state_file:',read_state_file
    ! cannot restore from restart and state file  at the same time
    if (trim(read_restart_file).ne."no".and.trim(read_state_file).ne."no") then
        call face_error('Cannot restore from restart file and state file simultaneously')
    endif

    ! restore from restart file?
    call print_milestone('Restoration from restart file : ' //trim(read_restart_file))
    if (trim(read_restart_file).eq."no") then

    elseif (read_restart_file.eq."yes") then
        restart_file=trim(path_folder)//'dsave.rst'
        call restore_restart(restart_file)
    else
        restart_file=read_restart_file
        call restore_restart(restart_file)
    endif

    ! restore from state file?
     call print_milestone('Restoration from restart file : ' //trim(read_state_file))
    if (trim(read_state_file).eq."no") then
    elseif (read_state_file.eq."yes") then
        state_file=trim(path_folder)//'face.state'
        call restore_state(state_file)
    else
        state_file=read_state_file
        call restore_state(state_file)
    endif

end subroutine restore

! ***** *****
subroutine store_final_state

    call store_state(final_state_file)
    call print_milestone('dumping file state completed')

end subroutine store_final_state

subroutine close_log()
   ! close log if connected to an open unit
    if (iout.ge.10) then
        close(iout)
    endif
end subroutine close_log





subroutine write_header_log
character*256::timestamp
if (verbose_header) then
call timestring ( timestamp )
write(iout,*) '# created      :  ', timestamp
write(iout,*) '# casename     :  ', trim(casename)
write(iout,*) '# path folder  :  ', trim(path_folder)
write(iout,*) '# data folder  :  ', trim(dat_folder)
endif
end subroutine write_header_log





subroutine write_grid

        integer j, ios,unit_grid
        character*256:: filename,myfmt2,myfmt3


        write(myfmt2,*) "(a4,a12)"
        write(myfmt3,*)"(i4.4, 12es3)"
        call set_unit(unit_grid)
        write (filename, '(a,a,a,a)') trim(path_folder),'/', trim(casename), '.grid'
        open (unit_grid, file=trim(filename), status='replace',iostat=ios)
        if (ios.eq.0) then
            write (unit_grid, '(a4,a12)') 'j','x'

            do j=0,ngrd
                write (unit_grid, '(i4.4, es12.3)')j,x(j)
            enddo

            close (unit_grid)
        else
            write (iout, '(2a)') ' *** error saving ', trim(filename)
        endif



    end subroutine write_grid

    subroutine print_formatted(str)
    character(*)::str
    write(iout,'(a5,a)') ' ',trim(str)
    end subroutine print_formatted

    subroutine print_section(str)
    character(*)::str
    write(iout,'(a5,a,a,a)') ' ','*************** ',trim(str),' ***************'
    end subroutine print_section

     subroutine print_line(str)
    character(*)::str
    write(iout,'(a5,a,a)') ' ','* ',trim(str)
    end subroutine print_line

    subroutine print_end_section(str)
character(*)::str
character*256::str2
integer :: i
     i=len_trim(str)
     str2 = REPEAT('*',i)
    write(iout,'(a5,a,a,a)') ' ','****************',trim(str2),'****************'
     write(iout,'(a)') ' '
    end subroutine print_end_section

    subroutine print_headline(str)
character(*)::str
character*256::str2,str3
integer :: i
     i=40
     str2 = REPEAT('>',i)
     str3 = REPEAT('<',i)
     write(iout,'(a)') ' '
    write(iout,'(a,a,a,a,a)') trim(str2),' ',trim(str),' ',trim(str3)
     write(iout,'(a)') ' '
    end subroutine print_headline

    subroutine print_timestep_info(convergence)
    logical,intent(in)::convergence
character*256::myfmt,str
real::tot=0d0,fin=0d0,fn=0d0
integer:: k

if (mod(iteration,Nprint_run_info).eq.0.or.(time.ge.end_time).or.(time.le.start_time)) then
if (convergence) then
write(myfmt,*) "('++ iter=', 1I6,' time=', es14.6, 's;   T_l=',es9.2, 'K ; T_r=',es9.2, 'K; "&
 ," |f|=',es9.2, ' iter_solver=',i3,' dt=',es9.2)"
else
write(myfmt,*) "('-- iter=', 1I6,' time=', es14.6, 's;   T_l=',es9.2, 'K ; T_r=',es9.2, 'K; "&
 ," |f|=',es9.2, ' iter_solver=',i3,' dt=',es9.2)"
endif
 write (str, myfmt) iteration,time, temp(ndt,0), temp(ndt,ngrd),normf,iter_solver,dt_face
  call print_formatted(str)
if (print_onthefly_inventory) then
write(myfmt,*) "(' - onthefly inventory: k=', 1I2,' src+des=', es10.3,' net_int_dens=',",&
"es10.3,' net_int_dsrf=', es10.3,' tot=', es10.3,' f-n=',es9.3,' f-src=',es9.3)"

do k=1,nspc
tot=-onthefly_int_des(k)+onthefly_int_src(k)&
-(onthefly_net_int_dens(k)+onthefly_net_int_dsrf(k))
fn=abs(tot)/abs(onthefly_net_int_dens(k)+onthefly_net_int_dsrf(k))
fin=abs(tot)/abs(onthefly_int_src(k))
write (str, myfmt) k,-onthefly_int_des(k)+onthefly_int_src(k),onthefly_net_int_dens(k),&
onthefly_net_int_dsrf(k),tot,fn,fin
  call print_formatted(str)
enddo
endif
endif

end subroutine print_timestep_info

subroutine print_reduction_timestep_info(count_loop)
integer,intent(in):: count_loop
character*256::myfmt,str


!write(myfmt,*) "('dt REDUCTION: iter=', 1I6,' time=', es9.2, 's;   T_l=',es9.2, 'K;   T_r=',es9.2, 'K;"&
! ,"dt=', es9.2, 's, |f|=',es9.2, ' iter_solver=',i3)"
 write(myfmt,*) "('dt REDUCTION: iter=', 1I6,' time=', es9.2, ' dt=', es9.2, 's iloop=',i3)"
 !write (str, myfmt) iteration,time, temp(ndt,0), temp(ndt,ngrd), dt_face,normf,iter_solver
  write (str, myfmt) iteration,time, dt_face,count_loop
  call print_formatted(str)

end subroutine print_reduction_timestep_info

subroutine print_steady_timestep_info
if (mod(iteration,Nprint_run_info).eq.0.or.(time.ge.end_time).or.(time.le.start_time)) then
  call print_formatted('>>>> dt steady due to non-convergence with larger timestep')
  endif

end subroutine print_steady_timestep_info

subroutine print_milestone(str)
character(*)::str
write(iout,"(a5,'--- ',a,' ---')")" ", str
write(iout,"(a)") " "
end subroutine print_milestone

subroutine print_vector(u,str)
real,intent(in) :: u(neq)
integer :: i
character(*) :: str
do i=1,neq
        write(iout,*) str,'(',i,')=',u(i)
            enddo
end subroutine print_vector
subroutine set_unit(ifile)
      integer ifile
      logical unit_open
      unit_open = .true.
      current_ifile=10
      Do While (unit_open.and.current_ifile.le.max_ifile)
         current_ifile=current_ifile + 1
         Inquire (Unit = current_ifile, Opened = Unit_open)
      End Do
      if (current_ifile.lt.max_ifile) then
      ifile=current_ifile
      else
      write(iout,*) "ERROR: cannot find available unit number <", max_ifile
      stop
      endif

      return


   end subroutine set_unit
end module IO
