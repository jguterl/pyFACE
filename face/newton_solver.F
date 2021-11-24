      module newton
      use SolverHeader
      use IOHeader
      use VerboseHeader
      use compute
      contains
      subroutine newton_solver(quick_convergence)
      implicit none
      integer i,iter_solver_max_
      integer  cntl, idx
      logical:: quick_convergence
      real:: u(neq), du(neq), f(neq), fdot(neq,neq)
      real:: unew(neq), unorm(neq)
      real:: norm, normnew
c      solver_status: 0: solver step not completed
c                   : -1: solver step completed with iter_solver=max
c                   : 1 : solver step completed iwht iter_solver<ideal
       quick_convergence=.false.
       iter_solver=0
       if (verbose_debug) write(iout,*) 'newton_solver'
       if (iteration>10) then
       iter_solver_max_=iter_solver_max
       else
       iter_solver_max_=iter_solver_max_first
       endif
      call build_vector(u,du)



      norm=compute_fnorm(u)

      call compute_f(u,f)

c      if (verbose_debug) then
c      do ii=1,neq
c      write(iout,*) 'i=',ii,'u=',u(ii),'f(i)=',f(ii)
c      enddo
c      endif

5199  call jac(u,f,fdot,norm)
      call dsolve(du,u,f,fdot)
c
c     --- newton step reduction ---
      norm=compute_fnorm(u)

      cntl=0
5198  idx=0

      do i=1,neq

       unew (i)=u(i)+du(i)

       if (u(i).eq.0.0) then
       unorm(i)=0.0
       else
        unorm(i)=du(i)/u(i)
c     if (verbose_step) write(iout,*) "i=",i," ;unorm=",unorm(i)
c     endif
       if (abs(unorm(i)) .lt. 1.d-30) then
c     write (*,*) '***warning: displacement vector is too small at '
c     +              , i
        unorm(i)=1d-30
       endif
       endif
       if ((-unorm(i) .gt. solver_udspl) .or.
     +     ( unorm(i) .gt. 1.d0/(1.d0-solver_udspl))) then
        idx=i
      if (verbose_step) write(iout,*) "unorm too large: idx=",idx
       endif
      enddo


      if (idx .ne. 0) then
       do i=1,neq
        du(i)=solver_fstp*du(i)
       enddo
       cntl=cntl+1
       if (verbose_step) write (*,*) cntl, ' var step reduction',
     +idx
       goto 5198
      endif


      normnew=compute_fnorm(unew)


      if ((normnew/norm-1.d0) .gt. solver_fdspl) then

       do i=1,neq
        du(i)=solver_fstp*du(i)
       enddo
       cntl=cntl+1
        if (verbose_step) then
       write (iout,*) cntl, ' norm step reduction', normnew
       endif
       goto 5198
      else
       do i=1,neq
        u(i)=unew(i)
       enddo
       norm=normnew
      endif

      call compute_f(u,f)

c     --- check convergence ---
      iter_solver=iter_solver+1
      if (iter_solver .lt. iter_solver_max_) then
      if ((verbose_step).OR.(verbose_debug)) then
      write (iout,*)'-- Newton iter# ',iter_solver,
     +' norm ', norm
      endif

       if (norm .gt. solver_eps) then
       if (verbose_debug) then
      write(iout,*) 'jump to 5199'
      endif
       goto 5199
      endif
       ! if norm<eps then do nothing and exit if at line 107
      else
       if (verbose_step) then
       write (iout,*) 'iter# ', iter_solver,
     * ' norm ', norm
       call face_warning('implicit precision was not reached!')
      endif

      endif

      ! update values of densities and temperature only if we are happy with the convergence...
      if (iter_solver.lt.iter_solver_max_) then
      quick_convergence=.true.
      endif

      if (quick_convergence) then
      call get_density_values(u)
      endif

      normf=norm

      end subroutine newton_solver
      end module