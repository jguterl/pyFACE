module solver
    use header
    implicit none
contains
    subroutine jac(u,f,fdot,norm)
        integer i, j
        real,intent(in):: u(:), f(:)
        real::  um(neq), up(neq), fm(neq), fp(neq)
        real:: fdot(neq,neq)
        real::  ddu, uud, amax, norm
        real:: ftran
        parameter (ftran=1.d-0)

        do i=1,neq
            um(i)=u(i)
            up(i)=u(i)
        enddo

        do i=1,neq
            ddu=u(i)*jac_eps
            um(i)=u(i)-ddu
            up(i)=u(i)+ddu

            call compute_f(um,fm)
            call compute_f(up,fp)

            uud=1.d0/(ddu+ddu)
            do j=1,neq
                fdot(j,i)=(fp(j)-fm(j))*uud
            enddo
            um(i)=u(i)
            up(i)=u(i)
        enddo
        !     --- pseudo-transient continuation ---
        if (steady_state) then
            do i=1,neq
                fdot(i,i)=fdot(i,i)-ftran*norm
            enddo
        endif
        !     --- check for singular Jacobian ---
        do i=1,neq
            amax=0.d0
            do j=1,neq
                if (abs(fdot(i,j)) .gt. amax) amax=abs(fdot(i,j))
            enddo
            if (amax .eq. 0.d0) then
                write (iout,*) '***warning: underflow in Jacobian raw ', i
                write (iout,*) '***fdot(i,i)= ', fdot(i,i),fp(i)
                call which_eq(i)
                fdot(i,i)=-1.d+99
            endif
        enddo

    end subroutine jac

    subroutine build_vector(u,du)

        real:: u(:), du(:)
        integer i,k,j
        i=0
if(compute_spc) then
        do k=1,nspc
            i=i+1
            u (i)=dsrfl(ndt-1,k)
            if (isnan(u (i))) call face_error(" dsrfl NaN")
            du(i)=0.d0
            do j=0,ngrd
                i=i+1
                u (i)=dens(ndt-1,j,k)
                if (isnan(u (i))) call face_error(" dens NaN")
                du(i)=0.d0
            enddo
            i=i+1
            u (i)=dsrfr(ndt-1,k)
            if (isnan(u (i))) call face_error(" dsrfr NaN")
            du(i)=0.d0
        enddo
        endif
        if (solve_heat_eq) then
            do j=0,ngrd
                i=i+1
                u (i)=temp(ndt-1,j)
                du(i)=0.d0
            enddo
        endif
        if (i.ne.neq) then
            call face_error("mismatch in vector size: neq=",neq, 'size(vec)=',i)

        endif
    call check_isNaN(u,'build_vector')
    end subroutine build_vector

    subroutine dsolve(du,u,f,fdot)

        real,intent(in)::  fdot(:,:), u(:), f(:)
        real,intent(out) :: du(:)
        real:: a(neq,neq), b(neq)
        integer indx(neq)
        integer n, d, code
        integer i, j

        n=neq
        do i=1,neq
            b(i)=-f(i)
            do j=1,neq
                a(i,j)=fdot(i,j)
            enddo
        enddo
        ! solve A.X=B where A (N,N), B(N) X(N). Result is returned as b
        call solve_linear_system(a,b,n)

        do i=1,neq
            if (u(i) .gt. 1.d-10) then
                du(i)=b(i)
            else
                if (b(i) .gt. 0.d0) then
                    du(i)=b(i)
                else
                    du(i)=0.d0
                endif
            endif
        enddo

    !      do i=1,neq
    !!!       c(i)=0.d0
    !!!       do j=1,neq
    !!!        c(i)=c(i)+fdot(i,j)*du(j)
    !!!       enddo
    !!!       c(i)=c(i)+f(i)
    !!!       if (abs(c(i)) .gt. 1.d-6) then
    !!!        write (*,*) '***warning: loss of direction precision ',
    !!!     +              i, c(i)
    !!!       endif
    !!!      enddo

    end subroutine dsolve


end module solver
