module maths
    use header
    implicit none
    contains
    subroutine solve_linear_system(a,b,n)
    integer,intent(in):: n
    real:: a(n,n)
    real,intent(inout):: b(n)
    integer :: indx(n), code,d
    !     Externals from the LAPACK library
      external dgesv
!       DGESV computes the solution to a real system of linear equations
!    A * X = B,
! where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
!
! The LU decomposition with partial pivoting and row interchanges is
! used to factor A as
!    A = P * L * U,
! where P is a permutation matrix, L is unit lower triangular, and U is
! upper triangular.  The factored form of A is then used to solve the
! system of equations A * X = B.
!Parameters
![in]    N
!          N is INTEGER
!          The number of linear equations, i.e., the order of the
!          matrix A.  N >= 0.
![in]    NRHS
!          NRHS is INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
![in,out]    A
!          A is DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the N-by-N coefficient matrix A.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
![in]    LDA
!          LDA is INTEGER
 !         The leading dimension of the array A.  LDA >= max(1,N).
![out]   IPIV
!          IPIV is INTEGER array, dimension (N)
 !         The pivot indices that define the permutation matrix P;
 !         row i of the matrix was interchanged with row IPIV(i).
![in,out]    B
 !         B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the N-by-NRHS matrix of right hand side matrix B.
!          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
![in]    LDB
!          LDB is INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
![out]   INFO
!          INFO is INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
!                has been completed, but the factor U is exactly
!                singular, so the solution could not be computed.
      !

#ifdef MLKSOLVER
!            if (DP.ne.8) call face_error("dgsev solver for double precision only")
            call dgesv(n, 1, a, n, indx, b, n, code)
                    if (code .ne. 0) then
            call face_error("Singular Jacobian! code=",code)
         endif
#else
    call DLUDCMP(a,n,indx,d,code)
    if (code .ne. 1) then
            call DLUBKSB(a,n,indx,b)
        else
            call face_error("Singular Jacobian!")

        endif
#endif



    end subroutine solve_linear_system


subroutine which_eq(j_eq)

        integer,intent(in):: j_eq
        integer j_idx,k_idx

        if (solve_heat_eq) then

            if (j_eq.ge.neq-ngrd+1.and.j_eq.le.neq) then

                j_idx=ngrd+1-(neq-j_eq)
                write(iout,*) "eq: heat eq at j=",j_idx

            elseif (j_eq.lt.neq-ngrd+1) then


                k_idx=floor(real(j_eq)/real(ngrd+3))
                j_idx=j_eq-k_idx*(ngrd+3)+1
                if  (j_idx.eq.1) then
                    write(iout,*) "eq: spc k = ",k_idx," at left surface "
                elseif  (j_idx.eq.ngrd+3) then
                    write(iout,*) "eq: spc k = ",k_idx," at right surface"
                elseif (j_idx.lt.ngrd+3.and.j_idx.gt.1) then
                    write(iout,*) "eq: volume spc k = ",k_idx," at j = ",j_idx-1
                else
                    call face_error("cannot find the equation ",j_eq)
                endif
            endif

        else

                k_idx=floor(real(j_eq)/real(ngrd+3))
                j_idx=j_eq-k_idx*(ngrd+3)+1
                if  (j_idx.eq.1) then
                    write(iout,*) "eq: spc k = ",k_idx," at left surface "
                elseif  (j_idx.eq.ngrd+3) then
                    write(iout,*) "eq: spc k = ",k_idx," at right surface"
                elseif (j_idx.lt.ngrd+3.and.j_idx.gt.1) then
                    write(iout,*) "eq: volume spc k = ",k_idx," at j = ",j_idx-1
                else
                    call face_error("cannot find the equation ;", j_eq,"/",neq,"; jdix=",j_idx)
                endif
            !else
             !   call face_error("j_eq too large j_eq=",j_eq,"; neq= ",neq)
            endif
    end subroutine which_eq

            ! ***************************************************************
    ! * Given an N x N matrix A, this routine replaces it by the LU *
    ! * decomposition of a rowwise permutation of itself. A and N *
    ! * are input. INDX is an output vector which records the row *
    ! * permutation effected by the partial pivoting; D is output *
    ! * as -1 or 1, depending on whether the number of row inter- *
    ! * changes was even or odd, respectively. This routine is used *
    ! * in combination with LUBKSB to solve linear equations or to *
    ! * invert a matrix. Return code is 1, if matrix is singular. *
    ! ***************************************************************
     Subroutine DLUDCMP(A,N,INDX,D,CODE)
     IMPLICIT NONE
     integer, parameter :: nmax = 10000
     real, parameter :: tiny = 1.d-30

     real, intent(inout), dimension(N,N) :: A
     integer, intent(in) :: N
     integer, intent(out) :: D, CODE
     integer, intent(out), dimension(N) :: INDX
     !f2py depend(N) A, indx

     real :: AMAX, DUM, SUMM, VV(NMAX)
     INTEGER :: i, j, k, imax

     D=1; CODE=0

     DO I=1,N
       AMAX=0.d0
       DO J=1,N
         IF (ABS(A(I,J)).GT.AMAX) AMAX=ABS(A(I,J))
       END DO ! j loop
       IF(AMAX.EQ.0.d0) THEN
         CODE = 1
         write(iout,*) i, amax
         RETURN

       END IF
       VV(I) = 1.d0 / AMAX
     END DO ! i loop

     DO J=1,N
       DO I=1,J-1
         SUMM = A(I,J)
         DO K=1,I-1
           SUMM = SUMM - A(I,K)*A(K,J)
         END DO ! k loop
         A(I,J) = SUMM
       END DO ! i loop
       AMAX = 0.d0
       DO I=J,N
         SUMM = A(I,J)
         DO K=1,J-1
           SUMM = SUMM - A(I,K)*A(K,J)
         END DO ! k loop
         A(I,J) = SUMM
         DUM = VV(I)*ABS(SUMM)
         IF(DUM.GE.AMAX) THEN
           IMAX = I
           AMAX = DUM
         END IF
       END DO ! i loop

       IF(J.NE.IMAX) THEN
         DO K=1,N
           DUM = A(IMAX,K)
           A(IMAX,K) = A(J,K)
           A(J,K) = DUM
         END DO ! k loop
         D = -D
         VV(IMAX) = VV(J)
       END IF

       INDX(J) = IMAX
       IF(A(J,J).EQ.0.d0) THEN
         if (verbose_maths) write(iout,*) 'Warning: Newton direction is not precise', j
         if (verbose_maths) call which_eq(j)
         A(J,J) = TINY
       END IF

       IF(J.NE.N) THEN
         DUM = 1.d0 / A(J,J)
         DO I=J+1,N
           A(I,J) = A(I,J)*DUM
         END DO ! i loop
       END IF
     END DO ! j loop

     RETURN
END subroutine DLUDCMP


! ******************************************************************
! * Solves the set of N linear equations A . X = B. Here A is *
! * input, not as the matrix A but rather as its LU decomposition, *
! * determined by the routine LUDCMP. INDX is input as the permuta-*
! * tion vector returned by LUDCMP. B is input as the right-hand *
! * side vector B, and returns with the solution vector X. A, N and*
! * INDX are not modified by this routine and can be used for suc- *
! * cessive calls with different right-hand sides. This routine is *
! * also efficient for plain matrix inversion. *
! ******************************************************************
 Subroutine DLUBKSB(A, N, INDX, B)
 integer, intent(in) :: N
 real, intent(in), dimension(N,N) :: A
 integer, intent(in), dimension(N) :: INDX
 real, intent(inout), dimension(N) :: B
 integer :: I,J,II,LL
 !f2py depend(N) A, INDX, B

 real SUMM

 II = 0

 DO I=1,N
   LL = INDX(I)
   SUMM = B(LL)
   B(LL) = B(I)
   IF(II.NE.0) THEN
     DO J=II,I-1
       SUMM = SUMM - A(I,J)*B(J)
     END DO ! j loop
   ELSE IF(SUMM.NE.0.d0) THEN
     II = I
   END IF
   B(I) = SUMM
 END DO ! i loop

 DO I=N,1,-1
   SUMM = B(I)
   IF(I < N) THEN
     DO J=I+1,N
       SUMM = SUMM - A(I,J)*B(J)
     END DO ! j loop
   END IF
   B(I) = SUMM / A(I,I)
 END DO ! i loop

 RETURN
END subroutine DLUBKSB
   end module maths
