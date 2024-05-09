MODULE MATRIX
  USE NRUTIL, ONLY : assert_eq,swap,outerand,nrerror,imaxloc
  USE toolbox

  
  IMPLICIT NONE

  INTERFACE Outer_Product
     MODULE PROCEDURE Outer_Product_r,Outer_Product_i
  END INTERFACE
  INTERFACE Add_Diagonal
     MODULE PROCEDURE Add_Diagonal_rv,Add_Diagonal_r,Add_Diagonal_iv,Add_Diagonal_i
  END INTERFACE
  INTERFACE Multiply_Diagonal
     MODULE PROCEDURE Multiply_Diagonal_rv,Multiply_Diagonal_r,Multiply_Diagonal_iv,Multiply_Diagonal_i
  END INTERFACE
  INTERFACE Put_Diagonal
     MODULE PROCEDURE Put_Diagonal_rv, Put_Diagonal_r,Put_Diagonal_iv, Put_Diagonal_i
  END INTERFACE
  INTERFACE Get_Diagonal
     MODULE PROCEDURE Get_Diagonal_r, Get_Diagonal_i
  END INTERFACE

  !PRIVATE :: Jacobi,ludcmp,lubksb
  
CONTAINS
  
  FUNCTION Outer_Product_r(a,b)
    REAL(8), DIMENSION(:), INTENT(IN) :: a,b
    REAL(8), DIMENSION(size(a),size(b)) :: Outer_Product_r
    Outer_Product_r = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION Outer_Product_r
  
  FUNCTION Outer_Product_i(a,b)
    INTEGER, DIMENSION(:), INTENT(IN) :: a,b
    INTEGER, DIMENSION(size(a),size(b)) :: Outer_Product_i
    Outer_Product_i = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION Outer_Product_i

  SUBROUTINE Add_Diagonal_rv(mat,diag)
    REAL(8), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(8), DIMENSION(:), INTENT(IN) :: diag
    INTEGER :: j,n
    n = assert_eq(size(diag),min(size(mat,1),size(mat,2)),'Add_Diagonal_rv')
    do j=1,n
       mat(j,j)=mat(j,j)+diag(j)
    end do
  END SUBROUTINE Add_Diagonal_rv

  SUBROUTINE Add_Diagonal_r(mat,diag)
    REAL(8), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(8), INTENT(IN) :: diag
    INTEGER :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=mat(j,j)+diag
    end do
  END SUBROUTINE Add_Diagonal_r

  SUBROUTINE Add_Diagonal_iv(mat,diag)
    INTEGER, DIMENSION(:,:), INTENT(INOUT) :: mat
    INTEGER, DIMENSION(:), INTENT(IN) :: diag
    INTEGER :: j,n
    n = assert_eq(size(diag),min(size(mat,1),size(mat,2)),'Add_Diagonal_iv')
    do j=1,n
       mat(j,j)=mat(j,j)+diag(j)
    end do
  END SUBROUTINE Add_Diagonal_iv

  SUBROUTINE Add_Diagonal_i(mat,diag)
    INTEGER, DIMENSION(:,:), INTENT(INOUT) :: mat
    INTEGER, INTENT(IN) :: diag
    INTEGER :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=mat(j,j)+diag
    end do
  END SUBROUTINE Add_Diagonal_i

  SUBROUTINE Multiply_Diagonal_rv(mat,diag)
    REAL(8), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(8), DIMENSION(:), INTENT(IN) :: diag
    INTEGER :: j,n
    n = assert_eq(size(diag),min(size(mat,1),size(mat,2)),'diagmult_rv')
    do j=1,n
       mat(j,j)=mat(j,j)*diag(j)
    end do
  END SUBROUTINE Multiply_Diagonal_rv
  
  SUBROUTINE Multiply_Diagonal_r(mat,diag)
    REAL(8), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(8), INTENT(IN) :: diag
    INTEGER :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=mat(j,j)*diag
    end do
  END SUBROUTINE Multiply_Diagonal_r

  SUBROUTINE Multiply_Diagonal_iv(mat,diag)
    INTEGER, DIMENSION(:,:), INTENT(INOUT) :: mat
    INTEGER, DIMENSION(:), INTENT(IN) :: diag
    INTEGER :: j,n
    n = assert_eq(size(diag),min(size(mat,1),size(mat,2)),'Multiply_Diagonal_rv')
    do j=1,n
       mat(j,j)=mat(j,j)*diag(j)
    end do
  END SUBROUTINE Multiply_Diagonal_iv
  
  SUBROUTINE Multiply_Diagonal_i(mat,diag)
    INTEGER, DIMENSION(:,:), INTENT(INOUT) :: mat
    INTEGER, INTENT(IN) :: diag
    INTEGER :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=mat(j,j)*diag
    end do
  END SUBROUTINE Multiply_Diagonal_i

  FUNCTION Get_Diagonal_r(mat)
    REAL(8), DIMENSION(:,:), INTENT(IN) :: mat
    REAL(8), DIMENSION(size(mat,1)) :: Get_Diagonal_r
    INTEGER :: j
    j=assert_eq(size(mat,1),size(mat,2),'Get_Diagonal_r')
    do j=1,size(mat,1)
       Get_Diagonal_r(j)=mat(j,j)
    end do
  END FUNCTION Get_Diagonal_r

  FUNCTION Get_Diagonal_i(mat)
    INTEGER, DIMENSION(:,:), INTENT(IN) :: mat
    INTEGER, DIMENSION(size(mat,1)) :: Get_Diagonal_i
    INTEGER :: j
    j=assert_eq(size(mat,1),size(mat,2),'Get_Diagonal_i')
    do j=1,size(mat,1)
       Get_Diagonal_i(j)=mat(j,j)
    end do
  END FUNCTION Get_Diagonal_i

  SUBROUTINE Put_Diagonal_rv(diagv,mat)
    REAL(8), DIMENSION(:), INTENT(IN) :: diagv
    REAL(8), DIMENSION(:,:), INTENT(INOUT) :: mat
    INTEGER :: j,n
    n=assert_eq(size(diagv),min(size(mat,1),size(mat,2)),'Put_Diagonal_rv')
    do j=1,n
       mat(j,j)=diagv(j)
    end do
  END SUBROUTINE Put_Diagonal_rv

  SUBROUTINE Put_Diagonal_r(scal,mat)
    REAL(8), INTENT(IN) :: scal
    REAL(8), DIMENSION(:,:), INTENT(INOUT) :: mat
    INTEGER :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=scal
    end do
  END SUBROUTINE Put_Diagonal_r

  SUBROUTINE Put_Diagonal_iv(diagv,mat)
    INTEGER, DIMENSION(:), INTENT(IN) :: diagv
    INTEGER, DIMENSION(:,:), INTENT(INOUT) :: mat
    INTEGER :: j,n
    n=assert_eq(size(diagv),min(size(mat,1),size(mat,2)),'Put_Diagonal_iv')
    do j=1,n
       mat(j,j)=diagv(j)
    end do
  END SUBROUTINE Put_Diagonal_iv

  SUBROUTINE Put_Diagonal_i(scal,mat)
    INTEGER, INTENT(IN) :: scal
    INTEGER, DIMENSION(:,:), INTENT(INOUT) :: mat
    INTEGER :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=scal
    end do
  END SUBROUTINE Put_Diagonal_i

  SUBROUTINE Identity(mat)
    REAL(8), DIMENSION(:,:), INTENT(OUT) :: mat
    INTEGER :: i,n
    n=min(size(mat,1),size(mat,2))
    mat(:,:)=0.0d0
    do i=1,n
       mat(i,i)=1.0d0
    end do
  END SUBROUTINE Identity

  FUNCTION Matrix_Inverse(ain)
    ! Inverts a square matrix by LU decomposition
    REAL(8), INTENT(IN) :: ain(:,:)
    REAL(8) :: Matrix_Inverse(SIZE(ain,1),SIZE(ain,2))
    INTEGER :: j,indx(SIZE(ain,1)),n,sing
    REAL(8) :: aux(SIZE(ain,1),SIZE(ain,2)),b(SIZE(ain,1)),lu(SIZE(ain,1),SIZE(ain,2)),d
    n=assert_eq(SIZE(ain,1),SIZE(ain,2),'Matrix_Inverse')
    aux=ain
    CALL ludcmp(aux,indx,d,sing)
    if (sing==1) call nrerror('singular matrix in Matrix Inverse')    
    DO j=1,n
       b=0.0d0
       b(j)=1.0d0
       CALL lubksb(aux,indx,b)
       Matrix_Inverse(:,j)=b
    END DO
  END FUNCTION Matrix_Inverse
  
  FUNCTION Matrix_Inverse_Jor(AIN)
    ! Inverts a Square Matrix
    IMPLICIT NONE
    REAL(8), DIMENSION(:,:), INTENT(IN) :: AIN
    INTEGER, DIMENSION(SIZE(AIN,1)) :: ipiv,indxr,indxc
    LOGICAL, DIMENSION(SIZE(AIN,1)) :: lpiv
    REAL(8) :: pivinv,A(SIZE(AIN,1),SIZE(AIN,2)),Matrix_Inverse_Jor(SIZE(AIN,1),SIZE(AIN,2))
    REAL(8), DIMENSION(SIZE(AIN,1)) :: dumc
    INTEGER, TARGET :: irc(2)
    INTEGER :: i,l,n
    INTEGER, POINTER :: irow,icol
    n=assert_eq(SIZE(ain,1),size(ain,2),'Matrix_Inverse_Jordan')
    A=AIN
    irow => irc(1)
    icol => irc(2)
    ipiv=0
    DO i=1,n
       lpiv = (ipiv == 0)
       irc=MAXLOC(ABS(A),outerand(lpiv,lpiv))
       ipiv(icol)=ipiv(icol)+1
       IF (ipiv(icol) > 1) CALL NRERROR('singular matrix(1): Matrix_Inverse_Jor')
       IF (irow /= icol) THEN
          CALL swap(A(irow,:),A(icol,:))
       END IF
       indxr(i)=irow
       indxc(i)=icol
       IF (a(icol,icol) == 0.0) CALL NRERROR('singular matrix (2): Matrix_Inverse_Jor')
       pivinv=1.0d0/A(icol,icol)
       A(icol,icol)=1.0d0
       A(icol,:)=A(icol,:)*pivinv
       dumc=A(:,icol)
       A(:,icol)=0.0
       A(icol,icol)=pivinv
       A(1:icol-1,:)=A(1:icol-1,:)-Outer_Product(dumc(1:icol-1),A(icol,:))
       A(icol+1:,:)=A(icol+1:,:)-Outer_Product(dumc(icol+1:),A(icol,:))
    END DO
    DO l=n,1,-1
       CALL swap(A(:,indxr(l)),A(:,indxc(l)))
    END DO
    Matrix_Inverse_Jor=A
  END FUNCTION Matrix_Inverse_Jor

  FUNCTION Cholesky(AIN)
    IMPLICIT NONE
    REAL(8), DIMENSION(:,:), INTENT(IN) :: AIN
    REAL(8) :: A(SIZE(AIN,1),SIZE(AIN,2)),p(SIZE(AIN,1)),Cholesky(SIZE(AIN,1),SIZE(AIN,2))
    INTEGER :: i,j,n
    REAL(8) :: summ
    n=assert_eq(size(ain,1),size(ain,2),'Cholesky')    
    A=AIN
    Cholesky=0.0d0
    DO i=1,n
       summ=A(i,i)-DOT_PRODUCT(A(i,1:i-1),A(i,1:i-1))
       IF (summ <= 0.0d0) CALL NRERROR('Cholesky failed')
       p(i)=SQRT(summ)
       A(i+1:n,i)=(A(i,i+1:n)-matmul(A(i+1:n,1:i-1),A(i,1:i-1)))/p(i)
       A(i,i)=p(i)
    END DO
    FORALL (i=1:n,j=1:n,j<=i) Cholesky(i,j)=A(i,j)
  END FUNCTION Cholesky

  FUNCTION Determinant(B)
    ! Subroutine for evaluating the determinant of a matrix using 
    ! the partial-pivoting Gaussian elimination scheme.
    IMPLICIT NONE
    REAL(8), INTENT (IN) :: B(:,:)
    INTEGER :: I,J,MSGN,N,INDX(SIZE(B,1))
    REAL(8) :: Determinant,A(SIZE(B,1),SIZE(B,1)),D
    A=B
    N=SIZE(B,1)
    CALL ELGS(A,N,INDX)
    D = 1.0
    DO I = 1, N
       D = D*A(INDX(I),I)
    END DO
    MSGN = 1
    DO I = 1, N
       DO WHILE (I.NE.INDX(I))
          MSGN = -MSGN
          J = INDX(I)
          INDX(I) = INDX(J)
          INDX(J) = J
       END DO
    END DO
    Determinant = MSGN*D
  CONTAINS
    SUBROUTINE ELGS (A,N,INDX)
      ! Subroutine to perform the partial-pivoting Gaussian elimination.
      ! A(N,N) is the original matrix in the input and transformed matrix
      ! plus the pivoting element ratios below the diagonal in the output.
      ! INDX(N) records the pivoting order.  Copyright (c) Tao Pang 2001.
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: N
      INTEGER :: I,J,K,ITMP
      INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
      REAL(8) :: C1,PI,PI1,PJ
      REAL(8), INTENT (INOUT), DIMENSION (N,N) :: A
      REAL(8), DIMENSION (N) :: C
      ! Initialize the index
      DO I = 1, N
         INDX(I) = I
      END DO
      ! Find the rescaling factors, one from each row
      DO I = 1, N
         C1= 0.0
         DO J = 1, N
            C1 = DMAX1(C1,ABS(A(I,J)))
         END DO
         C(I) = C1
      END DO
      ! Search the pivoting (largest) element from each column
      DO J = 1, N-1
         PI1 = 0.0
         DO I = J, N
            PI = ABS(A(INDX(I),J))/C(INDX(I))
            IF (PI.GT.PI1) THEN
               PI1 = PI
               K   = I
            ENDIF
         END DO
         ! Interchange the rows via INDX(N) to record pivoting order
         ITMP    = INDX(J)
         INDX(J) = INDX(K)
         INDX(K) = ITMP
         DO I = J+1, N
            PJ  = A(INDX(I),J)/A(INDX(J),J)
            ! Record pivoting ratios below the diagonal
            A(INDX(I),J) = PJ
            ! Modify other elements accordingly
            DO K = J+1, N
               A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
            END DO
         END DO
      END DO
    END SUBROUTINE ELGS
  END FUNCTION Determinant

  FUNCTION ols(x,y)
    REAL(8), intent(in) :: x(:,:), y(:)
    REAL(8) :: ols(size(x,2))
    REAL(8) :: xx(size(x,2),size(x,2)),d
    INTEGER :: index(size(x,2)),sing
    xx=matmul(transpose(x),x)
    ols=matmul(transpose(x),y)
    call ludcmp(xx,index,d,sing)
    if (sing==1) call nrerror('singular matrix in OLS')
    call lubksb(xx,index,ols)
  END FUNCTION ols
  
  INTEGER FUNCTION Rank(A)
    ! RETURNS THE RANK OF A MATRIX
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: A(:,:)
    REAL(8) :: AA(SIZE(A,2),SIZE(A,2)),EVAL(SIZE(A,2)),TOL,EVEC(SIZE(A,2),SIZE(A,2))
    INTEGER :: N,nrot
    TOL=1.0d-6
    N = SIZE(A,2)
    AA = MATMUL(TRANSPOSE(A),A)
    CALL JACOBI(AA,EVAL,EVEC,nrot)
    RANK = COUNT(EVAL>TOL)
  END FUNCTION Rank

  SUBROUTINE Colinear(X,IND,N)
    ! CHECKS FOR LINEAR INDEPENDENCE IN A MATRIX
    ! IT RETURNS THE INDEX OF N LINEARLY INDEPENDENT COLUMNS OF THE MATRIX IN IND(1:N)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: X(:,:)
    INTEGER, INTENT(OUT) :: IND(SIZE(X,2)),N
    INTEGER :: h,NR,NC,IN(SIZE(X,2)),R,R1,NC1,RS
    NR = SIZE(X,1)
    NC = SIZE(X,2)
    N = RANK(X)
    DO h = 1, NC
       IN(h) = h
    END DO
    IF (N==NC) THEN
       IND = IN
    ELSE
       RS = N
       NC1=0
       DO h = 1, NC
          R = NC - h
          R1 = RANK(X(:,IN(h+1:NC)))
          IF ((R1<RS).and.(R1<R)) THEN
             NC1 = NC1 + 1
             IND(NC1) = h
             RS=RS-1
          ELSEIF ((R1.ge.RS).and.(R1<R)) THEN
             RS=RS
          ELSE
             IND(NC1+1:N) = IN(h+1:NC)
             EXIT
          END IF
       END DO
    END IF
  END SUBROUTINE Colinear

  SUBROUTINE Perfect_Predict_Probit(D,X,CON,IND,N)
    ! CHECKS FOR VARIABLES THAT PREDICT PERFECTLY IN A PROBIT MODEL
    ! IT RETURNS THE INDEX OF N COLUMNS OF THE MATRIX IN IND(1:N) MINUS THE ONES THAT PREDICT PERFECTLY
    ! CON IS THE LOCATION OF THE CONSTANT. IF NO CONSTANT MAKE CON=0
    INTEGER, INTENT(IN) :: CON
    REAL(8), INTENT(IN) :: X(:,:),D(:)
    INTEGER, INTENT(OUT) :: IND(SIZE(X,2)),N
    INTEGER :: c,NR,NC,N_DUMMY,SUMA,CHECK(SIZE(X,1))
    NR = SIZE(X,1)
    NC = SIZE(X,2)
    N=0
    DO c = 1, NC
       IF (c==CON) THEN
          N = N + 1
          IND(N) = c
       ELSE
          CHECK=0
          WHERE (X(:,c)==0.0d0) CHECK=1
          WHERE (X(:,c)==1.0d0) CHECK=1
          IF (SUM(CHECK)==NR) THEN
             N_DUMMY = COUNT(X(:,c)==0.0d0)
             SUMA = INT(SUM(D,X(:,c)==0.0d0))
             IF ((.NOT.SUMA==0) .AND. (.NOT.SUMA==N_DUMMY)) THEN
                N_DUMMY = COUNT(X(:,c)==1.0d0)
                SUMA = INT(SUM(D,X(:,c)==1.0d0))
                IF ((.NOT.SUMA==0) .AND. (.NOT.SUMA==N_DUMMY)) THEN
                   N = N + 1
                   IND(N) = c
                END IF
             END IF
          ELSE
             N = N + 1
             IND(N) = c
          END IF
       END IF.
    END DO
  END SUBROUTINE Perfect_Predict_Probit
  
  SUBROUTINE Sort(Xin,xout,ind)
    !  Sorts xout into ascENDing order - Quicksort
    REAL (8), DIMENSION (:), INTENT (in) :: Xin
    INTEGER, INTENT(INOUT) :: ind(:)
    REAL (8), DIMENSION (SIZE(ind)), INTENT (Out) :: xout
    xout=XIN
    Call SUBSOR (xout, 1, Size (xout),ind)
    Call INSSORT (xout,ind)
  CONTAINS
    RECURSIVE SUBROUTINE SUBSOR (xout, IDEB1, IFIN1,ind)
      !  Sorts xout from IDEB1 to IFIN1
      REAL(8), DIMENSION (:), INTENT (InOut) :: xout
      INTEGER, INTENT (In) :: IDEB1, IFIN1
      INTEGER, INTENT(INOUT) :: ind(:)
      INTEGER, PARAMETER :: NINS = 16 ! Max for insertion sort
      INTEGER :: ICRS, IDEB, IDCR, IFIN, IMIL,IWRK,IPIV
      REAL(8) :: XPIV, XWRK
      IDEB = IDEB1
      IFIN = IFIN1
      !  IF we DOn't have enough values to make it worth while, we leave
      !  them unsorted, and the final insertion sort will take care of them
      IF ((IFIN - IDEB) > NINS) THEN
         IMIL = (IDEB+IFIN) / 2
         !  One chooses a pivot, median of 1st, last, and middle values
         IF (xout(IMIL) < xout(IDEB)) THEN
            XWRK = xout (IDEB)
            IWRK = IND (IDEB)
            xout (IDEB) = xout (IMIL)
            IND (IDEB) = IND (IMIL)
            xout (IMIL) = XWRK
            IND (IMIL) = IWRK
         END IF
         IF (xout(IMIL) > xout(IFIN)) Then
            XWRK = xout (IFIN)
            IWRK = IND(IFIN)
            xout (IFIN) = xout (IMIL)
            IND (IFIN) = IND (IMIL)
            xout (IMIL) = XWRK
            IND (IMIL) = IWRK
            IF (xout(IMIL) < xout(IDEB)) Then
               XWRK = xout (IDEB)
               IWRK = IND (IDEB)
               xout (IDEB) = xout (IMIL)
               IND (IDEB) = IND (IMIL)
               xout (IMIL) = XWRK
               IND (IMIL) = IWRK
            END IF
         END IF
         XPIV = xout (IMIL)
         IPIV = IND (IMIL)
         !  One exchanges values to put those > pivot in the END and
         !  those <= pivot at the beginning
         ICRS = IDEB
         IDCR = IFIN
         ECH2: DO
            DO
               ICRS = ICRS + 1
               IF (ICRS >= IDCR) Then
                  !  the first  >  pivot is IDCR
                  !  the last   <= pivot is ICRS-1
                  !  Note: IF one arrives here on the first iteration, then
                  !  the pivot is the maximum of the set, the last value is equal
                  !  to it, and one can reduce by one the size of the set to process,
                  !  as IF xout (IFIN) > XPIV
                  Exit ECH2
               END IF
               IF (xout(ICRS) > XPIV) Exit
            END DO
            DO
               IF (xout(IDCR) <= XPIV) Exit
               IDCR = IDCR - 1
               IF (ICRS >= IDCR) Then
                  !  The last value < pivot is always ICRS-1
                  Exit ECH2
               END IF
            END DO
            XWRK = xout (IDCR)
            IWRK = IND (IDCR)
            xout (IDCR) = xout (ICRS)
            IND (IDCR) = IND (ICRS)
            xout (ICRS) = XWRK
            IND (ICRS) = IWRK
         END DO ECH2
         !  One now sorts each of the two sub-intervals
         Call SUBSOR (xout, IDEB1, ICRS-1,IND)
         Call SUBSOR (xout, IDCR, IFIN1,IND)
      END IF
      RETURN
    END SUBROUTINE SUBSOR
    
    SUBROUTINE INSSORT (xout,IND)
      !  Sorts xout into increasing order (Insertion sort)
      REAL(8), DIMENSION (:), INTENT (InOut) :: xout
      INTEGER, INTENT(INOUT) :: ind(:)
      INTEGER :: ICRS, IDCR,IWRK
      REAL(8) :: XWRK
      DO ICRS = 2, Size (xout)
         XWRK = xout (ICRS)
         IWRK = IND (ICRS)
         IF (XWRK >= xout(ICRS-1)) Cycle
         xout (ICRS) = xout (ICRS-1)
         IND (ICRS) = IND (ICRS-1)
         DO IDCR = ICRS - 2, 1, - 1
            IF (XWRK >= xout(IDCR)) Exit
            xout (IDCR+1) = xout (IDCR)
            IND (IDCR+1) = IND (IDCR)
         END DO
         xout (IDCR+1) = XWRK
         IND (IDCR+1) = IWRK
      END DO
      RETURN
    END SUBROUTINE INSSORT
    
  END SUBROUTINE Sort
  
  SUBROUTINE ludcmp(a,indx,d,singular)
    IMPLICIT NONE
    REAL(8), DIMENSION(:,:), INTENT(INOUT) :: a
    INTEGER, DIMENSION(:), INTENT(OUT) :: indx
    INTEGER, INTENT(OUT) :: singular
    REAL(8), INTENT(OUT) :: d
    REAL(8), DIMENSION(size(a,1)) :: vv
    REAL(8), PARAMETER :: TINY=1.0d-20
    INTEGER :: j,n,imax
    singular=0
    n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
    d=1.0
    vv=maxval(abs(a),dim=2)
    if (any(vv == 0.0)) call nrerror('singular matrix in ludcmp')
    vv=1.0d0/vv
    do j=1,n
       imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
       if (j /= imax) then
          call swap(a(imax,:),a(j,:))
          d=-d
          vv(imax)=vv(j)
       end if
       indx(j)=imax
       if (a(j,j) == 0.0) THEN
           a(j,j)=TINY
           singular=1
        END IF
       a(j+1:n,j)=a(j+1:n,j)/a(j,j)
       a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-Outer_Product(a(j+1:n,j),a(j,j+1:n))
    end do
  END SUBROUTINE ludcmp
  
  SUBROUTINE lubksb(a,indx,b)
    IMPLICIT NONE
    REAL(8), DIMENSION(:,:), INTENT(IN) :: a
    INTEGER, DIMENSION(:), INTENT(IN) :: indx
    REAL(8), DIMENSION(:), INTENT(INOUT) :: b
    INTEGER :: i,n,ii,ll
    REAL(8) :: summ
    n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
    ii=0
    do i=1,n
       ll=indx(i)
       summ=b(ll)
       b(ll)=b(i)
       if (ii /= 0) then
          summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
       else if (summ /= 0.0) then
          ii=i
       end if
       b(i)=summ
    end do
    do i=n,1,-1
       b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
    end do
  END SUBROUTINE lubksb
    
  SUBROUTINE jacobi(ain,d,v,nrot)
    IMPLICIT NONE
    REAL(8), DIMENSION(:,:), INTENT(IN) :: ain
    REAL(8), INTENT(OUT) :: d(SIZE(AIN,1)),v(SIZE(AIN,1),SIZE(AIN,2))
    INTEGER, INTENT(OUT) :: nrot
    INTEGER :: i,j,ip,iq,n
    REAL(8) :: c,g,h,s,sm,t,tau,theta,tresh,a(SIZE(AIN,1),SIZE(AIN,2))
    REAL(8), DIMENSION(size(AIN,1)) :: b,z
    LOGICAL :: upper(SIZE(AIN,1),SIZE(AIN,2))
    n=assert_eq(size(AIN,1),size(AIN,2),'The routine only works for square symmetric matrices: jacobi') 
    A=AIN
    v=0.0d0
    upper=.FALSE.
    FORALL(i=1:n,j=1:n,i<j) upper(i,j)=.TRUE.
    FORALL(i=1:n) 
       v(i,i)=1.0D0
       b(i)=A(i,i)
    END FORALL
    d=b
    z=0.0d0
    nrot=0
    DO i=1,50
       sm=SUM(ABS(a),upper)
       IF (sm == 0.0d0) RETURN
       tresh=MERGE(0.2d0*sm/n**2,0.0d0, i < 4 )
       DO ip=1,n-1
          DO iq=ip+1,n
             g=100.0d0*abs(a(ip,iq))
             IF ((i > 4) .and. (abs(d(ip))+g == abs(d(ip))).and.(abs(d(iq))+g == abs(d(iq)))) THEN
                a(ip,iq)=0.0
             ELSEIF (abs(a(ip,iq)) > tresh) THEN
                h=d(iq)-d(ip)
                IF (abs(h)+g == abs(h)) THEN
                   t=a(ip,iq)/h
                ELSE
                   theta=0.5d0*h/a(ip,iq)
                   t=1.0d0/(abs(theta)+sqrt(1.0d0+theta**2))
                   IF (theta < 0.0) t=-t
                END IF
                c=1.0d0/dsqrt(1+t**2)
                s=t*c
                tau=s/(1.0d0+c)
                h=t*a(ip,iq)
                z(ip)=z(ip)-h
                z(iq)=z(iq)+h
                d(ip)=d(ip)-h
                d(iq)=d(iq)+h
                a(ip,iq)=0.0
                CALL jrotate(a(1:ip-1,ip),a(1:ip-1,iq))
                CALL jrotate(a(ip,ip+1:iq-1),a(ip+1:iq-1,iq))
                CALL jrotate(a(ip,iq+1:n),a(iq,iq+1:n))
                CALL jrotate(v(:,ip),v(:,iq))
                nrot=nrot+1
             END IF
          END DO
       END DO
       b=b+z
       d=b
       z=0.0
    END DO
    CALL NRERROR('too many iterations in jacobi')
  CONTAINS
    SUBROUTINE jrotate(a1,a2)
      REAL(8), DIMENSION(:), INTENT(INOUT) :: a1,a2
      REAL(8), DIMENSION(size(a1)) :: wk1
      wk1=a1
      a1=a1-s*(a2+a1*tau)
      a2=a2+s*(wk1-a2*tau)
    END SUBROUTINE jrotate
  END SUBROUTINE jacobi
  
  
  function eye(n)
    
    
        !##### INPUT/OUTPUT VARIABLES #############################################
        
        ! matrix of the system
        integer, intent(in) :: n
        
        ! the inverse of the matrix
        real*8 :: eye(n,n)
        
        ! other variable 
        integer :: i,j 
        
       Do i=1,n 
         Do  j=1,n 
           
           if (i==j) then 
               
           eye(i,j)=1.0d0
           
           else 
               
               eye(i,j)=0.0d0 
               
           end if 
           
           
           end Do            
           
           end Do 
        
  end function
  
  function kron(A,B) 
  
  ! matrix of the system
        real*8, intent(in) :: A(:,:), B(:,:) 
        
        ! the inverse of the matrix
        real*8, dimension(:), allocatable :: kron(:,:) 
        
        ! other variable 
        integer :: i,j,n1,k1,n2,k2, n,k 
        
        n1=size(A(:,1)) ! number of rows 
        k1=size(A(1,:)) ! number of columns
        
        n2=size(B(:,1)) ! number of rows 
        k2=size(B(1,:)) ! number of columns
        
        n=n1*n2 
        k=k1*k2 
        
        allocate(kron(n,k)) 
        
        Do i=1,n2 
            Do j=1,k2
              
                kron((i-1)*n1+1:i*n1,(j-1)*k1+1:j*k1)= A*B(i,j) 
        
            end Do 
            
        end Do       
        
        
  end function  
  
  subroutine exvar_mv(F,F0,lambda,Q,EY,Var) 
  
	
	implicit none
	
	! autoregression parameter
       real*8, intent(in) :: F(:,:)
    
       ! constant 
       real*8, intent(in) :: F0(:)
	   
      ! QlambdaQ' is structural variance-covariance matrix 
       real*8, intent(in) :: lambda(:)
    
       ! variance of the shock
       real*8, intent(in) :: Q(:,:)
	   
       ! longrun mean 
       real*8, intent(out) :: ey(:)
    
       ! variance-covariance 
       real*8, intent(out) :: Var(:, :)
	   
    ! other variables 
	   
	   
       integer :: i,j,kk
       REAL*8, dimension (:,:), allocatable :: QA,tempmtrx
       real*8, dimension(:), allocatable :: QB, vecSigma, lambda1
       
       
        kk=size(F0) ! number of variables         
        
        allocate(QA(kk,kk)) 
        allocate(tempmtrx(kk,kk))
        allocate(QB(kk)) 
        allocate(vecSigma(kk**2))
        allocate(lambda1(kk**2))
		
		
	    QA=matmul(matmul(transpose(Q),F),Q) 
        QB=matmul(transpose(Q),F0)
        tempmtrx=eye(kk)-QA
        EY=matmul(lu_invert(tempmtrx),QB) 
		
        ! vectorization of variance covariance matrix 
        
        lambda1=0.0d0

        Do i=1,kk 
        
        lambda1((kk+1)*(i-1)+1)= lambda(i)
        
         end Do 
            
        
         
        vecSIGMA    = matmul(matmul(lu_invert(eye(kk**2) - kron(QA,QA)),eye(kk**2)),lambda1);
        Var= reshape(vecSIGMA, (/kk, kk/));       
    
	
	end subroutine
	
  
  
  
END MODULE MATRIX
