Module Gradient 

implicit none 

contains 

FUNCTION GRADIENT1(x0,f0,func)
    ! Routine to get the numerical gradient of a function using forward differences
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x0(:),f0
    INTERFACE
    REAL(8) FUNCTION func(theta)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: theta(:)
       END FUNCTION func
    END INTERFACE
    REAL(8), PARAMETER :: EPS=1.49011611938477D-08 !SQRT(NR_EPS)
    REAL(8) :: GRADIENT1(SIZE(x0)),grdd(SIZE(x0)),dh(SIZE(x0)),xdh(SIZE(x0))
    REAL(8) :: arg(SIZE(x0),SIZE(x0))
    INTEGER :: i,k
    k=SIZE(x0)
    grdd = 0.0D0
    ! Computation of stepsize (dh) for gradient
    dh = EPS*ABS(x0)
    WHERE (dh<1.0d-13) dh=EPS
    xdh = x0+dh
    dh=xdh-x0
    arg = SPREAD(x0,DIM=2,NCOPIES=k)
    DO i = 1, k
       arg(i,i)=xdh(i)    
       grdd(i)=func(arg(:,i))
    END DO
    GRADIENT1 = (grdd-f0)/dh
  END FUNCTION GRADIENT1
  
  FUNCTION GRADIENT2(x0,func)
    ! Routine to get the numerical gradient of a function using central differences
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x0(:)
    INTERFACE
    REAL(8) FUNCTION func(theta)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: theta(:)
       END FUNCTION func
    END INTERFACE
    REAL(8), PARAMETER :: EPS=6.05545445239334D-06 !NR_EPS**(1.0d0/3.0d0)
    REAL(8) :: GRADIENT2(SIZE(x0)),grdd(SIZE(x0)),dh(SIZE(x0)),xdh1(SIZE(x0)),xdh0(SIZE(x0)),dhh(SIZE(x0))
    REAL(8) :: arg(SIZE(x0),SIZE(x0))
    INTEGER :: i,k
    k=SIZE(x0)
    grdd = 0.0D0
    ! Computation of stepsize (dh) for gradient
    dh = EPS*ABS(x0)
    WHERE (dh==0.0D0) dh=EPS
    xdh1 = x0+dh
    arg = SPREAD(x0,DIM=2,NCOPIES=k)
    DO i = 1, k
       arg(i,i)=xdh1(i)    
       grdd(i)=func(arg(:,i))
    END DO
    xdh0 = x0-dh
    dhh = xdh1-xdh0
    arg = SPREAD(x0,DIM=2,NCOPIES=k)
    DO i = 1, k
       arg(i,i)=xdh0(i)    
       grdd(i)=grdd(i)-func(arg(:,i))
    END DO
    GRADIENT2 = grdd/dhh
  END FUNCTION GRADIENT2

  FUNCTION GRADIENT3(x0,func)
    ! Routine to get the numerical gradient of a function using a 5 point stencil
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x0(:)
    INTERFACE
    REAL(8) FUNCTION func(theta) !This is not the likelihood itself, it just let the program know
         IMPLICIT NONE
         REAL(8), INTENT(IN) :: theta(:)
       END FUNCTION func
    END INTERFACE
    REAL(8), PARAMETER :: EPS=6.05545445239334D-06!NR_EPS**(1.0d0/3.0d0)
    REAL(8) :: GRADIENT3(SIZE(x0)),grdd(SIZE(x0)),dh(SIZE(x0)),xdh(SIZE(x0))
    REAL(8) :: arg(SIZE(x0),SIZE(x0))
    INTEGER :: i,k
    k=SIZE(x0)
    grdd = 0.0D0
    ! Computation of stepsize (dh) for gradient
    dh = EPS*ABS(x0)
    WHERE (dh==0.0D0) dh=EPS
    xdh = x0+dh
    dh = xdh-x0
    arg = SPREAD(x0,DIM=2,NCOPIES=k)
    DO i = 1, k
       arg(i,i)=xdh(i)    
       grdd(i)=func(arg(:,i))
    END DO
    xdh = x0-dh
    dh = x0-xdh
    arg = SPREAD(x0,DIM=2,NCOPIES=k)
    DO i = 1, k
       arg(i,i)=xdh(i)    
       grdd(i)=grdd(i)-func(arg(:,i))
    END DO
    GRADIENT3 = grdd/(2.0d0*dh)
  END FUNCTION GRADIENT3

  end Module 