module rootfinding
 
 implicit none 

    contains 
    
    subroutine newton_method(x,func,fderiv)
    
    integer :: niter 
    REAL(8), INTENT(INOUT) :: x
    real*8 :: threshold,residual,xNew, f, fprime
    INTEGER, PARAMETER :: maxiter=500,EPS=2.220446049250313D-9
  
    INTERFACE
       REAL(8) FUNCTION func(x) ! 
       IMPLICIT NONE 
         REAL(8), INTENT(IN) :: x
       END FUNCTION func
    END INTERFACE
	   
	   
    INTERFACE   
       REAL(8) FUNCTION fderiv(x) ! This is not the liklihood 
       IMPLICIT NONE 
         REAL(8), INTENT(IN) :: x
       END FUNCTION fderiv
    END INTERFACE	
	
    
    residual=1.0d0
    niter=1 
    threshold=eps
    
    Do while ((residual > threshold) .and. (nIter < maxIter))
        
    ! compute function value evaluated at x
    f=func(x)
    fprime=fderiv(x)

    !! Exit if f' is near or become zero
    if (abs(fprime) < 1.e-12) then
       print *, '[Error: newton_method] Function derivative becomes very close to zero or zero.'
       print *, 'f=',f, 'df/dx =',fprime
       print *, 'Aborting now in order to avoid division by zero.'
       stop
    end if

    !! Algorithm
    xNew = x - f/fprime
    residual=abs(xNew-x)
    
    x = xNew
    nIter = nIter + 1
    threshold=eps*(1+abs(x)) 
      
    !print*, niter,x,f,fprime
    end do 
    

end subroutine

 Subroutine bisection(f,x1,x2,eps,Root,flag)
!============================================================
! Solutions of equation f(x)=0 on [x1,x2] interval
! Method: Bisectional (closed domain) (a single root)
! Alex G. January 2010
!------------------------------------------------------------
! input ...
! f   - function - evaluates f(x) for any x in [x1,x2]
! x1  - left endpoint of initial interval
! x2  - right endpoint of initial interval
! eps - desired uncertainity of the root as |b-a|<eps
! output ...
! Root  - root of the equation f(x)=0
! flag  - indicator of success
!         >0 - a single root found, flag=number of iterations
!          0 - no solutions for the bisectional method
!         <0 - not a root but singularity, flag=number of iterations
!
! Comments: Function f(x) has to change sign between x1 and x2
!           Max number of iterations - 200 (accuracy (b-a)/2**200)
!====================================================================
implicit none

real*8, intent(in) :: x1, x2, eps 
real*8, intent(out) :: Root
real*8 :: a, b, c
integer :: i, flag
integer, parameter:: iter=200

INTERFACE
       REAL(8) FUNCTION f(x) ! 
       IMPLICIT NONE 
         REAL(8), INTENT(IN) :: x
       END FUNCTION f
    END INTERFACE

!* check the bisection condition
if(f(x1)*f(x2)>0.0) then
  flag = 0
  return
end if

!* initialize calculations
a=x1
b=x2

!* Iterative refining the solution 
do i=1,iter
  c=(b+a)/2.0
  if(f(a)*f(c).le.0.0) then
      b = c
    else
      a=c
  end if
! condition(s) to stop iterations)
  if(abs(b-a)<= eps) exit  
end do
Root=(b+a)/2.0

!* check if it is a root or singularity
if (abs(f(Root)) < 1.0) then
  flag=i
  else
  flag = -i
end if
end subroutine bisection

end Module 
