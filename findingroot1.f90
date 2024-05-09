module rootfinding
 
 implicit none 

    contains 
    
    subroutine newton_method(x,f,fprime)
    
    integer :: niter 
    REAL(8), INTENT(INOUT) :: x
    real*8 :: threshold,residual,xNew, f, fprime
    INTEGER, PARAMETER :: maxiter=200,EPS=2.220446049250313D-9
  
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
    threshold=eps*(1+abs(x))
    
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
  
end module 
