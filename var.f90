Module longrun_VAR

	implicit none 
	
	contains 
	
	subroutine exvar_mv(F,F0,lambda,Q,EY,V) 
	
	use matrix 
	
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
       real*8, intent(out) :: V(:, :)
	   
    ! other variables 
	   
	   
       integer :: i,j,kk
       REAL*8, dimension (:,:), allocatable :: QA,kron
       real*8, dimension(:), allocatable :: QB, vecSigma, QA1
       
       
        kk=size(F0) ! number of variables         
        
        allocate(QA(kk,kk)) 
        allocate(QB(kk)) 
        allocate(kron(kk**2,kk**2)) 
        allocate(vecSigma(kk**2)
        allocate(QA1(kk**2)
		
		
	    QA=matmul(matmul(transpose(Q),F),Q) 
        QB=matmul(transpose(Q),F0)
        EY=matmul(transpose(Q), matmul(lu_invert(eye(kk)-F),F0)) 
		
		
        QA1=reshape(QA,(/kk**2/)
        vecSIGMA    = matmul(matmul(lu_invert(eye(kk**2) - kron(QA, QA)),eye(kk**2)),QA1);
        V= reshape(vecSIGMA, (/k, k/));       
    
	
	end subroutine
	
	
	end Module
