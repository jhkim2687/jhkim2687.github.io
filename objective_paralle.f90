Module Objective 
    
    use GLOBALS  
    

IMPLICIT NONE 

CONTAINS 

REAL(8) FUNCTION obj(theta)

    use MPI 
    USE toolbox, only : lu_invert
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: theta(:)
    REAL(8) :: thetastar(N_COV),yfake(N_OBS),dist(N_COV),tempy, XXF(N_Cov,N_COV),XYF(N_COV),randomnum,tp,A,B,thetahat(N_cov), & 
                & btemp(N_cov),thetatemp(N_cov),mybtemp(N_cov), obj_temp
	Integer :: i, j
	
     
    thetahat=betahat(1:N_cov) 
    
    mybtemp=0.0d0
    btemp=0.0d0 
    
    
    Do j=myid+1,S,numprocs 
        
        Do i=1,N_OBS
            
            tempy=dot_product(x(i,:),theta)+ u(j,i)
            tp=tempy/0.009
		    A =exp(tp)
		    B =1+A
            yfake(i)=A/B
         
        end Do
        
        XXF=matmul(transpose(X),X)  
        XYF=matmul(transpose(X),yfake) 
        thetatemp=matmul(lu_invert(XXF),XYF) 
        mybtemp=mybtemp+thetatemp
        
    end do     
    
                     
    CALL MPI_REDUCE(mybtemp,btemp,N_COV,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr) ! Perform the partial sum
    CALL MPI_BCAST(btemp,N_COV,MPI_REAL8,0,MPI_COMM_WORLD,ierr) ! Make sure everyone has it
    
    
    thetastar=btemp/S
	dist=thetastar-thetahat
    obj_temp=dot_product(dist,dist) 
	CALL MPI_BCAST(obj_temp,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr) ! Make sure everyone has it
    
    obj=obj_temp 
    
	end Function 
	
end 

