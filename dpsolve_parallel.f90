program dpsolve_parallel

use GLOBALS
use toolbox
use MPI 
use grid_equal
use discretize_cshocks
use MATRIX
use rootfinding


implicit none

!!!! declare variables 

integer :: i,j,k,is,js,iter,ik
real*8 :: x1, x2
!real*8 :: m,step,Q(kk,kk),lambda(kk) 

CALL MPI_INIT(ierr ) ! Initialize the MPI execution environment
CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr ) !  Determines the rank of the calling process in the communicator, myid : rank of calling process 
CALL MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr ) ! Determines the size of the group associated with a communicator

call tic()
! initialize the grid point of debt 
        
        V(1,1)=V11 
        V(1,2)=V21 
        V(2,1)=V21 
        V(2,2)=V22  

        
        F(1,1)=rho11 
        F(1,2)=rho12 
        F(2,1)=rho21 
        F(2,2)=rho22
        
        F0=0.0d0         
        
! grid point of debt 
call linspace (B_l,B_h,B)        

! dicretize the continuous shocks process
call tauchen_AR2(F,F0,V,z,zz,zf,pi) 

    IF (myid==0) THEN ! process with myid==0 reads initial values   
         OPEN(6440,file='markov.txt',status="replace")
  DO i=1,NN**kk
     WRITE(6440,fmt='(16f20.8)') pi(i,:)
  END DO  
  
  DO i=1,NN**kk
     WRITE(6440,fmt='(2f20.8)') zf(i,:)
  END DO  
  
  CLOSE(6440) 
    end if 

    DO i=1,NN
     WRITE(*,fmt='(16f20.8)') z(i,:)
  END DO  
! initializa C, P, B 


do js=1,Nb+1
    do is=1,NN
        do ik=1,NN   
          c(js, is,ik) = exp(z(is,1))+rate*B(js)   
        Bp(js,is,ik)=B(js)*(1+rate)-c(js,is,ik)+exp(z(is,1))
        Pt(js,is,ik)=((1-omega)/omega)*(c(js,is,ik)/exp(z(ik,2)))**(eta+1) 
        end Do 
        end Do     
end Do 

Do iter=1,itermax

      !do is = 1, NS  
      !    
      !      call spline_interp(c(:, is), coeff_c(:, is))
      !enddo
    
    do js=1,Nb+1 
        do is = 1,NN 
            do ik=1,NN
       
                                
                Bp_temp(js,is,ik)=-kratio*(Pt(js,is,ik)*exp(z(ik,2))+exp(z(is,1)))
                c_temp(js,is,ik)=exp(z(is,1))+B(js)*(1+rate)-Bp_temp(js,is,ik)              
        end Do 
        end Do 
    end Do 

     do js=1,Nb+1 
        do is = 1,NN 
            do ik=1,NN
                icom1=is
                icom2=ik     
                
    checkU(js,is,ik)=period_utility(c_temp(js,is,ik),Bp_temp(js,is,ik))   
            end Do 
        end Do 
     end Do 
   
    do js=1,Nb+1 
        do is = 1,NN 
            do ik=1,NN       
                debt=B(js) 
                icom1=is
                icom2=ik
                x1=0.01d0
                x2=B(js)*(1+rate)+exp(zz(is,1))+20.0d0  
                
                call bisection(foc,x1,x2,1d-6,Root,flag)
                
                
                c_temp1(js,is,ik) = Root
                Bp_temp1(js,is,ik)=B(js)*(1+rate)-c_temp1(js,is,ik)+exp(z(is,1))
                                              
        end Do 
        end Do 
    end Do
     
      
      do js = 1,Nb+1
            do is = 1,NN 
                do ik=1, NN 
                
                
                
                if (checkU(js,is,ik) > 0) then ! contraint binds 
                   
                    Bp_new(js,is,ik)=Bp_temp(js,is,ik)
                    c_new(js,is,ik)=C_temp(js,is,ik)
                    Pt_new(js,is,ik)=((1-omega)/omega)*(c_new(js,is,ik)/exp(z(ik,2)))**(eta+1) 
                    else 
      
                Bp_new(js,is,ik)=Bp_temp1(js,is,ik)
                c_new(js,is,ik)=C_temp1(js,is,ik)
                Pt_new(js,is,ik)=((1-omega)/omega)*(c_new(js,is,ik)/exp(z(ik,2)))**(eta+1)
                    
                    end if 
                  
                   
    
    
                    end Do                 
      
            enddo
      enddo
      
      con_lev = maxval(abs((c_new(:,:,:) - c(:, :,:))/max(abs(c(:, :,:)), 1d-10)))
      price_lev = maxval(abs((Pt_new(:,:,:) - Pt(:, :,:))/max(abs(Pt(:, :,:)), 1d-10)))
      b_lev = maxval(abs((Bp_new(:,:,:) - Bp(:, :,:))/max(abs(Bp(:, :,:)), 1d-10)))
      
      
      !con_lev = maxval(abs(c_new(:,:,:) - c(:,:,:))/max(abs(c(:, :,:)), 1d-10))
      !price_lev =maxval(abs(Pt_new(:,:,:) - Pt(:, :,:))/max(abs(Pt(:, :,:)), 1d-10))
      !b_lev=maxval(abs(Bp_new(:, :,:) - Bp(:, :,:))/max(abs(Bp(:, :,:)), 1d-10))
      
      write(*,'(i5,2x,3f20.7)')iter, con_lev, price_lev, b_lev 
      
          

      
      criterion=max(con_lev,price_lev,b_lev) 
      
                          
        !check for convergence       
       
        c = c_new
        Pt=Pt_new 
        Bp=Bp_new
        
        if (criterion<tolerance) then  
           call toc()
           exit
        end if         
       
end do

      do js = 1,Nb+1
            do is = 1,NN 
                do ik=1, NN 
                 debt=B(js) 
                icom1=is
                icom2=ik
                
          Check_foc(js,is,ik)=foc(c_new(js,is,ik)) 
                end Do 
            end Do 
      end Do 
      
      
 IF (myid==0) THEN     
      
    OPEN(6440,file='check.txt',status="replace")
    DO i=1,Nb+1
     WRITE(6440,fmt='(16f20.8)') Check_foc(i,1,1),Check_foc(i,1,2),Check_foc(i,1,3),Check_foc(i,1,4), &
         & Check_foc(i,2,1),Check_foc(i,2,2),Check_foc(i,2,3),Check_foc(i,2,4), &
         & Check_foc(i,3,1),Check_foc(i,3,2),Check_foc(i,3,3),Check_foc(i,3,4), &
         & Check_foc(i,4,1),Check_foc(i,4,2),Check_foc(i,4,3),Check_foc(i,4,4)
    end Do 
    close(6440)
    
OPEN(6440,file='result.txt',status="replace")
    DO i=1,Nb+1
     WRITE(6440,fmt='(49f20.8)') B(i), c_new(i,1,1), c_new(i,1,2),c_new(i,1,3),c_new(i,1,4),c_new(i,2,1),c_new(i,2,2),c_new(i,2,3),c_new(i,2,4), & 
         & c_new(i,3,1), c_new(i,3,2),c_new(i,3,3),c_new(i,3,4),c_new(i,4,1),c_new(i,4,2),c_new(i,4,3),c_new(i,4,4), & 
         & Bp_new(i,1,1), Bp_new(i,1,2),Bp_new(i,1,3),Bp_new(i,1,4),Bp_new(i,2,1),Bp_new(i,2,2),Bp_new(i,2,3),Bp_new(i,2,4), &  
         & Bp_new(i,3,1), Bp_new(i,3,2),Bp_new(i,3,3),Bp_new(i,3,4),Bp_new(i,4,1),Bp_new(i,4,2),Bp_new(i,4,3),Bp_new(i,4,4), & 
         checkU(i,1,1),checkU(i,1,2),checkU(i,1,3),checkU(i,1,4),checkU(i,2,1),checkU(i,2,2),checkU(i,2,3),&
          checkU(i,2,4),checkU(i,3,1),checkU(i,3,2),checkU(i,3,3),checkU(i,3,4),checkU(i,4,1),checkU(i,4,2),checkU(i,4,3),checkU(i,4,4)
         END DO              
    close(6440) 

    
    
    
call plot(B(1:150),Bp_new(1:150,3,2), color='#a83655',legend='B')
call execplot(xlabel='Currnet Bond Holdings',legend='rs' )


call plot(B(1:150),C_new(1:150,3,2), color='#0025ad',legend='Consumption')
call execplot(xlabel='Currnet Bond Holdings',legend='rs' )
     


call plot(B(1:150),Pt_new(1:150,3,2), color='#a83655',legend='Pt')
call execplot(xlabel='Currnet Bond Holdings',legend='rs' )

end if 
     
CALL MPI_FINALIZE(ierr)

end program 

    
  