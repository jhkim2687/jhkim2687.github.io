Module discretize_cshocks 

IMPLICIT NONE 

	CONTAINS

  !##############################################################################
    ! SUBROUTINE discretize_AR
    !
    ! Discretizes an AR(1) process of the form z_j = zmu + rho*z_{j-1} + eps using
    !     the Tauchen method. Univariate variable. 
    !##############################################################################

    subroutine tauchen_AR(rho, mu, sig, z, pi)
    
    use toolbox 
    
    implicit none 
    
    !##### INPUT/OUTPUT VARIABLES #############################################
    
       ! autoregression parameter
       real*8, intent(in) :: rho
    
       ! unconditional mean of the process
       real*8, intent(in) :: mu
    
       ! st.dev of the shock
       real*8, intent(in) :: sig
    
       ! discrete shock values
       real*8, intent(out) :: z(:)
    
       ! transition matrix
       real*8, intent(out) :: pi(:, :)
    
       
       
        !##### OTHER VARIABLES ####################################################
    
       integer :: n, i,j
       real*8 :: m,step
       
        m=4d0 
        n = size(z) 
        z(1)=mu/(1-rho) - m*sqrt(sig**2/(1-rho**2))
        z(n)=mu/(1-rho) + m*sqrt(sig**2/(1-rho**2))
        step=(z(n)-z(1))/(n-1)
        
        Do i=2,(n-1) 
        
        z(i)=z(i-1)+step 
        
        end Do 
        
        Do i=1,N 
        
         Do j=1,N 
         
          if (j==1) then 
          
          pi(i,j)=normalcdf(((z(1)-mu-rho*z(i)+step/2)/sig),0.0d0,1.0d0)
          
          else if (j==n) then 
          
          pi(i,j)= 1-normalcdf(((z(n)-mu-rho*z(i)-step/2)/sig),0.0d0,1.0d0)
          
          else 
          ! i to j transition probability 
          
          pi(i,j)=normalcdf(((z(j) - mu - rho*z(i) + step/2) / sig),0.0d0,1.0d0) - &
          & normalcdf(((z(j) - mu - rho*z(i) - step/2) / sig),0.0d0,1.0d0);
         
         end if 
         
         end Do 
         
        end Do         
        
        
    
    end subroutine     
    
    
    
    subroutine tauchen_AR_G(rho, mu, sig, z, pi)
    
    ! yt = mu + rho*y(t-1) + eps 
    
    
    use toolbox 
    
    implicit none 
    
    !##### INPUT/OUTPUT VARIABLES #############################################
    
       ! autoregression parameter
       real*8, intent(in) :: rho
    
       ! unconditional mean of the process
       real*8, intent(in) :: mu
    
       ! variance of the shock
       real*8, intent(in) :: sig
    
       ! discrete shock values
       real*8, intent(in) :: z(:)
    
       ! transition matrix
       real*8, intent(out) :: pi(:, :)
    
       
       
        !##### OTHER VARIABLES ####################################################
    
       integer :: n, i,j
       real*8 :: step
       
        n = size(z) 
        step=(z(n)-z(1))/(n-1)
        
         
        
        Do i=1,N 
        
         Do j=1,N 
         
          if (j==1) then 
          
          pi(i,j)=normalcdf(((z(1)-mu-rho*z(i)+step/2)/sig),0.0d0,1.0d0)
          
          else if (j==n) then 
          
          pi(i,j)= 1-normalcdf(((z(n)-mu-rho*z(i)-step/2)/sig),0.0d0,1.0d0)
          
          else 
          ! i to j transition probability 
          
          pi(i,j)=normalcdf(((z(j) - mu - rho*z(i) + step/2) / sig),0.0d0,1.0d0) - &
          & normalcdf(((z(j) - mu - rho*z(i) - step/2) / sig),0.0d0,1.0d0);
         
         end if 
         
         end Do 
         
        end Do         
        
        
    
    end subroutine     
    
    
    subroutine tauchen_AR2(F, F0, sig, z,zz,zf,pi)
    ! yt=F0+F*y(t-1)+et 
    ! et~ normal(0,sig)
    ! et is corrlated 
    ! dicretize k shocks into n^k spte spaces (n:the number of nodes, same across the shocks) 
    ! zz is shocks state space and pi is transition matrix (markov chain) 
    
    use toolbox 
    use matrix
        
        
    implicit none 
    
    !##### INPUT/OUTPUT VARIABLES #############################################
    
       ! autoregression parameter
       real*8, intent(in) :: F(:,:)
    
       ! constant 
       real*8, intent(in) :: F0(:)
    
       ! variance of the shock
       real*8, intent(in) :: sig(:,:)
    
       ! discrete shock values
       real*8, intent(out) :: zz(:,:), z(:,:), zf(:,:) 
    
       ! transition matrix
       real*8, intent(out) :: pi(:, :)
    
       
       
        !##### OTHER VARIABLES ####################################################
    
       integer :: nd,ii,jj,kk,N,nrot,l
       REAL*8, dimension (:,:), allocatable :: Q, QA, lnp2, Var 
       real*8, dimension(:), allocatable :: step, lambda, QB, EY, vecSigma, QA1, p_temp, sigma,ztemp
       real*8 :: m, ylag
       
        m=2.5d0 
        kk=size(F0) 
        nd = size(z(:,1)) !number of nodes 
        N=size(zz(:,1))
        
        allocate(step(kk))
        allocate(Q(kk,kk)) 
        allocate(lambda(kk))
        allocate(sigma(kk))
        allocate(QA(kk,kk)) 
        allocate(QB(kk)) 
        allocate(EY(kk))         
        allocate(Var(kk,kk))
        allocate(p_temp(kk))
        allocate(ztemp(kk))
        
        
        call jacobi(sig,lambda,Q,nrot)
        
        call exvar_mv(F,F0,lambda,Q,EY,Var)
        
        QA=matmul(matmul(transpose(Q),F),Q)
        QB=matmul(transpose(Q),F0)
        
        Do ii=1,kk 
            
        z(1,ii)=EY(ii) - m*sqrt(Var(ii,ii))
        z(nd,ii)=EY(ii) + m*sqrt(Var(ii,ii))
            
        end Do 
        
        Do ii=1,kk
        step(ii)=(z(nd,ii)-z(1,ii))/(nd-1)
        end Do 
        
        Do jj=1,kk
        Do ii=2,(nd-1) 
        
        z(ii,jj)=z(ii-1,jj)+step(jj)  
        
        end Do 
        end Do             
        
        
      Do l=1,kk       
       Do ii=1,nd                         
       zz((ii-1)*nd**(kk-l)+1:ii*nd**(kk-l),l)= z(ii,l)       
       end Do
      end Do
      
      Do l=2,kk 
          Do ii=nd**(kk-l+1)+1,nd**kk
          zz(ii,l)=zz(ii-nd**(kk-l+1),l)           
          end Do
      end Do
      
   
   sigma=lambda**(.5); 
    
   Do ii=1,N 
       Do jj=1,N 
           Do l=1,kk
               
               ztemp=zz(ii,:) 
               ylag=dot_product(QA(l,:),ztemp)
               
               if (zz(jj,l)+step(l)/2 > maxval(z(:,l))) then 
                   p_temp(l)=1-normalcdf((QB(l)+zz(jj,l)-step(l)/2-ylag)/sigma(l))
           elseif (zz(jj,l)-step(l)/2 < minval(z(:,l))) then
                p_temp(l)=normalcdf((QB(l)+zz(jj,l)+step(l)/2-ylag)/sigma(l))      
                    
               else 
           p_temp(l)=normalcdf((QB(l)+zz(jj,l)+step(l)/2-ylag)/sigma(l))- normalcdf((QB(l)+zz(jj,l)-step(l)/2-ylag)/sigma(l))
               end if 
            end Do 
       
               
               pi(ii,jj)=p_temp(1);
               Do l=2,kk 
                   pi(ii,jj)=pi(ii,jj)*p_temp(l)
               end Do 
               
           end Do            
    end Do  
    
        z=transpose(matmul(Q,transpose(z)))
        
          Do l=1,kk       
       Do ii=1,nd                         
       zf((ii-1)*nd**(kk-l)+1:ii*nd**(kk-l),l)= z(ii,l)       
       end Do
      end Do
      
      Do l=2,kk 
          Do ii=nd**(kk-l+1)+1,nd**kk
          zf(ii,l)=zf(ii-nd**(kk-l+1),l)           
          end Do
      end Do
        
    end subroutine     
    
	
end Module
