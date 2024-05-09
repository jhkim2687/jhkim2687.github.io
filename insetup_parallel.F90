
Module GLOBALS   
  
IMPLICIT NONE 


! model parameters

real*8, parameter :: rho11=0.901d0, rho12=0.495d0 , rho21=-0.453d0, rho22=0.225d0
real*8, parameter :: rho1=0.53d0, rho2=0.61d0, sig1=0.058d0, sig2=0.057d0 
real*8, parameter :: V11=0.00219d0, V21=0.00162d0, V22=0.00167d0, mu1=0.0d0, mu2=0.0d0, alpha=0.99d0
real*8, parameter :: eta=0.2048d0 
real*8, parameter :: rate=0.04d0 
real*8, parameter :: sigma=2.0d0 
real*8, parameter :: kratio=0.27d0, kt=0.32d0  
real*8, parameter :: omega=0.31d0 
real*8, parameter :: beta=0.91d0, eps=1d-2
integer, parameter :: NN=4, Nb=200, kk=2, Ns=16 
INTEGER :: myid,numprocs,ierr,icom1,icom2, flag
integer, parameter :: itermax = 2000
real*8, parameter :: B_l=-1.2d0,B_h=1.5d0 
real*8 :: zz(Ns,kk), zf(Ns,kk), pi(Ns,Ns), B(Nb+1), Bp(Nb+1,NN,NN),C(Nb+1,NN,NN),Pt(Nb+1,NN,NN), V(kk,kk), mu(kk), F0(kk), F(kk,kk), z(NN,kk) 
real*8 :: C_temp(Nb+1,NN,NN), Bp_temp(Nb+1,NN,NN), C_new(Nb+1,NN,NN), Bp_new(Nb+1,NN,NN), Pt_new(Nb+1,NN,NN), c_in,con_lev,debt,checkU(Nb+1,NN,NN), criterion, price_lev, b_lev 
real*8 :: C_temp1(Nb+1,NN,NN), Bp_temp1(Nb+1,NN,NN), Check_foc(Nb+1,NN,NN), Root 
real*8, parameter :: tolerance = 1d-6
logical :: check

    contains 
    
    function foc(c_in) 
    
    use toolbox
    
    implicit none
    
        real*8, intent(in) :: c_in
        real*8 :: foc,bplus,cplus  
        real*8 :: yt,yn,ctotal,utility, ctotal_p, utility_p, varphi, varphiz1, varphiz2
        integer :: ip, il, ir, iq, izl1, izr1, izl2, izr2
        
        
        yt=exp(z(icom1,1))
        yn=exp(z(icom2,2)) 
        bplus=yt+debt*(1+rate)-c_in 
        
        ctotal=(omega*c_in**(-eta)+(1-omega)*(yn)**(-eta))**(-1/eta) 
        utility=(ctotal**(-sigma))*((ctotal**(-eta))**((-1/eta)-1))*omega*c_in**(-eta-1) 
        
         foc = 0d0
         
        do ip = 1,NN
            do iq=1,NN  
            
            yt=exp(z(ip,1))
            yn=exp(z(iq,2)) 
            
           if (bplus<=B_h .AND. bplus>=B_l) then
            call linint_Equi(bplus, B_l, B_h, Nb, il, ir, varphi)       
            il=il+1 
            ir=ir+1
            cplus=varphi*c(il,ip,iq)+(1d0-varphi)*c(ir,ip,iq)
            
           elseif (bplus > B_h) then
               
            cplus = c(Nb+1,ip,iq) + (c(Nb+1,ip,iq)-c(Nb,ip,iq))/(B(Nb+1)-B(Nb))*(bplus-B_h)
            
           else  
            
               cplus =  c(1,ip,iq)
                !
                 end if 
                
                        
            ctotal_p=(omega*cplus**(-eta)+(1-omega)*yn**(-eta))**(-1/eta) 
            utility_p=(ctotal_p**(-sigma))*((ctotal_p**(-eta))**((-1/eta)-1))*omega*cplus**(-eta-1) 
            foc = foc + pi(NN*(icom1-1)+icom2,NN*(ip-1)+iq)*utility_p 
            enddo
        end Do 
        
        
    foc = utility - (beta*(1+rate))*foc
    end function 
    
    
    function period_utility(x_in,b_in) 
    
    use toolbox
    
    implicit none
    
        real*8, intent(in) :: x_in,b_in
        real*8 :: period_utility  
        real*8 :: yt,yn,ctotal,utility, ctotal_p, utility_p, cplus, varphi,varphiz1,varphiz2
        integer :: ip, il, ir,izl1, izr1, izl2, izr2, iq  
        
        
        yt=exp(z(icom1,1))
        yn=exp(z(icom2,2))        
         
            
        ctotal=(omega*x_in**(-eta)+(1-omega)*(yn)**(-eta))**(-1/eta) 
        utility=(ctotal**(-sigma))*((ctotal**(-eta))**((-1/eta)-1))*omega*x_in**(-eta-1) 
        
         period_utility   = 0.0d0
         
        do ip = 1,NN
            do iq=1,NN   
            
            yt=exp(z(ip,1))
            yn=exp(z(iq,2)) 
            
            ! determine the way of interporlation outside of support 
            
            
            if (b_in<=B_h .AND. b_in>=B_l) then
            call linint_Equi(b_in, B_l, B_h, Nb, il, ir, varphi)   
            il=il+1 
            ir=ir+1
            cplus=varphi*c(il,ip,iq)+(1d0-varphi)*c(ir,ip,iq)
            
           elseif (b_in > B_h) then
               
            cplus = c(Nb+1,ip,iq) + (c(Nb+1,ip,iq)-c(Nb,ip,iq))/(B(Nb+1)-B(Nb))*(b_in-B_h)
           else  
            
                cplus =  c(1,ip,iq) 
           end if 
           
            ctotal_p=(omega*cplus**(-eta)+(1-omega)*yn**(-eta))**(-1/eta) 
            utility_p=(ctotal_p**(-sigma))*((ctotal_p**(-eta))**((-1/eta)-1))*omega*cplus**(-eta-1) 
            period_utility   = period_utility + pi(NN*(icom1-1)+icom2,NN*(ip-1)+iq)*utility_p 
            enddo
            
            end Do 
        
    period_utility   = utility - (beta*(1+rate))*period_utility
        
   
        
        
    end function 
    
end Module 

