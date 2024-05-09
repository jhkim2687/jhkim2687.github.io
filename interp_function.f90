MODULE interp_funcs

use toolbox

IMPLICIT NONE 

    CONTAINS 
    
    
 function three_sinterp(x1,x2, x3, yyy, xx1, xx2, xx3)
!
! Nk Nb Ns 
! linear interpolation over Nb 
!
    REAL(8), INTENT(IN) :: x1(:), x2(:), x3(:), yyy(:,:,:), xx1, xx2, xx3   
    real*8 :: varphi, three_sinterp 
    REAL(8), dimension (:,:), allocatable :: ff, coeff
    real*8, dimension(:), allocatable :: coeff1, ff2 
    integer :: il, ir 
    integer :: x1size, x2size, x3size
	integer :: i, j

    x1size = SIZE(x1,1)
	x2size = SIZE(x2,1)
    x3size = SIZE(x3,1)
    
    
    allocate (ff(x1size,x3size))  ! Nk by Ns matrix   
    allocate (coeff(x1size,x3size+2))  ! Nk by Ns matrix   
     allocate (coeff1(x2size+2)) 
      allocate (ff2(x2size)) 
    
    if (xx2.LE.x2(x2size) .AND. xx2.GE.x2(1)) then 
        
        Do i=1,x1size 
            Do j=1,x3size
        
         ff2 = yyy(i,:,j)       
        call spline_interp(ff2,coeff1)
        
               
             ff(i,j) = spline_eval(xx2, coeff1, x2(1), x2(x2size))        
        
            endDo 
        endDo 
        
      
     elseif (xx2 > x2(x2size)) then            
                 
                 

    ff = ((yyy(:,x2size,:)-yyy(:,x2size-1,:))/ (x2(x2size)-x2(x2size-1)) )  * (xx2-x2(x2size))  + yyy(:,x2size,:)          
                     
     
     else 
     
     ff = ((yyy(:,2,:)-yyy(:,1,:))/ (x2(2)-x2(1)) )  * (xx2-x2(1))  + yyy(:,1,:)               
         
     endif 
     
     
    call spline_coeff2(x1,x3, ff, coeff)
    
    three_sinterp = two_sinterp(x1,x3, ff, xx1, xx3)
    


end function    

function three_interp(x1,x2, x3, yyy, xx1, xx2, xx3)
!
! Nk Nb Ns 
! linear interpolation over Nb 
!
    REAL(8), INTENT(IN) :: x1(:), x2(:), x3(:), yyy(:,:,:), xx1, xx2, xx3   
    real*8 :: varphi, three_interp 
    REAL(8), dimension (:,:), allocatable :: ff, coeff 
    integer :: il, ir 
    integer :: x1size, x2size, x3size
	integer :: i, j

    x1size = SIZE(x1,1)
	x2size = SIZE(x2,1)
    x3size = SIZE(x3,1)
    
    
    allocate (ff(x2size,x3size))  ! Nk by Ns matrix   
    allocate (coeff(x2size,x3size+2))  ! Nk by Ns matrix   
    
    if (xx1.LE.x1(x1size) .AND. xx1.GE.x1(1)) then 
        
               call linint_Equi(xx1, x1(1), x1(x1size), x1size-1, il, ir, varphi)
            
             ff = yyy(il+1,:,:)*varphi + yyy(ir+1,:,:)*(1d0-varphi)          
        
      
     elseif (xx1 > x1(x1size)) then            
                 
                 

    ff = ((yyy(x1size,:,:)-yyy(x1size-1,:,:))/ (x1(x1size)-x1(x1size-1)) )  * (xx1-x1(x1size))  + yyy(x1size,:,:)          
                     
     
     else 
     
    ff = ((yyy(2,:,:)-yyy(1,:,:))/ (x1(2)-x1(1)) )  * (xx1-x1(1))  + yyy(1,:,:)           
         
     endif 
     
     
    call spline_coeff2(x2,x3, ff, coeff)
    
    three_interp = two_sinterp(x2,x3, ff, xx2, xx3)
    


end function 

    
subroutine spline_coeff2(x1,x2, yy, coeff_f)
! coeff_f is x1size x2size+2 
use toolbox

    REAL(8), INTENT(IN) :: x1(:), x2(:), yy(:,:)
    real*8 :: coeff_f(:,:)
    !REAL(8) :: three_interp
	REAL(8), dimension (:), allocatable :: yytemp, coeff,  yy_over2 
    integer :: x1size, x2size
	integer :: i, j
    
	! interpolation over x1 x2  
	! spline interpolation over x2 
	! and linear interpolation over x1 	
	
	
	x1size = SIZE(x1,1)
	x2size = SIZE(x2,1)
    
  allocate (coeff(x2size+2))	
  allocate (yytemp(x2size))
  allocate (yy_over2(x1size))
  
	
	
	Do i=1,x1size 
    yytemp = yy(i,:)
    call spline_interp(yytemp, coeff)
    coeff_f(i,:) = coeff 
    end Do 	
         
end subroutine 


function two_interp(x2,x3, yyy, xx2, xx3)

use toolbox

    REAL(8), INTENT(IN) :: x2(:), x3(:), yyy(:,:), xx2, xx3 
    REAL(8) :: two_interp, templ, tempr, templ1, templ2, tempr1, tempr2, ttempl, ttempr
	REAL(8), dimension (:), allocatable :: coeff, ff    
	REAL(8) :: varphi2 
	integer :: x2size, x3size, il2, ir2, il, ir 
	integer :: i, j, k 
    
	! interpolation over x1 x2 x3 
	! spline interpolation over x3 
	! and linear interpolation over x1 x2 	
	
	
	x2size = SIZE(x2,1)
	x3size = SIZE(x3,1)
    
 allocate (coeff(x3size+2))
  allocate (ff(x3size))
 
	          
if (xx3 .LE. x3(x3size) .AND. xx3.GE.x3(1)) then   
    
             
    if (xx2.LE.x2(x2size) .AND. xx2.GE.x2(1)) then 
                 
            call linint_Equi(xx2, x2(1), x2(x2size), x2size-1, il2, ir2, varphi2)
            
            il2=il2+1 
            ir2=ir2+1          
    
                !WRITE(*,fmt='(i5,i5)') il2, ir2 
 
   ff = yyy(il2,:)*varphi2 + yyy(ir2,:)*(1-varphi2)             
                     
                        
                 
                 
         elseif (xx2 > x2(x2size)) then            
                 
                 

    ff = ((yyy(x2size,:)-yyy(x2size-1,:))/ (x2(x2size)-x2(x2size-1)) )  * (xx2-x2(x2size))  + yyy(x2size,:)          
                     
                     
                           
                 
                 
             else 
                 
       
                             
    ff = ((yyy(2,:)-yyy(1,:))/ (x2(2)-x2(1)) )  * (xx2-x2(1))  + yyy(1,:)          
   
                       
  

             end if            
             
             
             
             
   
             call spline_interp(ff,coeff)
             
             two_interp =  spline_eval(xx3, coeff, x3(1), x3(x3size))
                     
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
             
    elseif (xx3 > x3(x3size)) then 

                 
    if (xx2.LE.x2(x2size) .AND. xx2.GE.x2(1)) then 
                 
            call linint_Equi(xx2, x2(1), x2(x2size), x2size-1, il2, ir2, varphi2)
            
            il2=il2+1 
            ir2=ir2+1          
    
                 
 
   ff = yyy(il2,:)*varphi2 + yyy(ir2,:)*(1-varphi2)             
                     
                        
                 
                 
         elseif (xx2 > x2(x2size)) then            
                 
                 

    ff = ((yyy(x2size,:)-yyy(x2size-1,:))/ (x2(x2size)-x2(x2size-1)) )  * (xx2-x2(x2size))  + yyy(x2size,:)          
                     
                     
                           
                 
                 
             else 
                 
       
                             
    ff = ((yyy(2,:)-yyy(1,:))/ (x2(2)-x2(1)) )  * (xx2-x2(1))  + yyy(1,:)          
   
                       
  

             end if 
   
                 
                 
           two_interp =  ((ff(x3size)-ff(x3size-1))/ (x3(x3size)-x3(x3size-1)) )  * (xx3-x3(x3size))  + ff(x3size)
           
             
        
             
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
                  
         else  ! else x3<x3(1)
         
         
         
                 if (xx2.LE.x2(x2size) .AND. xx2.GE.x2(1)) then 
                 
            call linint_Equi(xx2, x2(1), x2(x2size), x2size-1, il2, ir2, varphi2)
            
            il2=il2+1 
            ir2=ir2+1          
    
                 
 
   ff = yyy(il2,:)*varphi2 + yyy(ir2,:)*(1-varphi2)             
                     
                        
                 
                 
         elseif (xx2 > x2(x2size)) then            
                 
                 

    ff = ((yyy(x2size,:)-yyy(x2size-1,:))/ (x2(x2size)-x2(x2size-1)) )  * (xx2-x2(x2size))  + yyy(x2size,:)          
                     
                     
                           
                 
                 
             else 
                 
       
                             
    ff = ((yyy(2,:)-yyy(1,:))/ (x2(2)-x2(1)) )  * (xx2-x2(1))  + yyy(1,:)          
   
                       
  

             end if 
   
                 
                 
           two_interp =  ((ff(2)-ff(1))/ (x3(2)-x3(1)) )  * (xx3-x3(1))  + ff(1)
           
           endif 

end function 



function one_interp(x2, yyy, xx2)

use toolbox

    REAL(8), INTENT(IN) :: x2(:),yyy(:), xx2 
    real*8, dimension(:), allocatable :: coeff
    REAL(8) :: one_interp
	REAL(8) :: varphi 
    integer :: x2size
	integer :: il, ir
   
    x2size = SIZE(x2,1)
    call linint_Equi(xx2, x2(1), x2(x2size), x2size-1, il, ir, varphi)
    il=il+1 
    ir=ir+1          
   one_interp = varphi*yyy(il) + (1d0-varphi)*yyy(ir)                 
   

end function 






function one_sinterp(x2, yyy, xx2in)

use toolbox

    REAL(8), INTENT(IN) :: x2(:),yyy(:), xx2in
    real*8, dimension(:), allocatable :: coeff(:)
    real*8 :: xx2 
    REAL(8) :: one_sinterp
	integer :: x2size
	
    xx2 = xx2in
    
    x2size = SIZE(x2,1)
    
    allocate(coeff(x2size+2))
    
	! linear interpolation over x2
             
    if (xx2.LE.x2(x2size) .AND. xx2.GE.x2(1)) then 
                 
                      
    call spline_interp(yyy,coeff)
       one_sinterp = spline_eval(xx2, coeff, x2(1), x2(x2size))                 
                        
                 
                 
         elseif (xx2 > x2(x2size)) then            
                 
               

  
        one_sinterp =  (yyy(x2size) - yyy(x2size-1))/(x2(x2size)-x2(x2size-1))*(xx2-x2(x2size)) + yyy(x2size)
                     
                     
                           
                 
                 
             else 
                 
       
                             
       one_sinterp =  (yyy(2) - yyy(1))/(x2(2)-x2(1))*(xx2-x2(1)) + yyy(1)
                     
   


             end if 
             
   

end function 


function two_sinterp(x2,x3, yyy, xx2, xx3)

use toolbox

    REAL(8), INTENT(IN) :: x2(:), x3(:), yyy(:,:), xx2, xx3
    REAL(8) :: two_sinterp
	REAL(8), dimension (:), allocatable :: coeff3, coeff2, ff    
    integer :: x2size, x3size 
	integer :: i, j, k 
    
! two dimension interpolation over 

    
	x2size = SIZE(x2,1)
	x3size = SIZE(x3,1)
    
 allocate (coeff3(x3size+2))
 allocate (ff(x3size))
 allocate (coeff2(x2size+2))
 
	          
if (xx3 .LE. x3(x3size) .AND. xx3.GE.x3(1)) then   
    
             
    if (xx2.LE.x2(x2size) .AND. xx2.GE.x2(1)) then 
    
    
    
          Do i=1,x3size 
          
              
              call spline_interp(yyy(:,i), coeff2)   
              
              ff(i) = spline_eval(xx2, coeff2, x2(1), x2(x2size))
          
          
          end Do      
                              
                        
                 
                 
         elseif (xx2 > x2(x2size)) then            
                 
                 

    ff = ((yyy(x2size,:)-yyy(x2size-1,:))/ (x2(x2size)-x2(x2size-1)) )  * (xx2-x2(x2size))  + yyy(x2size,:)
                     
                     
                           
                 
                 
             else 
                 
       
                             
      
      ff = ((yyy(2,:)-yyy(1,:))/ (x2(2)-x2(1)) )  * (xx2-x2(1))  + yyy(1,:)
  
                     
   


             end if 
             
             call spline_interp(ff, coeff3)                               
                     
             two_sinterp =  spline_eval (xx3, coeff3, x3(1), x3(x3size))  
             
                 
                     
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
             
    elseif (xx3 > x3(x3size)) then 

                 
    if (xx2.le.x2(x2size) .AND. xx2.ge.x2(1)) then 
                 
           Do i=1,x3size 
          
              
              call spline_interp(yyy(:,i), coeff2)   
              
              ff(i) = spline_eval(xx2, coeff2, x2(1), x2(x2size))
          
          
          end Do           
          
                
                 
                 
    elseif (xx2 > x2(x2size)) then                  

     ff = ((yyy(x2size,:)-yyy(x2size-1,:))/ (x2(x2size)-x2(x2size-1)) )  * (xx2-x2(x2size))  + yyy(x2size,:)
                 
                 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11                     
                 
                 
    else 
        ! (xx2 .LE. x2(1))                
    ff = ((yyy(2,:)-yyy(1,:))/ (x2(2)-x2(1)) )  * (xx2-x2(1))  + yyy(1,:)
       
    
     endif                     
   
                 
                 
           two_sinterp =  ((ff(x3size)-ff(x3size-1)) / (x3(x3size)-x3(x3size-1)) )  * (xx3-x3(x3size))  + ff(x3size)     
         
             
        
             
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
                  
         else  ! else x3<x3(1)
         
         
         if (xx2.le.x2(x2size) .AND. xx2.ge.x2(1)) then 
                 
          Do i=1,x3size 
          
              
              call spline_interp(yyy(:,i), coeff2)   
              
              ff(i) = spline_eval(xx2, coeff2, x2(1), x2(x2size))
          
          
          end Do             
          
                
                 
                 
    elseif (xx2 > x2(x2size)) then                  

     ff = ((yyy(x2size,:)-yyy(x2size-1,:))/ (x2(x2size)-x2(x2size-1)) )  * (xx2-x2(x2size))  + yyy(x2size,:)
                 
                 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11                     
                 
                 
    else 
        ! (xx2 .LE. x2(1))                
    ff = ((yyy(2,:)-yyy(1,:))/ (x2(2)-x2(1)) )  * (xx2-x2(1))  + yyy(1,:)
       
    
     endif          
   
                 
                 
         two_sinterp =  ((ff(2)-ff(1))/ (x3(2)-x3(1)) )  * (xx3-x3(1))  + ff(1)    
         
             
         endif   

end function 



function one_sinterp2(x2, yyy, xx2)

use toolbox

    REAL(8), INTENT(IN) :: x2(:),yyy(:), xx2
    real*8, dimension(:), allocatable :: coeff
    REAL(8) :: one_sinterp2
	integer :: x2size
	
    
    x2size = SIZE(x2,1)
    allocate (coeff(x2size+2))
    
	! linear interpolation over x2
             
    if (xx2.LE.x2(x2size) .AND. xx2.GE.x2(1)) then 
                 
                      
       call spline_interp(yyy,coeff) 
       
       one_sinterp2 = spline_eval(xx2, coeff, x2(1), x2(x2size))                 
                        
                 
                 
         elseif (xx2 > x2(x2size)) then            
                 
               

  
        one_sinterp2 =  (yyy(x2size) - yyy(x2size-1))/(x2(x2size)-x2(x2size-1))*(xx2-x2(x2size)) + yyy(x2size)
                     
                     
                           
                 
                 
             else 
                 
       
                             
       one_sinterp2 =  (yyy(2) - yyy(1))/(x2(2)-x2(1))*(xx2-x2(1)) + yyy(1)
                     
   


             end if 
             
                 
                     
  

end function 


function one_sinterp3(x2, yyy, coeff, xx2)

use toolbox

    REAL(8), INTENT(IN) :: x2(:),yyy(:), coeff(:), xx2
    REAL(8) :: one_sinterp3
	integer :: x2size
	
    
    x2size = SIZE(x2,1)
    
	! linear interpolation over x2
             
    if (xx2.LE.x2(x2size) .AND. xx2.GE.x2(1)) then 
                 
                      
    
       one_sinterp3 = spline_eval(xx2, coeff, x2(1), x2(x2size))                 
                        
                 
                 
         elseif (xx2 > x2(x2size)) then            
                 
               

  
        one_sinterp3 =  (yyy(x2size) - yyy(x2size-1))/(x2(x2size)-x2(x2size-1))*(xx2-x2(x2size)) + yyy(x2size)
                     
                     
                           
                 
                 
             else 
                 
       
                             
       one_sinterp3 =   yyy(1)
                     
   


             end if 
             
                 
                     
  

end function 


function two_linterp(x2,x3, yyy, xx2, xx3)

use toolbox

    REAL(8), INTENT(IN) :: x2(:), x3(:), yyy(:,:), xx2, xx3 
    REAL(8) :: two_linterp, templ, tempr, templ1, templ2, tempr1, tempr2, ttempl, ttempr
	REAL(8), dimension (:), allocatable :: coeff, ff    
	REAL(8) :: varphi2
	integer :: x2size, x3size, il2, ir2, il, ir 
	integer :: i, j, k 
    
	! interpolation over x1 x2 x3 
	! spline interpolation over x3 
	! and linear interpolation over x1 x2 	
	
	
	x2size = SIZE(x2,1)
	x3size = SIZE(x3,1)
    
 allocate (coeff(x3size+2))
  allocate (ff(x3size))
 
	          
if (xx3 .LE. x3(x3size) .AND. xx3.GE.x3(1)) then   
    
             
    if (xx2.LE.x2(x2size) .AND. xx2.GE.x2(1)) then 
                 
            call linint_Equi(xx2, x2(1), x2(x2size), x2size-1, il2, ir2, varphi2)
            
            il2=il2+1 
            ir2=ir2+1          
    
                 
 
   ff = yyy(il2,:)*varphi2 + yyy(ir2,:)*(1-varphi2)             
                     
                        
                 
                 
         elseif (xx2 > x2(x2size)) then            
                 
                 

    ff = ((yyy(x2size,:)-yyy(x2size-1,:))/ (x2(x2size)-x2(x2size-1)) )  * (xx2-x2(x2size))  + yyy(x2size,:)          
                     
                     
                           
                 
                 
             else 
                 
       
                             
    ff = ((yyy(2,:)-yyy(1,:))/ (x2(2)-x2(1)) )  * (xx2-x2(1))  + yyy(1,:)          
   
                       
  

             end if 
             
   
             call linint_Equi(xx3, x3(1), x3(x3size), x3size-1, il2, ir2, varphi2)
            
            il2=il2+1 
            ir2=ir2+1          
    
                 
   
             two_linterp =  ff(il2)*varphi2 + ff(ir2)*(1d0-varphi2) 
                     
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
             
    elseif (xx3 > x3(x3size)) then 

                 
    if (xx2.LE.x2(x2size) .AND. xx2.GE.x2(1)) then 
                 
            call linint_Equi(xx2, x2(1), x2(x2size), x2size-1, il2, ir2, varphi2)
            
            il2=il2+1 
            ir2=ir2+1          
    
                 
 
   ff = yyy(il2,:)*varphi2 + yyy(ir2,:)*(1-varphi2)             
                     
                        
                 
                 
         elseif (xx2 > x2(x2size)) then            
                 
                 

    ff = ((yyy(x2size,:)-yyy(x2size-1,:))/ (x2(x2size)-x2(x2size-1)) )  * (xx2-x2(x2size))  + yyy(x2size,:)          
                     
                     
                           
                 
                 
             else 
                 
       
                             
    ff = ((yyy(2,:)-yyy(1,:))/ (x2(2)-x2(1)) )  * (xx2-x2(1))  + yyy(1,:)          
   
                       
  

             end if 
   
                 
                 
           two_linterp =  ((ff(x3size)-ff(x3size-1))/ (x3(x3size)-x3(x3size-1)) )  * (xx3-x3(x3size))  + ff(x3size)
           
             
        
             
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
                  
         else  ! else x3<x3(1)
         
         
         
                 if (xx2.LE.x2(x2size) .AND. xx2.GE.x2(1)) then 
                 
            call linint_Equi(xx2, x2(1), x2(x2size), x2size-1, il2, ir2, varphi2)
            
            il2=il2+1 
            ir2=ir2+1          
    
                 
 
   ff = yyy(il2,:)*varphi2 + yyy(ir2,:)*(1-varphi2)             
                     
                        
                 
                 
         elseif (xx2 > x2(x2size)) then            
                 
                 

    ff = ((yyy(x2size,:)-yyy(x2size-1,:))/ (x2(x2size)-x2(x2size-1)) )  * (xx2-x2(x2size))  + yyy(x2size,:)          
                     
                     
                           
                 
                 
             else 
                 
       
                             
    ff = ((yyy(2,:)-yyy(1,:))/ (x2(2)-x2(1)) )  * (xx2-x2(1))  + yyy(1,:)          
   
                       
  

             end if 
   
                 
                 
           two_linterp =  ((ff(2)-ff(1))/ (x3(2)-x3(1)) )  * (xx3-x3(1))  + ff(1)
           
           endif 

end function 

	
    end module
    
    
    
    
    !function two_interp(from, to, array)
!    REAL(8), INTENT(IN) :: from, to
!    REAL(8), DIMENSION(:), INTENT(OUT) :: array
!    REAL(8) :: range
!    integer :: nn, i
!    nn = size(array)
!    range = to - from
!
!    if (nn == 0) return
!
!    if (nn == 1) then
!        array(1) = from
!        return
!    end if
!
!
!    do i=1, nn
!        array(i) = from + range * (i - 1) / (nn - 1)
!    end do
!end function 

    
    
!    subroutine spline_coeff(x1,x2, x3, yy, coeff_x3)
!
!use toolbox
!
!    REAL(8), INTENT(IN) :: x1(:), x2(:), x3(:), yy(:,:,:)
!    real*8 :: coeff_x3(:,:,:)
!    !REAL(8) :: three_interp
!	REAL(8), dimension (:), allocatable :: yytemp, coeff 
!    REAL(8), dimension (:,:), allocatable :: yy_over3 
!    integer :: x1size, x2size, x3size
!	integer :: i, j, k 
!    
!	! interpolation over x1 x2 x3 
!	! spline interpolation over x3 
!	! and linear interpolation over x1 x2 	
!	
!	
!	x1size = SIZE(x1,1)
!	x2size = SIZE(x2,1)
!	x3size = SIZE(x3,1)
!    
!  allocate (coeff(x3size+2))	
!  !allocate (coeff_x3(x1size,x2size,x3size+2))
!  allocate (yytemp(x3size))
!  allocate (yy_over3(x1size,x2size))
!  
!	
!	
!	Do i=1,x1size 
!       
!       Do j = 1,x2size        
!           
!           Do k=1,x3size
!           
!           yytemp(k) = yy(i,j,k)
!
!           end Do 
!           
!           call spline_interp(yytemp, coeff)
!           
!           coeff_x3(i,j,:) = coeff 
!           
!         		 
!	   end Do 
!	   
!	end Do 
!	
!	
!         
!  end subroutine 
!    
!    
!    
!function three_interp(x1,x2, x3, coeff, yyy, xx1, xx2, xx3)
!
!use toolbox
!
!    REAL(8), INTENT(IN) :: x1(:), x2(:), x3(:), coeff(:,:,:), yyy(:,:,:), xx1, xx2, xx3 
!    REAL(8) :: three_interp, templ1, tempr1, templ2, tempr2, ttempl, ttempr, templ, tempr
!	REAL(8), dimension (:), allocatable :: coeff1, coeff2, coeff3, coeff4    
!	REAL(8) :: varphi2, varphi 
!	integer :: x1size, x2size, x3size, il2, ir2, il, ir 
!	integer :: i, j, k 
!    
!	! interpolation over x1 x2 x3 
!	! spline interpolation over x3 
!	! and linear interpolation over x1 x2 	
!	
!	
!	x1size = SIZE(x1,1)
!	x2size = SIZE(x2,1)
!	x3size = SIZE(x3,1)
!	
! allocate (coeff1(x3size+2))
! allocate (coeff2(x3size+2))
! allocate (coeff3(x3size+2))
! allocate (coeff4(x3size+2))
! 
!	          
!if (xx3 .LE. x3(x3size) .AND. xx3.GE.x3(1)) then   
!    
!             
!    if (xx2.LE.x2(x2size) .AND. xx2.GE.x2(1)) then 
!                 
!            call linint_Equi(xx2, x2(1), x2(x2size), x2size-1, il2, ir2, varphi2)
!            
!            il2=il2+1 
!            ir2=ir2+1
!            
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
!                 if (xx1 .LE. x1(x1size) .AND. xx1.ge.x1(1)) then 
!                     
!                        call linint_Equi(xx1, x1(1), x1(x1size), x1size-1, il, ir, varphi)                         
!                        
!                        
!                        il=il+1
!                        ir=ir+1
!                        
!                        coeff1 = coeff(il,il2,:)  
!                        
!                        coeff2 = coeff(ir,il2,:)
!                        
!                        coeff3 = coeff(il,ir2,:)
!                        
!                        coeff4 = coeff(ir,ir2,:)
!                        
!                        
!                        three_interp = varphi*varphi2*spline_eval(xx3, coeff1, x3(1), x3(x3size))  & 
!                            
!                             +  (1d0-varphi)*varphi2*spline_eval(xx3, coeff2, x3(1), x3(x3size))  & 
!                            
!                            +  varphi*(1d0-varphi2)*spline_eval(xx3, coeff3, x3(1), x3(x3size)) & 
!                            
!                            + (1d0-varphi)*(1d0-varphi2)*spline_eval(xx3, coeff4, x3(1), x3(x3size))                         
!                        
!                       
!            
!            
!                 elseif (xx1 > x1(x1size)) then 
!                     
!                     
!    coeff1 = coeff(x1size,il2,:)  
!    coeff2 = coeff(x1size-1,il2,:)
!    coeff3 = coeff(x1size,ir2,:)
!    coeff4 = coeff(x1size-1,ir2,:)
!                                         
!                     
!  templ =  (spline_eval(xx3, coeff1, x3(1), x3(x3size)) -  spline_eval(xx3, coeff2, x3(1), x3(x3size)))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + spline_eval(xx3, coeff1, x3(1), x3(x3size))
!  tempr =  (spline_eval(xx3, coeff3, x3(1), x3(x3size)) -  spline_eval(xx3, coeff4, x3(1), x3(x3size)))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + spline_eval(xx3, coeff3, x3(1), x3(x3size))
!  three_interp = varphi2*templ + (1d0-varphi2)*tempr 
!                 
!                 else 
! 
!    coeff1 = coeff(2,il2,:)  
!    coeff2 = coeff(1,il2,:)
!    coeff3 = coeff(2,ir2,:)
!    coeff4 = coeff(1,ir2,:)
!                                         
!                     
!  templ =  (spline_eval(xx3, coeff1, x3(1), x3(x3size)) -  spline_eval(xx3, coeff2, x3(1), x3(x3size)))/(x1(x1size)-x1(x1size-1))*(xx1-x1(1)) + spline_eval(xx3, coeff2, x3(1), x3(x3size))
!  tempr =  (spline_eval(xx3, coeff3, x3(1), x3(x3size)) -  spline_eval(xx3, coeff4, x3(1), x3(x3size)))/(x1(x1size)-x1(x1size-1))*(xx1-x1(1)) + spline_eval(xx3, coeff4, x3(1), x3(x3size))
!  three_interp = varphi2*templ + (1d0-varphi2)*tempr 
!         
!                     
!                     
!                 endif   
!                 
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11               
!                 
!                 
!         elseif (xx2 > x2(x2size)) then            
!                 
!                 
!if (xx1 .le. x1(x1size) .AND. xx1.ge.x1(1)) then 
!                            
!call linint_Equi(xx1, x1(1), x1(x1size), x1size-1, il, ir, varphi) 
!
!il=il+1
!ir=ir+1
!                     
!coeff1 = coeff(il,x2size,:)  
!coeff2 = coeff(il,x2size-1,:)
!coeff3 = coeff(ir,x2size,:)
!coeff4 = coeff(ir,x2size-1,:)                                       
!                     
!  templ =  (spline_eval(xx3, coeff1, x3(1), x3(x3size)) -  spline_eval(xx3, coeff2, x3(1), x3(x3size)))/(x2(x2size)-x2(x2size-1))*(xx2-x2(x2size)) + spline_eval(xx3, coeff1, x3(1), x3(x3size))
!  tempr =  (spline_eval(xx3, coeff3, x3(1), x3(x3size)) -  spline_eval(xx3, coeff4, x3(1), x3(x3size)))/(x2(x2size)-x2(x2size-1))*(xx2-x2(x2size)) + spline_eval(xx3, coeff3, x3(1), x3(x3size))
!  three_interp = varphi*templ + (1d0-varphi)*tempr                     
!
!
!                      elseif (xx1 > x1(x1size)) then 
!                     
!                     
!    coeff1 = coeff(x1size,x2size,:)  
!    coeff2 = coeff(x1size-1,x2size,:)
!    coeff3 = coeff(x1size,x2size-1,:)
!    coeff4 = coeff(x1size-1,x2size-1,:)
!                                         
!  ! value at x2=x2(x2size)                    
!  tempr =  (spline_eval(xx3, coeff1, x3(1), x3(x3size)) -  spline_eval(xx3, coeff2, x3(1), x3(x3size)))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + spline_eval(xx3, coeff1, x3(1), x3(x3size))
!  
!  ! value at x2=x2(x2size-1)                    
!  templ =  (spline_eval(xx3, coeff3, x3(1), x3(x3size)) -  spline_eval(xx3, coeff4, x3(1), x3(x3size)))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + spline_eval(xx3, coeff3, x3(1), x3(x3size))
!  
!  
!  three_interp =  (tempr -  templ)/(x2(x2size)-x2(x2size-1))*(xx2-x2(x2size)) + tempr
!  
!                 
!else 
! 
!    coeff1 = coeff(2,x2size,:)  
!    coeff2 = coeff(1,x2size,:)
!    coeff3 = coeff(2,x2size-1,:)
!    coeff4 = coeff(1,x2size-1,:)
!                                         
!  ! value at x2=x2(x2size)                    
!  tempr =  (spline_eval(xx3, coeff1, x3(1), x3(x3size)) -  spline_eval(xx3, coeff2, x3(1), x3(x3size)))/(x1(2)-x1(1))*(xx1-x1(1)) + spline_eval(xx3, coeff2, x3(1), x3(x3size))
!  
!  ! value at x2=x2(x2size-1)                    
!  templ =  (spline_eval(xx3, coeff3, x3(1), x3(x3size)) -  spline_eval(xx3, coeff4, x3(1), x3(x3size)))/(x1(2)-x1(1))*(xx1-x1(1)) + spline_eval(xx3, coeff4, x3(1), x3(x3size))
!  
!  
!  three_interp =  (tempr -  templ)/(x2(x2size)-x2(x2size-1))*(xx2-x2(x2size)) + tempr
!                     
!                     
!                 endif 
!                 
!                 
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11                     
!                 
!                 
!             else 
!                 
!       
!                             
!      if (xx1.le.x1(x1size) .AND. xx1.ge.x1(1)) then 
!                            
!call linint_Equi(xx1, x1(1), x1(x1size), x1size-1, il, ir, varphi) 
!
!il=il+1
!ir=ir+1
!                     
!coeff1 = coeff(il,2,:)  
!coeff2 = coeff(il,1,:)
!coeff3 = coeff(ir,2,:)
!coeff4 = coeff(ir,1,:)                                       
!                     
!  templ =  (spline_eval(xx3, coeff1, x3(1), x3(x3size)) -  spline_eval(xx3, coeff2, x3(1), x3(x3size)))/(x2(2)-x2(1))*(xx2-x2(1)) + spline_eval(xx3, coeff2, x3(1), x3(x3size))
!  tempr =  (spline_eval(xx3, coeff3, x3(1), x3(x3size)) -  spline_eval(xx3, coeff4, x3(1), x3(x3size)))/(x2(2)-x2(1))*(xx2-x2(1)) + spline_eval(xx3, coeff4, x3(1), x3(x3size))
!  three_interp = varphi2*templ + (1d0-varphi2)*tempr                     
!                
!
!
!        elseif (xx1 > x1(x1size)) then 
!                     
!                     
!    coeff1 = coeff(x1size,2,:)  
!    coeff2 = coeff(x1size-1,2,:)
!    coeff3 = coeff(x1size,1,:)
!    coeff4 = coeff(x1size-1,1,:)
!                                         
!  ! value at x2=x2(2)                    
!  tempr =  (spline_eval(xx3, coeff1, x3(1), x3(x3size)) -  spline_eval(xx3, coeff2, x3(1), x3(x3size)))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + spline_eval(xx3, coeff1, x3(1), x3(x3size))
!  
!  ! value at x2=x2(1)                    
!  templ =  (spline_eval(xx3, coeff3, x3(1), x3(x3size)) -  spline_eval(xx3, coeff4, x3(1), x3(x3size)))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + spline_eval(xx3, coeff3, x3(1), x3(x3size))
!  
!  
!  three_interp =  (tempr -  templ)/(x2(2)-x2(1))*(xx2-x2(1)) + templ
!  
!                 
!else 
! 
!    coeff1 = coeff(2,2,:)  
!    coeff2 = coeff(1,2,:)
!    coeff3 = coeff(2,1,:)
!    coeff4 = coeff(1,1,:)
!                                         
!  ! value at x2=x2(2)                    
!  tempr =  (spline_eval(xx3, coeff1, x3(1), x3(x3size)) -  spline_eval(xx3, coeff2, x3(1), x3(x3size)))/(x1(2)-x1(1))*(xx1-x1(1)) + spline_eval(xx3, coeff2, x3(1), x3(x3size))
!  
!  ! value at x2=x2(1)                    
!  templ =  (spline_eval(xx3, coeff3, x3(1), x3(x3size)) -  spline_eval(xx3, coeff4, x3(1), x3(x3size)))/(x1(2)-x1(1))*(xx1-x1(1)) + spline_eval(xx3, coeff4, x3(1), x3(x3size))
!  
!  
!   three_interp =  (tempr -  templ)/(x2(2)-x2(1))*(xx2-x2(1)) + templ
!  
!                     
!                     
!endif    
!
!
!             end if 
!             
!                 
!                     
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!             
!    elseif (xx3 > x3(x3size)) then 
!
!                 
!    if (xx2.le.x2(x2size) .AND. xx2.ge.x2(1)) then 
!                 
!            call linint_Equi(xx2, x2(1), x2(x2size), x2size-1, il2, ir2, varphi2)
!            
!            il2=il2+1
!            ir2=ir2+1
!            
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
!                       if (xx1.le.x1(x1size) .AND. xx1.ge.x1(1)) then 
!                     
!  call linint_Equi(xx1, x1(1), x1(x1size), x1size-1, il, ir, varphi) 
!  
!  il=il+1
!  ir=ir+1
!
!tempr1 = yyy(il,ir2,x3size)*varphi + yyy(ir,ir2,x3size)*(1d0-varphi) ! x1,x2=x2size, x3=x3size
!templ1 = yyy(il,il2,x3size)*varphi + yyy(ir,il2,x3size)*(1d0-varphi) ! x1,x2=x2size-1, x3=x3size
!tempr2 = yyy(il,ir2,x3size-1)*varphi + yyy(ir,ir2,x3size-1)*(1d0-varphi) ! x1,x2=x2size, x3=x3size-1
!templ2= yyy(il,il2,x3size-1)*varphi + yyy(ir,il2,x3size-1)*(1d0-varphi) ! x1,x2=x2size-1, x3=x3size-1
!
!
!    elseif (xx1 > x1(x1size)) then 
!    
!
!tempr1 = (yyy(x1size,ir2,x3size)-yyy(x1size-1,ir2,x3size))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,ir2,x3size)  ! x1,x2=x2size, x3=x3size
!templ1 = (yyy(x1size,il2,x3size)-yyy(x1size-1,il2,x3size))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,il2,x3size) ! x1,x2=x2size-1, x3=x3size
!tempr2 = (yyy(x1size,ir2,x3size-1)-yyy(x1size-1,ir2,x3size-1))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,ir2,x3size-1) ! x1,x2=x2size, x3=x3size-1
!templ2= (yyy(x1size,il2,x3size-1)-yyy(x1size-1,il2,x3size-1))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,il2,x3size-1) ! x1,x2=x2size-1, x3=x3size-1
!
!    
!    
!    
!    else 
! 
!tempr1 = (yyy(2,ir2,x3size)-yyy(1,ir2,x3size))/ (x1(2)-x1(1))*(xx1-x1(1)) + yyy(1,ir2,x3size)  ! x1,x2=x2size, x3=x3size
!templ1 = (yyy(2,il2,x3size)-yyy(1,il2,x3size))/(x1(2)-x1(1))*(xx1-x1(1)) + yyy(1,il2,x3size) ! x1,x2=x2size-1, x3=x3size
!tempr2 = (yyy(2,ir2,x3size-1)-yyy(1,ir2,x3size-1))/(x1(2)-x1(1))*(xx1-x1(1))+ yyy(1,ir2,x3size-1) ! x1,x2=x2size, x3=x3size-1
!templ2= (yyy(2,il2,x3size-1)-yyy(1,il2,x3size-1))/(x1(2)-x1(1))*(xx1-x1(1)) + yyy(1,il2,x3size-1) ! x1,x2=x2size-1, x3=x3size-1
!
!  
!    
!                    
!                     
!    endif   
!    
!    ttempr = templ1*varphi2 + tempr1*(1d0-varphi2)
!
!    ttempl = templ2*varphi2 + tempr2*(1d0-varphi2)
!                 
!                 
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11               
!                 
!                 
!    elseif (xx2 > x2(x2size)) then                  
!                 
!        if (xx1.le.x1(x1size) .AND. xx1.ge.x1(1)) then 
!                     
!  call linint_Equi(xx1, x1(1), x1(x1size), x1size-1, il, ir, varphi) 
!  
!  il=il+1
!  ir=ir+1
!
!tempr1 = yyy(il,x2size,x3size)*varphi + yyy(ir,x2size,x3size)*(1d0-varphi) ! x1,x2=x2size, x3=x3size
!templ1 = yyy(il,x2size-1,x3size)*varphi +yyy(ir,x2size-1,x3size)*(1d0-varphi) ! x1,x2=x2size-1, x3=x3size
!tempr2 = yyy(il,x2size,x3size-1)*varphi + yyy(ir,x2size,x3size-1)*(1d0-varphi) ! x1,x2=x2size, x3=x3size-1
!templ2= yyy(il,x2size-1,x3size-1)*varphi + yyy(ir,x2size-1,x3size-1)*(1d0-varphi) ! x1,x2=x2size-1, x3=x3size-1
!
!
!    elseif (xx1 > x1(x1size)) then 
!    
!
!tempr1 = (yyy(x1size,x2size,x3size)-yyy(x1size-1,x2size,x3size))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,x2size,x3size)  ! x1,x2=x2size, x3=x3size
!templ1 = (yyy(x1size,x2size-1,x3size)-yyy(x1size-1,x2size-1,x3size))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,x2size-1,x3size) ! x1,x2=x2size-1, x3=x3size
!tempr2 = (yyy(x1size,x2size,x3size-1)-yyy(x1size-1,x2size,x3size-1))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,x2size,x3size-1) ! x1,x2=x2size, x3=x3size-1
!templ2= (yyy(x1size,x2size-1,x3size-1)-yyy(x1size-1,x2size-1,x3size-1))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,x2size-1,x3size-1) ! x1,x2=x2size-1, x3=x3size-1
!
!    
!    
!    
!    else 
! 
!tempr1 = (yyy(2,x2size,x3size)-yyy(1,x2size,x3size))/(x1(2)-x1(1))*(xx1-x1(1)) + yyy(1,x2size,x3size)  ! x1,x2=x2size, x3=x3size
!templ1 = (yyy(2,x2size-1,x3size)-yyy(1,x2size-1,x3size))/(x1(2)-x1(1))*(xx1-x1(1)) + yyy(1,x2size-1,x3size) ! x1,x2=x2size-1, x3=x3size
!tempr2 = (yyy(2,x2size,x3size-1)-yyy(1,x2size,x3size-1))/(x1(2)-x1(1))*(xx1-x1(1))+ yyy(1,x2size,x3size-1) ! x1,x2=x2size, x3=x3size-1
!templ2= (yyy(2,x2size-1,x3size-1)-yyy(1,x2size-1,x3size-1))/(x1(2)-x1(1))*(xx1-x1(1)) + yyy(1,x2size-1,x3size-1) ! x1,x2=x2size-1, x3=x3size-1
!
!  
!    
!                    
!                     
!    endif   
!    
!    ttempr = ((tempr1-templ1)/ (x2(x2size)-x2(x2size-1)) )  * (xx2-x2(x2size))  + tempr1 
!
!    ttempl = ((tempr2-templ2)/ (x2(x2size)-x2(x2size-1)) )  * (xx2-x2(x2size))  + tempr2
!                 
!                 
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11                     
!                 
!                 
!    else 
!        ! (xx2 .LE. x2(1))                
!       
!                             
!       if (xx1.le.x1(x1size) .AND. xx1.ge.x1(1)) then 
!                     
!  call linint_Equi(xx1, x1(1), x1(x1size), x1size-1, il, ir, varphi) 
!  
!  il=il+1
!  ir=ir+1
!
!tempr1 = yyy(il,2,x3size)*varphi + yyy(ir,2,x3size)*(1d0-varphi) ! x1,x2=2, x3=x3size
!templ1 = yyy(il,1,x3size)*varphi + yyy(ir,1,x3size)*(1d0-varphi) ! x1,x2=1, x3=x3size
!tempr2 = yyy(il,2,x3size-1)*varphi + yyy(ir,2,x3size-1)*(1d0-varphi) ! x1,x2=2, x3=x3size-1
!templ2= yyy(il,1,x3size-1)*varphi + yyy(ir,1,x3size-1)*(1d0-varphi) ! x1,x2=1, x3=x3size-1
!
!
!    elseif (xx1 > x1(x1size)) then 
!    
!
!tempr1 = (yyy(x1size,2,x3size)-yyy(x1size-1,2,x3size))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,2,x3size)  ! x1,x2=2, x3=x3size
!templ1 = (yyy(x1size,1,x3size)-yyy(x1size-1,1,x3size))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,1,x3size) ! x1,x2=1, x3=x3size
!tempr2 = (yyy(x1size,2,x3size-1)-yyy(x1size-1,2,x3size-1))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,2,x3size-1) ! x1,x2=2, x3=x3size-1
!templ2= (yyy(x1size,1,x3size-1)-yyy(x1size-1,1,x3size-1))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,1,x3size-1) ! x1,x2=1, x3=x3size-1
!
!    
!    
!    
!    else 
! 
!tempr1 = (yyy(2,2,x3size)-yyy(1,2,x3size))/  (x1(2)-x1(1))  *  (xx1-x1(1)) + yyy(1,2,x3size)  ! x1,x2=2, x3=x3size
!templ1 = (yyy(2,1,x3size)-yyy(1,1,x3size))/  (x1(2)-x1(1))  *  (xx1-x1(1))  + yyy(1,1,x3size) ! x1,x2=1, x3=x3size
!tempr2 = (yyy(2,2,x3size-1)-yyy(1,2,x3size-1))/  (x1(2)-x1(1))  *  (xx1-x1(1))  + yyy(1,2,x3size-1) ! x1,x2=2, x3=x3size-1
!templ2= (yyy(2,1,x3size-1)-yyy(1,1,x3size-1))/  (x1(2)-x1(1))  *  (xx1-x1(1)) + yyy(1,1,x3size-1) ! x1,x2=1, x3=x3size-1
!    
!                    
!                     
!    endif   
!    
!    ttempr = ((tempr1-templ1)/ (x2(2)-x2(1)) )  * (xx2-x2(1))  + templ1 
!
!    ttempl = ((tempr2-templ2)/ (x2(2)-x2(1)) )  * (xx2-x2(1))  + templ2
!       
!    
!     endif                     
!   
!                 
!                 
!           three_interp =  ((ttempr-ttempl)/ (x3(x3size)-x3(x3size-1)) )  * (xx3-x3(x3size))  + ttempr 
!         
!             
!        
!             
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!                  
!         else 
!         
!         
!         
!             if (xx2.le.x2(x2size) .AND. xx2.ge.x2(1)) then 
!                 
!            call linint_Equi(xx2, x2(1), x2(x2size), x2size-1, il2, ir2, varphi2)
!            
!            il2=il2+1
!            ir2=ir2+1
!            
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
!                       if (xx1.le.x1(x1size) .AND. xx1.ge.x1(1)) then 
!                     
!  call linint_Equi(xx1, x1(1), x1(x1size), x1size-1, il, ir, varphi) 
!  
!  il=il+1
!  ir=ir+1
!
!tempr1 = yyy(il,ir2,2)*varphi + yyy(ir,ir2,2)*(1d0-varphi) ! x1,x2=x2size, x3=x3size
!templ1 = yyy(il,il2,2)*varphi + yyy(ir,il2,2)*(1d0-varphi) ! x1,x2=x2size-1, x3=x3size
!tempr2 = yyy(il,ir2,1)*varphi + yyy(ir,ir2,1)*(1d0-varphi) ! x1,x2=x2size, x3=x3size-1
!templ2= yyy(il,il2,1)*varphi + yyy(ir,il2,1)*(1d0-varphi) ! x1,x2=x2size-1, x3=x3size-1
!
!
!    elseif (xx1 > x1(x1size)) then 
!    
!
!tempr1 = (yyy(x1size,ir2,2)-yyy(x1size-1,ir2,2))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,ir2,2)  ! x1,x2=x2size, x3=x3size
!templ1 = (yyy(x1size,il2,2)-yyy(x1size-1,il2,2))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,il2,2) ! x1,x2=x2size-1, x3=x3size
!tempr2 = (yyy(x1size,ir2,1)-yyy(x1size-1,ir2,1))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,ir2,1) ! x1,x2=x2size, x3=x3size-1
!templ2= (yyy(x1size,il2,1)-yyy(x1size-1,il2,1))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,il2,1) ! x1,x2=x2size-1, x3=x3size-1
!
!    
!    
!    
!    else 
! 
!tempr1 = (yyy(2,ir2,2)-yyy(1,ir2,2))/(x1(2)-x1(1))*(xx1-x1(1)) + yyy(1,ir2,2)  ! x1,x2=x2size, x3=x3size
!templ1 = (yyy(2,il2,2)-yyy(1,il2,2))/(x1(2)-x1(1))*(xx1-x1(1)) + yyy(1,il2,2) ! x1,x2=x2size-1, x3=x3size
!tempr2 = (yyy(2,ir2,1)-yyy(1,ir2,1))/(x1(2)-x1(1))*(xx1-x1(1))+ yyy(1,ir2,1) ! x1,x2=x2size, x3=x3size-1
!templ2= (yyy(2,il2,1)-yyy(1,il2,1))/(x1(2)-x1(1))*(xx1-x1(1)) + yyy(1,il2,1) ! x1,x2=x2size-1, x3=x3size-1
!
!  
!    
!                    
!                     
!    endif   
!    
!    ttempr = templ1*varphi2 + tempr1*(1d0-varphi2)
!
!    ttempl = templ2*varphi2 + tempr2*(1d0-varphi2)
!                 
!                 
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11               
!                 
!                 
!    elseif (xx2 > x2(x2size)) then                  
!                 
!        if (xx1.le.x1(x1size) .AND. xx1.ge.x1(1)) then 
!                     
!  call linint_Equi(xx1, x1(1), x1(x1size), x1size-1, il, ir, varphi) 
!  
!  il=il+1
!  ir=ir+1
!
!tempr1 = yyy(il,x2size,2)*varphi + yyy(ir,x2size,2)*(1d0-varphi) ! x1,x2=x2size, x3=x3size
!templ1 = yyy(il,x2size-1,2)*varphi + yyy(ir,x2size-1,2)*(1d0-varphi) ! x1,x2=x2size-1, x3=x3size
!tempr2 = yyy(il,x2size,1)*varphi + yyy(ir,x2size,1)*(1d0-varphi) ! x1,x2=x2size, x3=x3size-1
!templ2= yyy(il,x2size-1,1)*varphi + yyy(ir,x2size-1,1)*(1d0-varphi) ! x1,x2=x2size-1, x3=x3size-1
!
!
!    elseif (xx1 > x1(x1size)) then 
!    
!
!tempr1 = (yyy(x1size,x2size,2)-yyy(x1size-1,x2size,2))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,x2size,2)  ! x1,x2=x2size, x3=x3size
!templ1 = (yyy(x1size,x2size-1,2)-yyy(x1size-1,x2size-1,2))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,x2size-1,2) ! x1,x2=x2size-1, x3=x3size
!tempr2 = (yyy(x1size,x2size,1)-yyy(x1size-1,x2size,1))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,x2size,1) ! x1,x2=x2size, x3=x3size-1
!templ2= (yyy(x1size,x2size-1,1)-yyy(x1size-1,x2size-1,1))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,x2size-1,1) ! x1,x2=x2size-1, x3=x3size-1
!
!    
!    
!    
!    else 
! 
!tempr1 = (yyy(2,x2size,2)-yyy(1,x2size,2))/(x1(2)-x1(1))*(xx1-x1(1)) + yyy(1,x2size,2)  ! x1,x2=x2size, x3=x3size
!templ1 = (yyy(2,x2size-1,2)-yyy(1,x2size-1,2))/(x1(2)-x1(1))*(xx1-x1(1)) + yyy(1,x2size-1,2) ! x1,x2=x2size-1, x3=x3size
!tempr2 = (yyy(2,x2size,1)-yyy(1,x2size,1))/(x1(2)-x1(1))*(xx1-x1(1))+ yyy(1,x2size,1) ! x1,x2=x2size, x3=x3size-1
!templ2= (yyy(2,x2size-1,1)-yyy(1,x2size-1,1))/(x1(2)-x1(1))*(xx1-x1(1)) + yyy(1,x2size-1,1) ! x1,x2=x2size-1, x3=x3size-1
!
!  
!    
!                    
!                     
!    endif   
!    
!    ttempr = ((tempr1-templ1)/ (x2(x2size)-x2(x2size-1)) )  * (xx2-x2(x2size))  + tempr1 
!
!    ttempl = ((tempr2-templ2)/ (x2(x2size)-x2(x2size-1)) )  * (xx2-x2(x2size))  + tempr2
!                 
!                 
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11                     
!                 
!                 
!    else 
!        ! (xx2 .LE. x2(1))                
!       
!                             
!       if (xx1.le.x1(x1size) .AND. xx1.ge.x1(1)) then 
!                     
!  call linint_Equi(xx1, x1(1), x1(x1size), x1size-1, il, ir, varphi) 
!  
!  il=il+1
!  ir=ir+1
!
!tempr1 = yyy(il,2,2)*varphi + yyy(ir,2,2)*(1d0-varphi) ! x1,x2=2, x3=x3size
!templ1 = yyy(il,1,2)*varphi + yyy(ir,1,2)*(1d0-varphi) ! x1,x2=1, x3=x3size
!tempr2 = yyy(il,2,1)*varphi + yyy(ir,2,1)*(1d0-varphi) ! x1,x2=2, x3=x3size-1
!templ2= yyy(il,1,1)*varphi + yyy(ir,1,1)*(1d0-varphi) ! x1,x2=1, x3=x3size-1
!
!
!    elseif (xx1 > x1(x1size)) then 
!    
!
!tempr1 = (yyy(x1size,2,2)-yyy(x1size-1,2,2))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,2,2)  ! x1,x2=2, x3=x3size
!templ1 = (yyy(x1size,1,2)-yyy(x1size-1,1,2))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,1,2) ! x1,x2=1, x3=x3size
!tempr2 = (yyy(x1size,2,1)-yyy(x1size-1,2,1))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,2,1) ! x1,x2=2, x3=x3size-1
!templ2= (yyy(x1size,1,1)-yyy(x1size-1,1,1))/(x1(x1size)-x1(x1size-1))*(xx1-x1(x1size)) + yyy(x1size,1,1) ! x1,x2=1, x3=x3size-1
!
!    
!    
!    
!    else 
! 
!tempr1 = (yyy(2,2,2)-yyy(1,2,2))/  (x1(2)-x1(1))  *  (xx1-x1(1)) + yyy(1,2,2)  ! x1,x2=2, x3=x3size
!templ1 = (yyy(2,1,2)-yyy(1,1,2))/  (x1(2)-x1(1))  *  (xx1-x1(1))  + yyy(1,1,2) ! x1,x2=1, x3=x3size
!tempr2 = (yyy(2,2,1)-yyy(1,2,1))/  (x1(2)-x1(1))  *  (xx1-x1(1))  + yyy(1,2,1) ! x1,x2=2, x3=x3size-1
!templ2= (yyy(2,1,1)-yyy(1,1,1))/  (x1(2)-x1(1))  *  (xx1-x1(1)) + yyy(1,1,1) ! x1,x2=1, x3=x3size-1
!    
!                    
!                     
!    endif   
!    
!    ttempr = ((tempr1-templ1)/ (x2(2)-x2(1)) )  * (xx2-x2(1))  + templ1 
!
!    ttempl = ((tempr2-templ2)/ (x2(2)-x2(1)) )  * (xx2-x2(1))  + templ2
!       
!    
!     endif                     
!   
!                 
!                 
!           three_interp =  ((ttempr-ttempl)/ (x3(2)-x3(1)) )  * (xx3-x3(1))  + ttempl 
!             
!         
!             
!         endif   
!
!end function    