MODULE countline

IMPLICIT NONE 

CONTAINS 

Function nlines(filename)
  implicit none 
  character(len=*),intent(in) :: filename
  integer :: io, nlines 
  
  open(10,file=filename, iostat=io, status='old')
  if (io/=0) stop 'Cannot open file! '

  nlines = 0
  
  Do 
  read(10,*,iostat=io)
  if (io/=0) exit
  nlines = nlines + 1
  
  end do
  
  close(10)

end function nlines

end MODULE 

Module GLOBALS 

IMPLICIT NONE 

INTEGER, PARAMETER :: N_VAR=4,N_COV=3,S=5
REAL(8) :: beta(N_COV),betahat(N_COV)
INTEGER :: myid,numprocs, N_OBS,ierr ! myid: determine worker ! numprocs : number of processors
REAL, dimension (:,:), allocatable :: x,dataset,u 
REAL, dimension(:), allocatable :: y

contains 

SUBROUTINE SETSIZE 

USE COUNTLINE
implicit none 

N_OBS=nlines('data_ps1.raw')

allocate (x(N_OBS,N_COV))
allocate(y(N_OBS))

allocate (dataset(N_OBS,N_VAR+1))
allocate (u(S,N_OBS))

END SUBROUTINE 

end Module 

 

MODULE DATAIN 
USE GLOBALS
USE toolbox, only : simulate_normal, lu_invert
use MPI
IMPLICIT NONE 

CONTAINS

subroutine READDATA(id) 

implicit none 
INTEGER, INTENT(IN) :: id
INTEGER :: i,j
REAL(8) :: mu,var,XX(N_COV,N_COV),XY(N_COV),utemp(S) 
CALL SETSIZE

mu=0.0d0 
var=1.0d0

IF (id==0) open(10, file = 'data_ps1.raw')

Do i=1,N_OBS

IF (id==0) then 
    Read(10,*) dataset(i,1:N_VAR)
    dataset(i,N_VAR+1)=1.0d0
end if 

CALL MPI_BCAST(dataset(i,:),N_VAR+1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

end Do 

IF (id==0) close(10)


y=dataset(:,4)
x(:,2)=dataset(:,2)
x(:,3)=dataset(:,3)
x(:,1)=dataset(:,5) 


!! instantly message passing. 

Do i=1,N_OBS

IF (id==0) then 
    call simulate_normal(utemp,mu,var) 
    u(:,i)= utemp 
end if

CALL MPI_BCAST(u(:,i),S,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

end Do



!Do i=1,N_OBS

!IF (id==0) then 
 !   call simulate_normal(utemp,mu,var) 
!end if
 !   u(:,i)= utemp 
!end Do


!Do j=1,S
!CALL MPI_BCAST(u(j,:),N_OBS,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!end Do 

XX=matmul(transpose(X),X)  
XY=matmul(transpose(X),y) 
betahat=matmul(lu_invert(XX),XY) 

IF (id==0) WRITE(*,*) "dataset read"

end subroutine 

END 
