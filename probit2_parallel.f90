
program probit2_parallel

use GLOBALS
use DATAIN
use Objective
use toolbox
USE MINIMIZATION
use MPI 

implicit none

!!!! declare variables 

integer :: i,j
REAL(8) :: theta(N_COV)


CALL MPI_INIT(ierr ) ! Initialize the MPI execution environment
CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr ) !  Determines the rank of the calling process in the communicator, myid : rank of calling process 
CALL MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr ) ! Determines the size of the group associated with a communicator

write(*,fmt='(a,i4,a,i4)') 'I am  ', myid, ' of ', numprocs   
CALL READDATA(myid)


IF (myid==0) THEN ! process with myid==0 reads initial values   
OPEN(6424,file='initial.txt')
  DO i=1,N_COV
     READ(6424,fmt=*) theta(i)
  END DO
  CLOSE(6424) 
  
  OPEN(6440,file='betahat.txt')
  DO i=1,N_COV
     write(6440,fmt=*) betahat(i)
  END DO
  CLOSE(6440) 
  
WRITE(*,fmt='(f20.2)') theta
WRITE(*,fmt='(f20.2)') betahat
  
End if 

If (myid==0) then 
    write(*,fmt='(5f20.2)') (u(1,i),i=N_obs-3,N_OBS)
end if 
    
CALL MPI_BCAST(theta,N_COV,MPI_REAL8,0,MPI_COMM_WORLD,ierr)   

!*****************************************************************************************************
! broadcast a message from the process with rank root to all processes of the group, itself included *
! MPI_Bcast(buffer, count, datatype, root, comm)                                                     * 
! INOUT buffer starting address of buffer (choice)                                                   *
! IN count number of entries in buffer                                                               *
! IN datatype data type of buffer (handle)                                                           *
! IN root rank of broadcast root (integer)                                                           *
! IN comm communicator (handle)   
! Summary :passing message to processor with theta (input value)   
!******************************************************8**********************************************

CALL BFGS(theta,1.0d-6,1.0d-6,obj,myid,.TRUE.)

IF (myid==0) THEN 
WRITE(*,fmt='(f20.2)') theta
call plot(theta,betahat) 
call execplot( )
end if 

CALL MPI_FINALIZE(ierr)

end program 
