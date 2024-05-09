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
