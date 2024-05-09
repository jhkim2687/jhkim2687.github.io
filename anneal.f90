MODULE ANNEALING
  USE nrutil, ONLY : assert_eq,imaxloc,iminloc,swap
  IMPLICIT NONE
  PRIVATE 
  PUBLIC :: Simulated_Annealing

  INTEGER(4), SAVE :: kiss1=123456789,kiss2=362436069,kiss3=521288629,kiss4=380116160,seed_anneal,zigjz_anneal
  LOGICAL, SAVE :: initialized=.FALSE.
CONTAINS

  SUBROUTINE Simulated_Annealing(temp0,tfactor,nsteps,x,func,max_jiter,id,k1)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: id
    REAL(8), INTENT(IN) :: temp0,tfactor
    INTEGER, INTENT(IN) :: nsteps
    INTEGER, INTENT(IN) :: max_jiter
    INTEGER, INTENT(IN), OPTIONAL :: k1(4)
    REAL(8), INTENT(INOUT) :: x(:)
    REAL(8), PARAMETER :: FTOL = 1.0d-8
    INTERFACE
       REAL(8) FUNCTION FUNC(x)
         IMPLICIT NONE
         REAL(8), INTENT(IN) :: x(:)
       END FUNCTION FUNC
    END INTERFACE
    LOGICAL :: write 
    INTEGER :: i,j,iter,jiter,nit,ndim,ndimp1,tin,tout,hz
    REAL(8) :: yb,ybb,p(SIZE(x,1)+1,SIZE(x,1)),y(SIZE(x,1)+1),pb(SIZE(x,1))
    REAL(8) :: temp,elapsed_time,ver
    IF (.NOT.INITIALIZED) THEN
       IF (PRESENT(k1)) THEN
          CALL SET_SEED_ANNEAL(k1(1),k1(2),k1(3),k1(4))
       ELSE
          CALL SET_SEED_ANNEAL()
       END IF
       INITIALIZED=.TRUE.
    END IF
    write=.FALSE.
    ndim=SIZE(x,1)
    ndimp1=ndim+1
    ! initialize the simplex, seting P[1][:]=x, P[i][:]=x+ei, where {ei} is the standard basis for Rn
    DO i=2,ndimp1
       DO j=1,ndim
          IF (i-1==j) THEN
             ver=ABS(x(j))*0.5d0
             IF (ver<1.0d-10) ver=1.0d0
             pb(j) = x(j) + ver
             p(i,j) = x(j) + ver
          ELSE
             pb(j) = x(j)
             p(i,j) = x(j)                
          END IF
       END DO
       y(i)=FUNC(pb)
    END DO
    DO j=1,ndim
       pb(j)=x(j)
       p(1,j)=x(j)
    END DO
    y(1)=FUNC(pb)
    ybb=y(1)
    yb=y(1)
    
    nit=0
    temp=temp0
    DO jiter=1,max_jiter
       CALL SYSTEM_CLOCK(count_rate=hz) 
       CALL SYSTEM_CLOCK(count=tin) 
       iter=nsteps
       temp = temp*tfactor
       CALL amebsa(p,y,pb,yb,FTOL,func,iter,temp)
       nit = nit + nsteps - iter
       IF (yb< ybb) THEN
          ybb=yb
          IF (id==0) THEN
             OPEN(6424,file="point.txt",status="replace")
             DO j = 1, ndim
                WRITE(6424,'(F32.16)') pb(j)
             END DO
             CLOSE(6424)	
          END IF
       END IF
       CALL SYSTEM_CLOCK(count=tout) 
       tout=tout-tin
       elapsed_time=DBLE(tout)/DBLE(hz)
       IF (id==0.0d0) THEN
          write(*,'(A,I10,A,F16.6,A,F16.6,A,F16.6)') 'iterations:',jiter,' temp:',temp, &
            ' value:',yb,' time:',elapsed_time
          CALL FLUSH(6)
       END IF
       IF (iter>0) THEN
          IF (id==0) WRITE(*,*) 'FTOL Convergence'
          x=pb
          RETURN
       END IF
    END DO
    
    x=pb
    IF (id==0) THEN
       IF (jiter>=MAX_JITER) WRITE(*,'(A,I8,A)') 'No FTOL convergence, lowered temperature',jiter-1,'times'
       WRITE(*,*) 'Lowest y values',yb
       WRITE(*,*) 'Normal Completion'
    END IF
    RETURN
  END SUBROUTINE Simulated_Annealing
  
  SUBROUTINE amebsa(p,y,pb,yb,ftol,func,iter,temptr)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: iter
    REAL(8), INTENT(INOUT) :: yb
    REAL(8), INTENT(IN) :: ftol,temptr
    REAL(8), DIMENSION(:), INTENT(INOUT) :: y,pb
    REAL(8), DIMENSION(:,:), INTENT(INOUT) :: p
    INTERFACE
       REAL(8) FUNCTION FUNC(x)
         IMPLICIT NONE
         REAL(8), INTENT(IN) :: x(:)
       END FUNCTION FUNC
    END INTERFACE
    INTEGER, PARAMETER :: NMAX=200
    INTEGER :: ihi,ndim
    REAL(8) :: yhi
    REAL(8), DIMENSION(size(p,2)) :: psum
    CALL amebsa_private
  CONTAINS
    
    SUBROUTINE amebsa_private
      INTEGER :: i,ilo,inhi
      REAL(8) :: rtol,ylo,ynhi,ysave,ytry
      REAL(8), DIMENSION(size(y)) :: yt,harvest
      ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,size(pb),'amebsa')
      psum(:)=sum(p(:,:),dim=1)
      DO
         harvest=Sample_Uniform_Anneal()
         yt(:)=y(:)-temptr*log(harvest)
         ilo=iminloc(yt(:))
         ylo=yt(ilo)
         ihi=imaxloc(yt(:))
         yhi=yt(ihi)
         yt(ihi)=ylo
         inhi=imaxloc(yt(:))
         ynhi=yt(inhi)
         rtol=2.0d0*abs(yhi-ylo)/(abs(yhi)+abs(ylo))
         if (rtol < ftol .or. iter < 0) then
            call swap(y(1),y(ilo))
            call swap(p(1,:),p(ilo,:))
            RETURN
         end if
         ytry=amotsa(-1.0d0)
         iter=iter-1
         if (ytry <= ylo) then
            ytry=amotsa(2.0d0)
            iter=iter-1
         else if (ytry >= ynhi) then
            ysave=yhi
            ytry=amotsa(0.5d0)
            iter=iter-1
            if (ytry >= ysave) then
               p(:,:)=0.5d0*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
               do i=1,ndim+1
                  if (i /= ilo) y(i)=func(p(i,:))
               end do
               iter=iter-ndim
               psum(:)=sum(p(:,:),dim=1)
            end if
         end if
      end do
    END SUBROUTINE amebsa_private
    
    FUNCTION amotsa(fac)
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: fac
      REAL(8) :: amotsa
      REAL(8) :: fac1,fac2,yflu,ytry,harv
      REAL(8), DIMENSION(size(p,2)) :: ptry
      fac1=(1.0d0-fac)/ndim
      fac2=fac1-fac
      ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
      ytry=func(ptry)
      if (ytry <= yb) then
         pb(:)=ptry(:)
         yb=ytry
      end if
      harv=Sample_Uniform_Anneal()
      yflu=ytry+temptr*log(harv)
      if (yflu < yhi) then
         y(ihi)=ytry
         yhi=yflu
         psum(:)=psum(:)-p(ihi,:)+ptry(:)
         p(ihi,:)=ptry(:)
      end if
      amotsa=yflu
    END FUNCTION amotsa
  END SUBROUTINE amebsa
  
  SUBROUTINE Set_Seed_Anneal(k1,k2,k3,k4)
    IMPLICIT NONE
    INTEGER(4), INTENT(IN), OPTIONAL :: k1,k2,k3,k4
    INTEGER :: i
    IF (PRESENT(k1)) THEN
       kiss1=k1
       kiss2=k2
       kiss3=k3
       kiss4=k4
    ELSE
       CALL SYSTEM_CLOCK(seed_anneal)
       kiss1=SHR3_ANNEAL()
       kiss2=SHR3_ANNEAL()
       kiss3=SHR3_ANNEAL()
       kiss4=SHR3_ANNEAL()
    END IF
  END SUBROUTINE Set_Seed_Anneal

  ! Generate random 32-bit integers
  INTEGER FUNCTION shr3_Anneal( )
    zigjz_anneal = seed_anneal
    seed_anneal = IEOR(seed_anneal,ISHFT(seed_anneal,13))
    seed_anneal = IEOR(seed_anneal,ISHFT(seed_anneal,-17))
    seed_anneal = IEOR(seed_anneal,ISHFT(seed_anneal,5))
    shr3_anneal = zigjz_anneal + seed_anneal
  END FUNCTION shr3_anneal
  
  INTEGER(4) FUNCTION KISS_ANNEAL ()
    ! The  KISS (Keep It Simple Stupid) random number generator. Combines:
    ! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
    ! (2) A 3-shift shift-register generator, period 2^32-1,
    ! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
    !  Overall period>2^123;
    kiss1 = 69069 * kiss1 + 1327217885
    kiss2 = m (m (m (kiss2, 13), - 17), 5)
    kiss3 = 18000 * iand (kiss3, 65535) + ishft (kiss3, - 16)
    kiss4 = 30903 * iand (kiss4, 65535) + ishft (kiss4, - 16)
    kiss_anneal = kiss1 + kiss2 + ishft (kiss3, 16) + kiss4
  CONTAINS
    INTEGER FUNCTION m(k, n)
      INTEGER :: k, n
      m = ieor (k, ishft (k, n) )
    END FUNCTION m
  END FUNCTION KISS_ANNEAL

  REAL(8) FUNCTION Sample_Uniform_Anneal()
    Sample_Uniform_anneal = (DBLE(KISS_ANNEAL())+2147483649.0d0)/4294967297.0d0
  END FUNCTION Sample_Uniform_Anneal

END MODULE ANNEALING

