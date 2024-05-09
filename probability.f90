MODULE PROBABILITY
  USE NRUTIL, ONLY : NRERROR,NR_GAMMLN,NR_GAMMINC,NR_SMALL,NR_BIG
  IMPLICIT NONE  
CONTAINS
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PDF's !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(8) FUNCTION PDF_Normal(z,mu,varin)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: z,mu,varin
    REAL(8), PARAMETER :: logroot2pi = 0.918938533204672780563271317078D0
    REAL(8) :: arg,std,var
    var=varin
    IF (var<0.0d0) CALL NRERROR('variance has to be positive: PDF_Normal')
    IF (var<NR_SMALL) THEN
       IF (z==mu) THEN
          PDF_NORMAL=NR_BIG
       ELSE 
          PDF_NORMAL=0.0d0
       END IF
    END IF
    arg = -0.5d0*(z-mu)*(z-mu)/var
    arg = -logroot2pi - 0.5*LOG(var) + arg
    PDF_Normal = EXP(arg)
  END FUNCTION PDF_Normal
  
  REAL(8) FUNCTION PDF_Truncated_Normal(z,mu,varin,limit,infinity)
    IMPLICIT NONE ! If infinity =.TRUE. (limit,inf) if infinity=.FALSE. (inf,limit)
    REAL(8), INTENT(IN) :: z,mu,varin,limit
    LOGICAL, INTENT(IN) :: infinity
    REAL(8) :: var,num,den
    var=varin
    IF (var<0.0d0) CALL NRERROR('variance has to be positive: PDF_Truncated_Normal')
    IF ((infinity).AND.(z<limit)) THEN
       PDF_Truncated_Normal = 0.0d0
       RETURN
    END IF
    IF ((.NOT.infinity).AND.(z>limit)) THEN
       PDF_Truncated_Normal = 0.0d0
       RETURN
    END IF
    IF (SQRT(var)<NR_SMALL) var=NR_SMALL
    num = PDF_Normal(z,mu,var)
    IF (num==0.0d0) THEN
       PDF_Truncated_Normal = 0.0d0
       RETURN
    END IF
    den = CDF_Normal(limit,mu,var)
    IF (infinity) den = 1.0 -den
    IF (den==0.0d0) THEN
       PDF_Truncated_Normal=num
       RETURN
    END IF
    PDF_Truncated_Normal=num/den
  END FUNCTION PDF_Truncated_Normal
  
  REAL(8) FUNCTION PDF_Double_Truncated_Normal(z,mu,varin,lo,hi)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: z,mu,varin,lo,hi
    REAL(8) :: var,num,den
    var=varin
    IF (var<0.0d0) CALL NRERROR('variance has to be positive: PDF_Double_Truncated_Normal')
    IF (lo>hi) CALL NRERROR('lower limit lower than upper limit: PDF_Double_Truncated_Normal')
    IF ((z<lo) .OR. (z>hi)) THEN
       PDF_Double_Truncated_Normal = 0.0d0
       RETURN
    END IF
    IF (SQRT(var)<NR_SMALL) var=NR_SMALL
    den = CDF_Normal(hi,mu,var) - CDF_Normal(lo,mu,var)
    num = PDF_Normal(z,mu,var)
    PDF_Double_Truncated_Normal=num/den
  END FUNCTION PDF_Double_Truncated_Normal

  REAL(8) FUNCTION PDF_Bivariate_Normal(x,mu,var,rho)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x(2),mu(2),var(2),rho
    REAL(8), PARAMETER :: twopi = 6.28318530717959D0
    REAL(8) :: arg,std(2),z(2),or2
    INTEGER :: j
    IF (ANY(var<0.0d0)) CALL NRERROR('variance has to be positive: PDF_Bivariate_Normal')
    or2=1.0d0-rho*rho
    DO j=1,2
       std(j) = SQRT(var(j))
       z(j) = (x(j)-mu(j))/std(j)
    END DO
    arg = -0.5d0*((z(1)**2)+(z(2)**2)-(2.0d0*rho*z(1)*z(2)))/or2
    PDF_Bivariate_Normal = EXP(arg)/(twopi*std(1)*std(2)*SQRT(or2))
  END FUNCTION PDF_Bivariate_Normal

  REAL(8) FUNCTION PDF_Gamma(x,a,b)
    ! Returns the pdf of a gamma distribution evaluated at x. a is shape parameter and b inverse scale. 
    ! such that mean=a/b and var=a/(b*b)
    ! That is P(x) = [(b^a)/Gamma(a)] * [x^(a-1)] * exp(-x*b) for x>0 and a>=1.
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x,a,b
    IF (a<0.0d0) CALL NRERROR('a has to be positive: PDF_Gamma')
    IF (b<0.0d0) CALL NRERROR('b has to be positive: PDF_Gamma')
    IF (x>0) THEN
       PDF_Gamma=(a-1.0d0)*LOG(x)-(x*b)-NR_gammln(a)+a*LOG(b)
       PDF_Gamma=EXP(PDF_Gamma)
    ELSE IF ((x==0.0d0).AND.(a<1.0d0)) THEN
       PDF_Gamma=NR_BIG
    ELSE IF ((x==0.0d0).AND.(a==1.0d0)) THEN
       PDF_Gamma=b
    ELSE IF ((x==0.0d0).AND.(a>1.0d0)) THEN
       PDF_Gamma=0.0d0
    ELSE IF (x<0.0d0) THEN
       PDF_Gamma=0.0d0
    ELSE
       CALL NRERROR('PDF_Gamma failed')
    END IF
  END FUNCTION PDF_Gamma
  
  REAL(8) FUNCTION PDF_Mixture_Gamma(x,a,b,p)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x,a(:),b(:),p(:)
    REAL(8) :: out
    INTEGER :: j,n
    IF (ANY(a<0.0d0)) CALL NRERROR('a has to be positive: PDF_Mixture_Gamma')
    IF (ANY(b<0.0d0)) CALL NRERROR('b has to be positive: PDF_Mixture_Gamma')
    out=0.0d0
    DO j=1,SIZE(p)
       out=out + p(j)*PDF_Gamma(x,a(j),b(j))
    END DO
    PDF_MIXTURE_GAMMA=out
  END FUNCTION PDF_Mixture_Gamma
  
  REAL(8) FUNCTION PDF_Dirichlet(x,a)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: a(:),x(:)
    REAL(8) :: num,den,aux
    INTEGER :: j
    num=EXP(NR_gammln(SUM(a)))
    den=1.0d0
    aux=1.0d0
    DO j = 1, SIZE(a)
       den=den*EXP(NR_gammln(a(j)))
       aux = aux*(x(j)**(a(j)-1.0d0))
    END DO
    PDF_Dirichlet = num*aux/den
  END FUNCTION PDF_Dirichlet

  REAL(8) FUNCTION PDF_Exponential(x,a)
    ! PDF OF EXPONENTIAL DISTRIBUTION P(x|a) = (1/a)*exp(-x/a)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x,a
    IF (a<0.0d0) CALL NRERROR('a has to be positive: PDF_Exponential')
    PDF_Exponential=EXP(-x/a)/a
  END FUNCTION PDF_Exponential

  REAL(8) FUNCTION PDF_Chi2(x,df)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x,df
    IF (df<0.0d0) CALL NRERROR('Degrees of freedom have to be positive PDF_Chi2')
    PDF_Chi2=PDF_Gamma(x,0.5d0*df,0.5d0)
  END FUNCTION PDF_Chi2

  REAL(8) FUNCTION PDF_Mixture_Normal(x,mu,var,p)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: p(:),x,mu(:),var(:)
    INTEGER :: j
    IF (ANY(var<0.0d0)) CALL NRERROR('variance has to be positive: PDF_Mixture_Normal')
    IF ((SUM(p)>1.000001d0).OR.(SUM(p)<0.999999d0)) CALL NRERROR('weights have to add to one: PDF_Mixture_Normal')
    IF (ANY(p<0.0d0)) CALL NRERROR('weights have to be positive: PDF_Mixture_Normal')
    PDF_Mixture_Normal = 0.0d0
    DO j = 1, SIZE(p)
       PDF_Mixture_Normal = PDF_Mixture_Normal + p(j)*PDF_Normal(x,mu(j),var(j))
    END DO
  END FUNCTION PDF_Mixture_Normal

  REAL(8) FUNCTION PDF_Truncated_Mixture_Normal(x,mu,var,p,limit,infinity)
    ! IF infinity=.TRUE. (limit,infinity) otherwise (-inifinity,limit)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: p(:),x,mu(:),var(:),limit
    LOGICAL, INTENT(IN) :: infinity
    REAL(8) :: w(SIZE(p))
    INTEGER :: j,n
    IF (ANY(var<0.0d0)) CALL NRERROR('variance has to be positive: PDF_Mixture_Normal')
    IF ((SUM(p)>1.000001d0).OR.(SUM(p)<0.999999d0)) CALL NRERROR('weights have to add to one: PDF_Mixture_Normal')
    IF (ANY(p<0.0d0)) CALL NRERROR('weights have to be positive: PDF_Mixture_Normal')
    n=SIZE(p)
    IF (infinity) THEN
       IF (x<limit) THEN
          PDF_TRUNCATED_Mixture_Normal=0.0d0
          RETURN
       END IF
       DO j=1,n
          w(j)=(1.0d0-CDF_NORMAL(limit,mu(j),var(j)))*p(j)
       END DO
    ELSE
       IF (x>limit) THEN
          PDF_TRUNCATED_Mixture_Normal=0.0d0
          RETURN
       END IF
       DO j=1,n
          w(j)=CDF_NORMAL(limit,mu(j),var(j))*p(j)
       END DO
    END IF
    w=p/SUM(w)
    PDF_Truncated_Mixture_Normal = 0.0d0
    DO j = 1, SIZE(p)
       PDF_Truncated_Mixture_Normal = PDF_Truncated_Mixture_Normal + w(j)*PDF_Normal(x,mu(j),var(j))
    END DO
  END FUNCTION PDF_Truncated_Mixture_Normal

  REAL(8) FUNCTION PDF_DOUBLE_Truncated_Mixture_Normal(x,mu,var,p,lo,hi)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: p(:),x,mu(:),var(:),lo,hi
    REAL(8) :: w(SIZE(p))
    INTEGER :: j,n
    IF (ANY(var<0.0d0)) CALL NRERROR('variance has to be positive: PDF_Mixture_Normal')
    IF ((SUM(p)>1.000001d0).OR.(SUM(p)<0.999999d0)) CALL NRERROR('weights have to add to one: PDF_Mixture_Normal')
    IF (ANY(p<0.0d0)) CALL NRERROR('weights have to be positive: PDF_Mixture_Normal')
    IF (hi<lo) CALL NRERROR('lower limit higher than upper limit')
    IF ((x<lo).OR.(x>hi)) THEN
       PDF_DOUBLE_TRUNCATED_Mixture_Normal=0.0d0
       RETURN
    END IF
    n=SIZE(p)
    DO j=1,n
       w(j)=( CDF_NORMAL(hi,mu(j),var(j))-CDF_NORMAL(lo,mu(j),var(j)) )*p(j)
    END DO
    w=p/SUM(w)
    PDF_Double_Truncated_Mixture_Normal = 0.0d0
    DO j = 1, SIZE(p)
       PDF_Double_Truncated_Mixture_Normal = PDF_Double_Truncated_Mixture_Normal + w(j)*PDF_Normal(x,mu(j),var(j))
    END DO
  END FUNCTION PDF_DOUBLE_Truncated_Mixture_Normal
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOG OF PDFS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(8) FUNCTION LN_PDF_Normal(z,mu,varin)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: z,mu,varin
    REAL(8), PARAMETER :: logroot2pi = 0.918938533204672780563271317078D0
    REAL(8) :: arg,var
    var=varin
    IF (var<0.0d0) CALL NRERROR('variance has to be positive: LN_PDF_Normal')
    IF (SQRT(var)<NR_SMALL) var=NR_SMALL
    arg = -0.5d0*(z-mu)*(z-mu)/var
    LN_PDF_Normal = -logroot2pi - 0.5*LOG(var) + arg
  END FUNCTION LN_PDF_Normal
  
  REAL(8) FUNCTION LN_PDF_Gamma(x,a,b)
    ! Returns the pdf of a gamma distribution evaluated at x. a is shape parameter and b inverse scale.
    ! such that mean=a/b and var=a/(b*b)
    ! That is P(x) = [(b^a)/Gamma(a)] * [x^(a-1)] * exp(-x*b) for x>0 and a>=1.
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x,a,b
    REAL(8), PARAMETER :: LN_BIG = 575.646273248511421004497863671d0
    IF (a<0.0d0) CALL NRERROR('a has to be positive: LN_PDF_Gamma')
    IF (b<0.0d0) CALL NRERROR('b has to be positive: LN_PDF_Gamma')
    IF (x>0) THEN
       LN_PDF_Gamma=(a-1.0d0)*LOG(x)-(x*b)-NR_gammln(a)+a*LOG(b)
    ELSE IF ((x==0.0d0).AND.(a<1.0d0)) THEN
       LN_PDF_Gamma=LN_BIG
    ELSE IF ((x==0.0d0).AND.(a==1.0d0)) THEN
       LN_PDF_Gamma=LOG(b)
    ELSE IF ((x==0.0d0).AND.(a>1.0d0)) THEN
       LN_PDF_Gamma=-NR_BIG
    ELSE IF (x<0.0d0) THEN
       LN_PDF_Gamma=-NR_BIG
    ELSE
       CALL NRERROR('LN_PDF_Gamma failed')
    END IF
  END FUNCTION LN_PDF_Gamma
  
  REAL(8) FUNCTION LN_PDF_Exponential(x,a)
    ! PDF OF EXPONENTIAL DISTRIBUTION P(x|a) = (1/a)*exp(-x/a)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x,a
    IF (a<0.0d0) CALL NRERROR('a has to be positive: LN_PDF_Exponential')
    LN_PDF_Exponential=(-x/a) - LOG(a)
  END FUNCTION LN_PDF_Exponential
  
  REAL(8) FUNCTION LN_PDF_Dirichlet(x,a)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: a(:),x(:)
    REAL(8) :: aux
    INTEGER :: j
    aux=NR_gammln(SUM(a))
    DO j = 1, SIZE(a)
       aux = aux - NR_gammln(a(j)) + (a(j)-1.0d0)*LOG(x(j))
    END DO
    LN_PDF_Dirichlet = aux
  END FUNCTION LN_PDF_Dirichlet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CDF's !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(8) FUNCTION CDF_Normal(zin,mu,var)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: zin
    REAL(8), INTENT(IN) :: mu,var
    REAL(8) :: zabs,p,arg,logpdf,z,std
    REAL(8), PARAMETER :: p0=220.2068679123761D0,p1=221.2135961699311D0,p2=112.0792914978709D0, &
         p3 = 33.91286607838300D0,p4 = 6.373962203531650D0,p5 = .7003830644436881D0, &
         p6 = .3526249659989109D-01,q0 = 440.4137358247522D0,q1 = 793.8265125199484D0, &
         q2 = 637.3336333788311D0,q3 = 296.5642487796737D0,q4 = 86.78073220294608D0, &
         q5=16.06417757920695D0,q6=1.755667163182642D0,q7=.8838834764831844D-1,cutoff = 7.071D0, &
         logroot2pi = 0.918938533204672780563271317078D0
    IF (var<0.0d0) CALL NRERROR('variance has to be positive: CDF_Normal')
    z=zin-mu
    zabs=ABS(z)  
    IF (zabs<NR_SMALL) THEN
       CDF_Normal=0.5d0 !Is this right?
       RETURN
    END IF
    std = SQRT(var)
    IF (std<NR_SMALL) THEN
       IF (zin-mu>=0.0d0) THEN
          CDF_Normal = 1.0d0
       ELSE IF (zin-mu<0.0d0) THEN
          CDF_Normal = 0.0d0
       END IF
    END IF
    IF (z > 37.0D0) THEN
       CDF_Normal = 1.0D0
       RETURN
    ELSE IF (z < -37.0D0) THEN
       CDF_Normal = 0.0D0
       RETURN
    END IF
    zabs=zabs/std
    arg = -0.5D0*zabs*zabs
    logpdf = -logroot2pi - LOG(std) + arg
    IF (zabs < cutoff) THEN
       p = arg + LOG(((((((p6*zabs + p5)*zabs + p4)*zabs + p3)*zabs + &
            p2)*zabs + p1)*zabs + p0)) - LOG((((((((q7*zabs + q6)*zabs + &
            q5)*zabs + q4)*zabs + q3)*zabs + q2)*zabs + q1)*zabs + &
            q0))
    ELSE
       p = logpdf - LOG((zabs + 1.0D0/(zabs + 2.0D0/(zabs + 3.0D0/(zabs + 4.0D0/ &
            (zabs + 0.65D0))))))
    END IF
    p = EXP(p)
    IF (z < 0.0D0) THEN
       CDF_Normal=p
       RETURN
    ELSE
       CDF_Normal = 1.0D0 - p
       RETURN
    END IF
    RETURN
  END FUNCTION CDF_Normal

  REAL(8) FUNCTION CDF_Truncated_Normal(x,mu,varin,limit,infinity)
    IMPLICIT NONE ! If infinity =.TRUE. (limit,inf) if infinity=.FALSE. (inf,limit)
    REAL(8), INTENT(IN) :: x,mu,varin,limit
    LOGICAL, INTENT(IN) :: infinity
    REAL(8) :: aux,var,num,den,cdf
    IF (varin<0.0) CALL NRERROR('variance has to be positive: CDF_Truncated_Normal')
    IF ((infinity).AND.(x<limit)) THEN
       CDF_Truncated_Normal=0.0d0
       RETURN
    ELSE IF ((.NOT.infinity).AND.(x>limit)) THEN
       CDF_Truncated_Normal=1.0d0
       RETURN
    ELSE
       var = varin
       IF (SQRT(varin)<NR_SMALL) var = NR_SMALL
       den = CDF_Normal(limit,mu,var)
       cdf = den
       num = CDF_Normal(x,mu,var)
       IF (infinity) THEN
          den = 1.0 - cdf 
          num = num - cdf
       END IF
       CDF_Truncated_Normal=(num/den)
       RETURN
    END IF
  END FUNCTION CDF_Truncated_Normal

  REAL(8) FUNCTION CDF_Double_Truncated_Normal(x,mu,varin,lo,hi)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x,mu,varin,lo,hi
    REAL(8) :: aux,var,num,den,cdf
    IF (varin<0.0) CALL NRERROR('variance has to be positive: CDF_Double_Truncated_Normal')
    IF (lo>hi) CALL NRERROR('lower limit larger than upper limit: CDF_Double_Truncated_Normal')
    IF (x<lo) THEN
       CDF_Double_Truncated_Normal = 0.0d0
       RETURN
    ELSE IF (x>hi) THEN
       CDF_Double_Truncated_Normal = 0.0d0
       RETURN
    ELSE 
       var=varin
       IF (SQRT(varin)<NR_SMALL) var=NR_SMALL
       cdf = CDF_Normal(lo,mu,var)
       den = CDF_Normal(hi,mu,var)-cdf
       num = CDF_Normal(x,mu,var)-cdf
       CDF_Double_Truncated_Normal = num/den
       RETURN
    END IF
  END FUNCTION CDF_Double_Truncated_Normal
  
  REAL(8) FUNCTION CDF_Truncated_Mixture_Normal(x,mu,var,p,limit,infinity)
    ! IF infinity=.TRUE. (limit,infinity) otherwise (-inifinity,limit)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: p(:),x,mu(:),var(:),limit
    LOGICAL, INTENT(IN) :: infinity
    REAL(8) :: w(SIZE(p)),cdf
    INTEGER :: j,n
    IF (ANY(var<0.0d0)) CALL NRERROR('variance has to be positive: PDF_Mixture_Normal')
    IF ((SUM(p)>1.000001d0).OR.(SUM(p)<0.999999d0)) CALL NRERROR('weights have to add to one: PDF_Mixture_Normal')
    IF (ANY(p<0.0d0)) CALL NRERROR('weights have to be positive: PDF_Mixture_Normal')
    n=SIZE(p)
    IF (infinity) THEN
       IF (x<limit) THEN
          CDF_TRUNCATED_Mixture_Normal=0.0d0
          RETURN
       END IF
       DO j=1,n
          w(j)=(1.0d0-CDF_NORMAL(limit,mu(j),var(j)))*p(j)
       END DO
    ELSE
       IF (x>limit) THEN
          CDF_TRUNCATED_Mixture_Normal=1.0d0
          RETURN
       END IF
       DO j=1,n
          w(j)=CDF_NORMAL(limit,mu(j),var(j))*p(j)
       END DO
    END IF
    w=p/SUM(w)
    CDF_Truncated_Mixture_Normal = 0.0d0
    DO j = 1, SIZE(p)
       IF (infinity) THEN
          cdf=CDF_Normal(x,mu(j),var(j))-CDF_Normal(limit,mu(j),var(j))
       ELSE
          cdf=CDF_Normal(x,mu(j),var(j))
       END IF
       CDF_Truncated_Mixture_Normal = CDF_Truncated_Mixture_Normal + w(j)*cdf
    END DO
  END FUNCTION CDF_Truncated_Mixture_Normal

  REAL(8) FUNCTION CDF_Gamma(x,a,b)
    ! Gives me the gamma cdf
    ! such that mean=a/b and var=a/(b*b)
    ! That is F(x) = integral(0,x) of [(b^a)/Gamma(a)] * [x^(a-1)] * exp(-x*b) for x>0 and a>=1.
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x,a,b
    IF (a<0.0d0) CALL NRERROR('a has to be positive CDF_Gamma')
    IF (b<0.0d0) CALL NRERROR('b has to be positive CDF_Gamma')
    CDF_Gamma=NR_gamminc(x*b,a)
    IF (CDF_Gamma>1.0d0) CDF_Gamma=1.0d0 
  END FUNCTION CDF_Gamma

  REAL(8) FUNCTION CDF_Exponential(x,a)
    ! CDF OF EXPONENTIAL DISTRIBUTION P(x|a) = (1/a)*exp(-x/a)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x,a
    IF (a<0.0d0) CALL NRERROR('a has to be positive: CDF_Exponential')
    CDF_Exponential=1.0d0-EXP(-x/a)
  END FUNCTION CDF_Exponential

  REAL(8) FUNCTION CDF_Chi2(x,df)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x,df
    IF (df<0.0d0) CALL NRERROR('Degrees of freedom have to be positive: CDF_Chi2')
    CDF_Chi2=CDF_Gamma(x,0.5d0*df,0.5d0)
  END FUNCTION CDF_Chi2

  REAL(8) FUNCTION CDF_Mixture_Normal(x,mu,var,p)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: p(:),x,mu(:),var(:)
    INTEGER :: j
    IF (ANY(var<0.0d0)) CALL NRERROR('variance has to be positive: CDF_Mixture_Normal')
    IF ((SUM(p)>1.000001d0).OR.(SUM(p)<0.999999d0)) CALL NRERROR('weights have to add to one: CDF_Mixture_Normal')
    IF (ANY(p<0.0d0)) CALL NRERROR('weights have to be positive: CDF_Mixture_Normal')
    CDF_Mixture_Normal = 0.0D0
    DO j = 1, SIZE(p)
       CDF_Mixture_Normal = CDF_Mixture_Normal + p(j)*CDF_Normal(x,mu(j),var(j))
    END DO
  END FUNCTION CDF_Mixture_Normal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INVERSE CDF's !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(8) FUNCTION CDF_Normal_Inverse(P,mu,varin)
    !	Produces the normal deviate Z corresponding to a given lower
    !	tail area of P; Z is accurate to about 1 part in 10**16.
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: P,mu,varin
    REAL(8) :: Q,R,var
    REAL(8), PARAMETER :: ZERO=0.0D0, ONE = 1.0D0, HALF = 0.5D0,SPLIT1 = 0.425D0, SPLIT2 = 5.D0, &
         CONST1 = 0.180625D0, CONST2 = 1.6D0
    !	Coefficients for P close to 0.5
    REAL(8), PARAMETER :: A0 = 3.3871328727963666080D0,A1=1.3314166789178437745D+2, &
         A2=1.9715909503065514427D+3,A3=1.3731693765509461125D+4,A4=4.5921953931549871457D+4, &
         A5=6.7265770927008700853D+4,A6=3.3430575583588128105D+4,A7=2.5090809287301226727D+3, &
         B1=4.2313330701600911252D+1,B2=6.8718700749205790830D+2,B3=5.3941960214247511077D+3, &
         B4=2.1213794301586595867D+4,B5=3.9307895800092710610D+4,B6=2.8729085735721942674D+4, &
         B7=5.2264952788528545610D+3
    !	Coefficients for P not close to 0, 0.5 or 1.
    REAL(8), PARAMETER :: C0=1.42343711074968357734D0,C1=4.63033784615654529590D0, &
         C2=5.76949722146069140550D0,C3=3.64784832476320460504D0,C4=1.27045825245236838258D0, &
         C5=2.41780725177450611770D-1,C6=2.27238449892691845833D-2,C7=7.74545014278341407640D-4, &
         D1=2.05319162663775882187D0,D2=1.67638483018380384940D0,D3=6.89767334985100004550D-1, &
         D4=1.48103976427480074590D-1,D5=1.51986665636164571966D-2,D6=5.47593808499534494600D-4, &
         D7=1.05075007164441684324D-9
    !	CoefficientsforP near 0 or 1.
    REAL(8), PARAMETER :: E0=6.65790464350110377720D0,E1=5.46378491116411436990D0, &
         E2=1.78482653991729133580D0,E3=2.96560571828504891230D-1,E4=2.65321895265761230930D-2, &
         E5=1.24266094738807843860D-3,E6=2.71155556874348757815D-5,E7=2.01033439929228813265D-7, &
         F1=5.99832206555887937690D-1,F2=1.36929880922735805310D-1,F3=1.48753612908506148525D-2, &
         F4=7.86869131145613259100D-4,F5=1.84631831751005468180D-5,F6=1.42151175831644588870D-7, &
         F7=2.04426310338993978564D-15
    IF ((P>1.0D0).OR.(P<0.0D0)) CALL NRERROR('P has to be between zero and one: CDF_Normal_Inverse')
    IF (varin<0.0d0) CALL NRERROR('variance has to be positive: CDF_Normal_Inverse')
    var=varin
    Q = P - HALF
    IF (ABS(Q) .LE. SPLIT1) THEN
       R = CONST1 - Q * Q
       CDF_Normal_Inverse=Q*(((((((A7*R+A6)*R+A5)*R+A4)*R+A3)*R+A2)*R+A1)*R+ A0) / &
            (((((((B7 * R + B6) * R + B5) * R + B4) * R + B3)* R + B2) * R + B1) * R + ONE)
       CDF_Normal_Inverse=(SQRT(var)*CDF_Normal_Inverse)+mu
       RETURN
    ELSE
       IF (Q .LT. ZERO) THEN
          R = P
       ELSE
          R = ONE - P
       END IF
       R = SQRT(-LOG(R))
       IF (R .LE. SPLIT2) THEN
          R = R - CONST2
          CDF_Normal_Inverse=(((((((C7*R + C6) * R + C5) * R + C4) * R + C3)*R+C2)*R+ C1) * R + C0) / &
               (((((((D7 * R + D6) * R + D5) * R + D4) * R + D3)* R + D2) * R + D1) * R + ONE)
       ELSE
          R = R - SPLIT2
          CDF_Normal_Inverse=(((((((E7*R + E6) * R + E5) * R + E4) * R + E3)*R+E2)*R + E1) * R + E0) / &
               (((((((F7 * R + F6) * R + F5) * R + F4) * R + F3)* R + F2) * R + F1) * R + ONE)
       END IF
       IF (Q .LT. ZERO) CDF_Normal_Inverse = - CDF_Normal_Inverse
       CDF_Normal_Inverse=(SQRT(var)*CDF_Normal_Inverse)+mu
       RETURN
    END IF
  END FUNCTION CDF_Normal_Inverse

  REAL(8) FUNCTION CDF_Gamma_Inverse(p,a,b)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: p,a,b
    REAL(8) :: xk,eps,mn,v,temp,sigma,mu,h,xnew
    INTEGER :: j
    IF ((p<0.0d0).OR.(p>1.0d0)) CALL NRERROR('P has to be between 0 and 1: CDF_Gamma_Inverse')
    IF (a<0.0d0) CALL NRERROR('a has to be positive: CDF_Gamma_Inverse')
    IF (b<0.0d0) CALL NRERROR('b has to be positive: CDF_Gamma_Inverse')
    IF (p==0.0d0) THEN
       CDF_Gamma_Inverse=0.0d0
       RETURN
    ELSE IF (p==1.0d0) THEN
       CDF_Gamma_Inverse=NR_Big
       RETURN
    ELSE
       eps=2.220446049250313D-016
       mn=a*b
       v=mn*b
       temp=LOG(v+(mn**2)) 
       mu=2.0d0*LOG(mn)-0.5d0*temp
       sigma=-2.0d0*LOG(mn)+temp
       xk=EXP(CDF_Normal_Inverse(p,mu,sigma))
       h = 1.0d0
       j=0
       DO WHILE ((ABS(h)>SQRT(eps)*ABS(xnew)).AND.(ABS(h)>SQRT(eps)))!.AND.(j<200))
          j=j+1     
          h=(CDF_Gamma(xk,a,b) - p) / PDF_Gamma(xk,a,b)
          xnew = xk - h
          IF (xnew<0.0d0) THEN
             xnew = xk/10.d0
             h = xk-xnew
          END IF
          xk=xnew
       END DO
       CDF_Gamma_Inverse=xk
    END IF
  END FUNCTION CDF_Gamma_Inverse
  
  REAL(8) FUNCTION CDF_Exponential_Inverse(p,a)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: p,a
    IF ((p<0.0d0).OR.(p>1.0d0)) CALL NRERROR('P has to be between 0 and 1: CDF_Exponential_Inverse')
    IF (a<0.0d0) CALL NRERROR('a has to be positive: CDF_Exponential_Inverse')
    CDF_Exponential_Inverse=-LOG(1.0d0-p)*a
  END FUNCTION CDF_Exponential_Inverse

  REAL(8) FUNCTION CDF_Chi2_Inverse(p,df)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: p,df
    IF ((p<0.0d0).OR.(p>1.0d0)) CALL NRERROR('P has to be between 0 and 1: CDF_Chi2_Inverse')
    IF (df<0.0d0) CALL NRERROR('Degrees of freedom have to be positive: CDF_Chi2_Inverse')
    CDF_Chi2_Inverse=CDF_Gamma_Inverse(p,0.5d0*df,0.5d0)
  END FUNCTION CDF_Chi2_Inverse

  REAL(8) FUNCTION CDF_Mixture_Normal_Inverse(P,mu,var,prob)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: P,mu(:),var(:),prob(:)
    REAL(8) :: x1,x2,f1,f2,xm,rtbis,dx
    INTEGER :: n,j
    n=SIZE(mu)
    IF ((p<0.0d0).OR.(p>1.0d0)) CALL NRERROR('P has to be between 0 and 1: CDF_Mixture_Normal_Inverse')
    IF (ANY(var<0.0d0)) CALL NRERROR('variance has to be positive: CDF_Mixture_Normal_Inverse')
    IF (p==0.0d0) THEN
       CDF_Mixture_Normal_Inverse=0.0d0
       RETURN
    ELSE IF (p==1.0d0) THEN
       CDF_Mixture_Normal_Inverse=NR_Big
       RETURN
    ELSE
       x1 = -NR_Big
       x2 = 0.0d0
       f1 = CDF_Mixture_Normal(x1,mu,var,prob) - p
       f2 = CDF_Mixture_Normal(x2,mu,var,prob) - p
       IF (f2==0.0d0) THEN
          CDF_Mixture_Normal_Inverse=x2
          RETURN
       END IF
       IF (f1*f2>0.0d0) THEN
          f1=f2
          x1 = 0.0d0
          x2 = NR_Big
       END IF
       IF (f1<0.0d0) THEN
          CDF_Mixture_Normal_Inverse=x1
          dx=x2-x1
       ELSE
          CDF_Mixture_Normal_Inverse=x2
          dx=x1-x2
       END IF
       DO 
          dx=dx*0.5d0
          x2=CDF_Mixture_Normal_Inverse+dx
          f2=CDF_Mixture_Normal(x2,mu,var,prob) - p
          IF (f2<=0.0d0) CDF_Mixture_Normal_Inverse=x2
          IF ((ABS(dx)<0.0000000001d0).OR.(f2==0.0d0)) RETURN
       END DO
    END IF
  END FUNCTION CDF_Mixture_Normal_Inverse
  
END MODULE PROBABILITY
