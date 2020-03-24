!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!	program correlation_celegans.f90
!	The network contains 248 neurons and 511 gap junctions. 
! 	Calculates the correlation matrix and the total order parameters (r_total, dpsi_av) of network.
! 	Fixed parameters: force intensity, frequency of external force and internal coupling.
!	Choose a distribution: idist = 1 (uniform); idist = 2 (triangular); idist = 3 (Lorenztian); idist = 4 (Gaussian).
!	
!	INPUT:	
!	Reads the weighted adjacency matrix (netname) and communities (communities).
!
!	OUTPUT
! 	1. filenamer: order parameters for each module.
!	2. matrixcorr: correlation matrix (binary), ordered by modules; 248x248.
!
!	Carolina A. Moreira e Marcus A. M. de Aguiar - 22/10/2018.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! Module defining global variables
MODULE globals
IMPLICIT REAL*8(A-H,O-Z)
INTEGER(4), SAVE :: nosc, nfrac
REAL*8, SAVE :: avk,lambdahat,sigma,fhat, sigmahat,pi
REAL*8, ALLOCATABLE, SAVE :: omega(:),tstep(:)
INTEGER, ALLOCATABLE, SAVE :: a(:,:),vk(:)
INTEGER, ALLOCATABLE, SAVE :: g(:,:),modnode(:),modsize(:)
INTEGER, SAVE :: modmax,tmax,ntime,nmod,modforce
REAL*8, ALLOCATABLE, SAVE :: zx(:,:),zy(:,:),corr(:,:),corrdiag(:)
REAL*8, ALLOCATABLE, SAVE :: rmi(:,:),dpsimi(:,:),dpsierri(:,:) ! variable inter modulos
REAL*8, ALLOCATABLE, SAVE :: avdegreesize(:)
REAL*8, ALLOCATABLE, SAVE :: r_total,dpsi_av,dpsi_av_err
END MODULE

PROGRAM kuramoto
USE globals
IMPLICIT REAL*8(A-H,O-Z)
INTEGER :: iseed(12),idist
REAL*8, ALLOCATABLE :: dy(:),y(:),yscal(:),yp(:)
REAL*8 adist,xf,lambda,f,fi,ff
 CHARACTER*30 netname,communities
 CHARACTER*60 filenamer,matrixcorr
pi = dacos(-1.0d0)
pi2 = 2.0*pi

!------ PARAMETERS -------- 
 idist = 4			! Gaussian distribution
 xf = 100				! Total time
 ftempo = 0.5			! Fraction of total time
 nosc = 248			! Number of nodes
 modmax = 20			! Max number of modules
 tmax = 10000			! Max size of vectors
 lambda = 100.0		! Internal Kuramoto coupling (fixed)
 f = 17.0				! External forcing (fixed)
 modforce = 1			! Index of forced module

 adist = 1.0d0			! Width of gaussian distribution
 sigma = 3.0			! Frequency of external forcing (fixed)
 anosc = dfloat(nosc)
 eps=1.E-07

! input
 netname = 'EJ248_adj.txt'
 communities = 'communities.dat'

! output
 filenamer = 'namer.dat'
 matrixcorr =  'namecorr.dat'

ALLOCATE (dy(nosc),y(nosc),yscal(nosc),yp(nosc))
ALLOCATE (omega(nosc))
ALLOCATE (a(nosc,nosc), vk(nosc))

130 FORMAT(249(f6.3,1x))
131 FORMAT(2(f6.3,2x),i2)

! reads adjacency matrix
 OPEN(UNIT=11,FILE=netname,STATUS='OLD')
    do i=1,nosc
        read(11,900) (a(i,j),j=1,nosc)
    end do
 CLOSE(11)
900 FORMAT(2000(i2,1x))

OPEN(UNIT=56,file=matrixcorr,status='unknown')

avk = dfloat(sum(a))/anosc
vk = SUM(a,DIM=1)                         ! Weighted degree vector.
kmax = maxval(vk)

! Read communities
OPEN(UNIT=11,file=communities,status='old')
ALLOCATE (g(modmax,nosc),modnode(nosc),modsize(modmax))
im = 0
ic = 0
DO WHILE (ic < nosc)
    im = im + 1
    READ(11,*) modsize(im)
    READ(11,*) (g(im,i),i=1,modsize(im))
    DO i=1,modsize(im)
        modnode(g(im,i)) = im
    END DO
    ic = ic + modsize(im)
END DO
 CLOSE(11)
nmod = im

! Average degree of modules
ALLOCATE (avdegreesize(nmod))
avdegreesize = 0.d0
DO i=1,nmod
    DO j=1,modsize(i)
        DO k=1,nosc
            avdegreesize(i) = avdegreesize(i) + a(g(i,j),k)
        END DO
    END DO
    avdegreesize(i) = avdegreesize(i)/float(modsize(i))
END DO
write(6,*) 'average degree of modules - including weight'
write(6,*) avdegreesize
write(6,*)
write(6,*) 'global average degree'
write(6,*) float(sum(a))/float(nosc)


IF(idist==1) THEN
         OPEN(UNIT=10,FILE='uforced.dat',STATUS='UNKNOWN')      ! Uniform distribution.
ELSE IF (idist==2) THEN
         OPEN(UNIT=10,FILE='tforced.dat',STATUS='UNKNOWN')      ! Triangular distribution.
ELSE IF (idist==3) THEN
         OPEN(UNIT=10,FILE='lforced_l1e1.dat',STATUS='UNKNOWN') ! Lorentizian distribution.
ELSE
         OPEN(UNIT=50,FILE= filenamer, STATUS='UNKNOWN')        ! Gaussian distribution.
END IF

! Initialize random number generator
OPEN(UNIT=20,FILE="seed.in",STATUS='OLD')
    READ(20,*) iseed
CLOSE(20)
CALL RANDOM_SEED(put=iseed)
CALL RANDOM_NUMBER(aux)

lambdahat = lambda/adist
sigmahat = sigma/adist

    sumomega = 0.0
    DO i=1,nosc
        CALL RANDOM_NUMBER(aux)
        IF(idist==1) THEN
            CALL uniform(adist,aux,ya)
        ELSE IF (idist==2) THEN
            CALL triangle(adist,aux,ya)
        ELSE IF (idist==3) THEN
            CALL lorentz(adist,aux,ya)
        ELSE
            CALL gaussian(adist,aux,ya)
        END IF
        omega(i) = ya
        sumomega = sumomega + omega(i)
    END DO

    sumomega = sumomega/anosc

    DO i=1,nosc
        omega(i) = omega(i) - sumomega - sigmahat
    END DO

    omegak=0.0
    do i = 1,nosc
        omegak = omegak + vk(i)*(omega(i) + sigmahat)
    end do
    omegak = omegak/(anosc*avk)

    xft = (1.0d0 - ftempo)*xf
    
    ALLOCATE(tstep(tmax),zx(modmax,tmax),zy(modmax,tmax))
    ALLOCATE (rmi(nmod,nmod),dpsimi(nmod,nmod),dpsierri(nmod,nmod))
    ALLOCATE(corr(nosc,nosc),corrdiag(nosc))
    
    fhat = f/adist
    
    x = 0.0d0
        
         DO i=1,nosc
             CALL RANDOM_NUMBER(aux)
             y(i)=aux*pi2
        END DO
    
        htry=0.05
        it = 0
        zx = 0.D0
        zy = 0.D0
        corr = 0.D0
        corrdiag = 0.D0

   loop_time: DO WHILE (x < xf) 
            DO j=1,nosc
                yscal(j)=y(j)
                yp(j)=y(j)
                IF(yscal(j) == 0.0) yscal(j)=0.01
            END DO
            CALL DERIV(X,Y,DY)
            CALL RKQS(y,dy,nosc,x,htry,eps,yscal,hdid,hnext,DERIV)
            IF(X+HNEXT.EQ.X) THEN
                print*, 'Stepsize not significant in RKDUMB.'
            END IF
            HTRY=HNEXT

            IF(x > xft) THEN
                it = it+1
                tstep(it) = hdid
                
                CALL ZCALC(y,it)
                CALL CORRCALC(y,yp,it)
            
            END IF
    END DO loop_time
    CLOSE(32)
    ntime = it

    CALL RPSIPONTOINTER

     ! Divide by total time
     corrdiag = corrdiag/sum(tstep)
     corr = corr/sum(tstep)

	! Non-normalized correlation
	DO i=1,nosc
        DO j=i,nosc
            corr(i,j) = corr(i,j) - corrdiag(i)*corrdiag(j)
            corr(j,i) = corr(i,j)
        END DO
     END DO

	! Normalization
	DO i=1,nosc
        DO j=i+1,nosc
            corr(i,j) = corr(i,j)/dsqrt(corr(i,i)*corr(j,j))
            corr(j,i) = corr(i,j)
        END DO
        corr(i,i) = 1.0D0
    END DO


    DO i=1,nmod
        DO j=1,modsize(i)
            write(56,101) ((corr(g(i,j),g(il,jl)),jl=1,modsize(il)),il=1,nmod )
        END DO
    END DO 
    
    write(6,*) 'force = ', f

DO i1=1,nmod
    write(50,100) (rmi(i1,l),l=1,nmod)
END DO

CALL RANDOM_SEED(get=iseed)
OPEN(UNIT=20,FILE="seed.in",STATUS='OLD', POSITION='REWIND')
WRITE (20,*) iseed
CLOSE(20)

100 FORMAT(248(F8.4,1x))
101 FORMAT(248(f6.3,1x))
200 FORMAT(F6.2, 1x, 22(F6.2, 1x))
END PROGRAM kuramoto


! Calculates zx, zy for each module
SUBROUTINE ZCALC(y,it)
USE globals
IMPLICIT REAL*8(A-H,O-Z)
REAL*8 y(nosc)
DO i=1,nosc
    k = modnode(i)
    zx(k,it) = zx(k,it) + cos(y(i))
    zy(k,it) = zy(k,it) + sin(y(i))
END DO
END

! Calculates correlation between nodes i and j
SUBROUTINE CORRCALC(y,yp,it)
USE globals
IMPLICIT REAL*8(A-H,O-Z)
REAL*8 y(nosc),yp(nosc),ppp(nosc)
DO i=1,nosc
    dpsi =  y(i)-yp(i)
    IF(dpsi > 2.0D0)  dpsi = dpsi - pi
    IF(dpsi < -2.0D0) dpsi = dpsi + pi
    ppp(i) = dpsi/tstep(it)
END DO
DO i=1,nosc
	DO j=i,nosc
	    corr(i,j) = corr(i,j) + ppp(i)*ppp(j)*tstep(it)
	    corr(j,i) = corr(i,j)
	END DO
    corrdiag(i) = corrdiag(i) + ppp(i)*tstep(it)
END DO
END

! Calculates average order parameters inter-modules
SUBROUTINE RPSIPONTOINTER
USE globals
IMPLICIT REAL*8(A-H,O-Z)
REAL*8 vdpsi(ntime),psi_tot(ntime),dpsi_tot(ntime)
REAL*8 zx_tot(ntime),zy_tot(ntime),r_tot(ntime)
rmi = 0.0D0
dpsimi = 0.0D0
dpsierri = 0.0D0
aux = dfloat(ntime-1)
antime = dfloat(ntime)
anosc = dfloat(nosc)
DO i1=1,nmod
    DO i2=i1,nmod
        vdpsi = 0.0D0
        psiold = datan( (zy(i1,1)+zy(i2,1))/(zx(i1,1)+zx(i2,1)) )
        DO j=2,ntime
            rmi(i1,i2) = rmi(i1,i2) + dsqrt((zx(i1,j)+zx(i2,j))**2 + (zy(i1,j)+zy(i2,j))**2 )
            psi = datan((zy(i1,j)+zy(i2,j))/(zx(i1,j)+zx(i2,j)))
            dpsi = psi - psiold
            IF(dpsi > 2.0D0)  dpsi = dpsi - pi
            IF(dpsi < -2.0D0) dpsi = dpsi + pi
            dpsi = dpsi/tstep(j)
            dpsimi(i1,i2) = dpsimi(i1,i2) + dpsi
            vdpsi(j) = dpsi
            psiold = psi
        END DO
        rmi(i1,i2) = rmi(i1,i2)/(aux*dfloat(modsize(i1)+modsize(i2)))
        
        dpsimi(i1,i2) = dpsimi(i1,i2)/aux
        DO j=2,ntime
            dpsierri(i1,i2) = dpsierri(i1,i2) + (vdpsi(j) - dpsimi(i1,i2))**2
        END DO
        dpsierri(i1,i2) = dsqrt(dpsierri(i1,i2)/aux)
        rmi(i2,i1) = rmi(i1,i2)
        dpsimi(i2,i1) = dpsimi(i1,i2)
        dpsierri(i2,i1) = dpsierri(i1,i2)
    END DO
END DO

! r_total and dpsi_total
zx_tot = sum(zx,dim=1)/anosc    		 ! sum over all modules
zy_tot = sum(zy,dim=1)/anosc
r_tot = dsqrt( zx_tot**2 + zy_tot**2)  	 ! r_tot(t)
r_total = sum(r_tot)/antime              ! r_tot averaged over time

psi_tot = datan(zy_tot/zx_tot)  		 ! psi_tot(t)
dpsi_tot = 0.0D0                 		 ! psi_ponto_tot(t)
dpsi_av = 0.0D0                 		 ! psi_ponto averaged over time
DO j=2,ntime
    dpsi = psi_tot(j) - psi_tot(j-1)
    IF(dpsi > 2.0D0)  dpsi = dpsi - pi
    IF(dpsi < -2.0D0) dpsi = dpsi + pi
    dpsi = dpsi/tstep(j)
    dpsi_tot(j) = dpsi
    dpsi_av = dpsi_av + dpsi
END DO
dpsi_av = dpsi_av/aux
dpsi_av_err = 0.0D0              		 ! error of psi_ponto_tot

DO j=2,ntime
    dpsi_av_err = dpsi_av_err + (dpsi_tot(j)-dpsi_av)**2
END DO
dpsi_av_err = dsqrt( dpsi_av_err/aux )

100 FORMAT(248(F8.4,1x))
write(6,*)
write(6,*) 'Average of r = ', r_total
write(6,*) 'Average of dot-psi = ', dpsi_av

END

! Equations of Kuramoto model
SUBROUTINE DERIV(x,y,dy)
USE globals
IMPLICIT REAL*8(A-H,O-Z)
REAL*8 x,y(nosc),dy(nosc)
REAL*8 couplings(nosc),aux1
INTEGER iforce(nosc)

iforce = 0
DO i=1,nosc
    IF( modnode(i) == modforce) iforce(i) = 1
END DO

couplings = 0.0
DO i=1,nosc
    aux1 = lambdahat/dfloat(vk(i))
    DO j=1,nosc
        couplings(i) = couplings(i) + a(i,j)*sin(y(j)-y(i))
    END DO
    dy(i) = omega(i) + aux1*couplings(i) - fhat*sin(y(i))*iforce(i)
END DO

RETURN
END

! Uniform: idist = 1
SUBROUTINE uniform(a,y,x)
IMPLICIT REAL*8(A-H,O-Z)
REAL*8 a,y,x
x = a*(2.0*y-1.0)
RETURN
END

! Triangle: idist = 2
SUBROUTINE triangle(a,y,x)
IMPLICIT REAL*8(A-H,O-Z)
REAL*8 a,y,x
IF(y < 0.5) THEN
           x = -a + a*sqrt(2.0*y)
    ELSE
           x = a - a*sqrt(2.0*(1.0-y))
END IF
RETURN
END

! Lorentz: idist = 3
SUBROUTINE lorentz(a,y,x)
IMPLICIT REAL*8(A-H,O-Z)
REAL*8 a,y,x
pi = dacos(-1.0d0)
x = a/(tan(pi*y))
RETURN
END

! Gaussian: idist = 4
SUBROUTINE gaussian(a,y,x)
IMPLICIT REAL*8(A-H,O-Z)
REAL*8 a,y,x,error,b
REAL*8 x0,delta_x,error_x
REAL*8 F,rho,anorm

error=0.001
pi = dacos(-1.0d0)
b = 1.0/(2.0*a*a)
anorm = dsqrt(b/pi)
x0 = 0.0
F = 0.5*(1+erf(x0*dsqrt(b))) - y
rho = anorm*exp(-b*x0*x0)

DO WHILE(abs(F) > error)
    delta_x = - F/rho
    x = x0 + delta_x
    F = 0.5*(1+erf(x*sqrt(b))) - y
    rho = anorm*exp(-b*x*x)
    x0 = x
END DO

RETURN
END

! Subroutine from numerical recipies
      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      IMPLICIT REAL*8(A-H,O-Z)
      dimension dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=2000)
      dimension yerr(NMAX),ytemp(NMAX)
      PARAMETER (SAFETY=0.9D0,PGROW=-.2D0,PSHRNK=-.25D0,ERRCON=1.89D-4)
      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.D0
      do i=1,n
         errmax=max(errmax,abs(yerr(i)/yscal(i)))
      end do
      errmax=errmax/eps
      if(errmax.gt.1.D0)then
         h=SAFETY*h*(errmax**PSHRNK)
         if(h.lt.0.1D0*h)then
            h=.1D0*h
         endif
         xnew=x+h
         if(xnew.eq.x) then
            write(6,*) 'stepsize underflow in rkqs'
            stop
         end if
         goto 1
      else
         if(errmax.gt.ERRCON)then
            hnext=SAFETY*h*(errmax**PGROW)
         else
            hnext=5.D0*h
         endif
         hdid=h
         x=x+h
         do i=1,n
            y(i)=ytemp(i)
         end do
 20      return
      endif
      END

! Subroutine from numerical recipies
      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
      IMPLICIT REAL*8(A-H,O-Z)
      dimension dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=2000)
      dimension ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),ytemp(NMAX)
      PARAMETER (A2=.2D0,A3=.3D0,A4=.6D0,A5=1.D0,A6=.875D0,&
     & B21=.2D0,B31=3.D0/40.D0,B32=9.D0/40.D0,B41=.3D0,B42=-.9D0,&
     & B43=1.2D0,B51=-11.D0/54.D0,B52=2.5D0,B53=-70.D0/27.D0,&
     & B54=35.D0/27.D0,B61=1631.D0/55296.D0,B62=175.D0/512.D0,&
     & B63=575.D0/13824.D0,B64=44275.D0/110592.D0,B65=253.D0/4096.D0,&
     & C1=37.D0/378.D0,C3=250.D0/621.D0,C4=125.D0/594.D0,&
     & C6=512.D0/1771.D0,DC1=C1-2825.D0/27648.D0,&
     & DC3=C3-18575.D0/48384.D0,DC4=C4-13525.D0/55296.D0,&
     & DC5=-277.D0/14336.D0,DC6=C6-.25D0)
      do i=1,n
         ytemp(i)=y(i)+B21*h*dydx(i)
      end do
      call derivs(x+A2*h,ytemp,ak2)
      do i=1,n
         ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
      end do
      call derivs(x+A3*h,ytemp,ak3)
      do i=1,n
         ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
      end do
      call derivs(x+A4*h,ytemp,ak4)
      do i=1,n
         ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
      end do
      call derivs(x+A5*h,ytemp,ak5)
      do i=1,n
         ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+ B65*ak5(i))
      end do
      call derivs(x+A6*h,ytemp,ak6)
      do i=1,n
         yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
      end do
      do i=1,n
         yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
      end do
      return
      END
