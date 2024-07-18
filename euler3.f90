
! *************************************************************************************
! Lab for Msc course 2016
! Programmer V.A. Titarev/P. Tsoutsanis
! *************************************************************************************

! Finite-volume scheme with TVD Runge-Kutta time stepping
! for one-dimensional compressible Euler equations
! Use third order TVD Runge-Kutta in time

! Reference solution for test problem one is given in ref1.out


 ! Declaration of variables
 
 IMPLICIT NONE

 ! Spatial order
 Integer :: SpatialOrder ! 1 or 2
 
 ! Flux type
 Integer FluxType

 ! Courant number
 Real  CFL

 ! Type of initial condition, number of spatial cells
 Integer  IvType, N

 ! Vector of conservative variables CSV(1:4,-6+n+5), first index is the Runge-Kutta stage
 Real, ALLOCATABLE::  CSV(:,:,:) 
 ! Primitive variables, no RK stage, W = (rho,u,P)
 Real, ALLOCATABLE::  PV(:,:,:)
 ! Intercell fluxes
 Real, ALLOCATABLE::  IntercellFlux(:,:,:)
 

 ! other variables
 Integer  ::  OutFreq = 25, Rkstage=0
 Real     LB, RB, h,Time
 Real, ALLOCATABLE::  CELL_CENTERS(:)
 Real, parameter :: GM=1.4
 Real DT,T_
 Integer IT
 
 ! ---------------------------------- START THE PROGRAM----------------------

 
 ! initialise of variables, grid etc
  CALL INITALL
  
 ! START TIME CYCLE
 IT=0
 do  !
   ! Compute a stable time step
   CALL ComputeTimeStep(dt)
   ! Use third-order TVD Runge-Kutta
   Call ThirdOrderTVD
   ! advance time and time counter
   T_ = T_  + DT
   IT = IT+1

   If ( MOD(IT,OutFreq) ==0) PRINT*,' it= ',it,'  t= ',T_
   If (mod(it,2) .eq. 0)  Call OutputTecplot 
   
   ! check whether we reached the output time
   If ( ABS(T_ - TIME)/TIME .LE. 1D-8) GOTO 101
 enddo

 101 CONTINUE
 Call Output
 Call OutputTecplot 
 close(111) ! close the movie file
 print*,' Job finished.'
 Print*,' Number of time steps : ',it
 PRINT*,' The end' 


 ! //////// code's subroutines ///////////////
 
 Contains 

 
 !%%%%%%%%%% initialization of the run %%%%%%%%%%
 Subroutine InitAll
  Integer I
  Real X,U1,U2,U3

 ! read the input file
  Open(1,file='euler.ini')
   Read(1,*) SpatialOrder
   Read(1,*) IvType
   Read(1,*) N
   Read(1,*) CFL
   Read(1,*) FluxType
  Close(1)


 SELECT Case(IVTYPE)
 Case(1) 
  LB = 0.
  RB = +1. 
  TIME = 0.2
 Case(2)
  LB = -5.
  RB = +5.
 Time = 5.d0   
 Case Default
  print*,' Wrong test problem number. Stop the code.'
  stop
 END SELECT

 ! Spatial Cell size
 h = (RB-LB)/N

 ! Vector of conservative variables QC(1:4,-6+n+5), first index is the Runge-Kutta stage
 ALLOCATE(CSV(1:4,3,-7:n+6),PV(4,1:4,-7:n+6),InterCellFlux(1:4,3,-1:n+1))
 ALLOCATE(CELL_CENTERS(1:N))

 ! calculate cell centers
 do I=1,n
   CELL_CENTERS(I) =  LB+I*H - H/2
 enddo

 ! initialise the vector of conserved quantities
  do I=1,N
    X = CELL_CENTERS(I)
    CALL  U0(x,CSV(1,:,i))
  enddo
 
  ! calculate primitive variables from conservative
  Rkstage = 1
  do I=1,N
   PV(Rkstage,1,i) = CSV(Rkstage,1,I)
   PV(Rkstage,2,I) = CSV(Rkstage,2,I)/CSV(Rkstage,1,I)
   PV(Rkstage,3,I) = (GM-1)*( CSV(Rkstage,3,I) - 0.5*CSV(Rkstage,2,I)*PV(Rkstage,2,I))
   PV(Rkstage,4,I) = sqrt(GM*PV(Rkstage,3,I)/PV(Rkstage,1,I))
  enddo
  
  ! set flow time to zero
   T_=0.D0
   
  OPEN(UNIT = 111, FILE = 'movie.dat', STATUS = 'UNKNOWN')
  WRITE(111,*)'TITLE="Solution" '
  WRITE(111,*)'VARIABLES="X" "rho" "u" "p"'
  Call OutputTecplot 
   
 End subroutine


 !%%%%%%%%%%%%%% Set up boundary conditions for given stage of the Runge Kutta marching %%%%%%%%%%
 
 Subroutine SetBC(Rkstage)
   Integer k,i,Rkstage
   
   ! set up ghost cells 
   
   Do i=-5,0
    do k=1,3
    CSV(Rkstage,k,i)  = CSV(Rkstage,k,abs(i)+1)
	enddo
    do k=1,4
     PV(Rkstage,k,i)  = PV(Rkstage,k,abs(i)+1)
	enddo
   Enddo

   Do i=1,6
    do k=1,3
     CSV(Rkstage,k,N+i)  = CSV(Rkstage,k,N-1-I)    
	enddo
    do k=1,4
     PV(Rkstage,k,N+i)  = PV(Rkstage,k,N-1-I)    
	enddo
   Enddo
 End subroutine   


 !%%%%%%%%%%% Compute initial data at t=0 for given spatial position 'x' %%%%%%%%%%%%%
 Subroutine U0(X,Q)
  ! U1 = RHO, U2 = RHOU, U3 = E
   Real X,U1,U2,U3,Q(3)
   Real DL,DR,UL,UR,PL,PR,x0
   Real :: pi= 3.141592653589793
   ! IvType = 1 : Sod' Shock Tube  Problem
   ! IvType = 2    Shock - turbulence interaction

  SELECT Case(IVTYPE)
   Case(1) 
    DL=1.0 ;   UL=0.0 ;  PL=1.0 
    DR=0.125 ;  UR=0.0 ;  PR=0.1 

    If (X .LE. 0.4)  THEN
     U1 = DL
     U2 = DL*UL
     U3 = PL/(GM-1) + 0.5*DL*UL**2
    Else
     U1 = DR
     U2 = DR*UR
     U3 = PR/(GM-1) + 0.5*DR*UR**2
    Endif
    
 Case(2)    
  ! Long time shock/turbulence interaction
  ! Mach number 1.1, S=1.5
  DL=   1.51569506726457     
  UL=   0.523345519274197     
  PL=   1.80500000000000     
  
  IF (X .LE. -4.5)  THEN
   U1 = DL
   U2 = DL*UL
   U3 = PL/(GM-1) + 0.5*DL*UL**2
  ELSE
   U1 = 1 + 0.1d0*sin(20*pi*x)
   U2 = 0.
   U3 = 1./(GM-1)  ! U = 0
  ENDIF
   
 End Select

 Q(1) = U1
 Q(2) = U2
 Q(3) = U3
end subroutine


!%%%%%%%%%%%%%%%% Time marching algorithm, which uses third order TVD RK method  %%%%%%%
!%%%%%%  Jiang G.S. and Shu C.W. Efficient Implementation of  weighted ENO schemes //J. Comput. Phys. 1996.  V. 126.  pp.202-212.

  Subroutine  ThirdOrderTVD
   Integer i,k

  ! loop stages from 1 to 3
   do Rkstage=1,3
     ! set up boundary conditions
     Call SetBc(Rkstage)
     ! calculate intercell fluxes
     Call ComputeFlux(Rkstage)
     ! perform the update
     CALL Update(Rkstage)
     Do i=1,n
       PV(Rkstage+1,1,i) = CSV(Rkstage+1,1,i)
       PV(Rkstage+1,2,i) = CSV(Rkstage+1,2,i)/CSV(Rkstage+1,1,i)
       PV(Rkstage+1,3,i) = (gm-1)*( CSV(Rkstage+1,3,i) - 0.5*PV(Rkstage+1,1,i)*PV(Rkstage+1,2,i)**2)
       PV(Rkstage+1,4,I) = sqrt(GM*PV(Rkstage+1,3,I)/PV(Rkstage+1,1,I))
     Enddo
   enddo
  
   ! re-assign the flow variables to stage 1 of RK method
   Do i=1,n
   do k=1,3
     CSV(1,k,i) = CSV(4,k,i)
   enddo
   do k=1,4
     PV(1,k,i) = PV(4,k,i)
   enddo
   Enddo

  End subroutine

 !%%%%%%%%%%% Solution update for each stage of TVD RK method %%%%%%%%%%
  Subroutine UPDATE(Rkstage)
   Integer I,K,Rkstage

   SELECT Case(Rkstage)
    Case(1)
     do i=1,n
      do K=1,3
       CSV(2,K,i)  =  CSV(1,K,i)  - (Dt/H)*(InterCellFlux(1,K,i) - InterCellFlux(1,K,i-1))
      enddo
     enddo

    Case(2)
     do i=1,n
      do K=1,3
       CSV(3,k,i) =  0.75*Csv(1,k,i) + 0.25*CSV(2,k,i)     - (0.25*Dt/H)*(InterCellFlux(2,k,i) - InterCellFlux(2,k,i-1))
      enddo
     enddo

    Case(3)
     do i=1,n
      do k=1,3
       CSV(4,K,i) =  (1./3)*CSV(1,K,i) + (2./3)*CSV(3,K,i)     - (2./3)*(Dt/H)*(InterCellFlux(3,K,i) - InterCellFlux(3,K,i-1))
	  enddo
     enddo
    END SELECT
  end subroutine


 !%%%%%%%%% write the output file %%%%%%%%%
  Subroutine    Output
  Integer i
   202 format(6(2x,e11.4))
   open(1,file='2ndORDER-HLL400.dat')
   WRITE(1,*)'TITLE="Solution" '
  WRITE(1,*)'VARIABLES="X" "rho""u""P" "T"' ! "u" "p"
   WRITE(1,*)'ZONE ',',I=',n, ',F="POINT"'
   do i=1,n   
    write(1,202) cell_centers(i),CSV(1,1,i),PV(1,2,I),PV(1,3,I),PV(1,3,I)/PV(1,1,I)
   enddo
   close(1)
  End subroutine   


  !%%%%%%%%%% Evaluation of the physical flux function from the conserved vector CDS =(rho,rho*u,E)
  Subroutine FluEval(CDS,Flux)
    Real cds(3),p,u,flux(3)
    
    u = cds(2)/cds(1)
    p = (gM-1)*(Cds(3) - 0.5*cds(1)*u**2)
    Flux(1) = cds(2)
    Flux(2) = cds(2)*u + p
    Flux(3) = (cds(3)+p)*u 
    
   End subroutine


  !%%%%%%%%%% Calculation of a stable time step %%%%%
  Subroutine ComputeTimeStep(dt)
   Integer i
   Real Umax,dt, a
  
   umax  = 0.0 
   Do i=1,n
    ! compute the sound speed
    a = ComputeSoundSpeed(CSV(1,:,i))
    umax = max(umax, a + abs(PV(1,2,i)))
   Enddo

   ! reduce the time step for first 10 time steps
   If ( IT<10) THEN
    dt = MIN(0.1*H/UMAX, TIME-T_)
   Else
    dt = MIN(CFL*H/UMAX, TIME-T_)
   Endif
  End subroutine

 
 !%%%%%%%%%%%%%%%% Calculation of the sound speed  on the conserved vector CDS
  Real function ComputeSoundSpeed(cds)
    Real cds(3),p,u  
    u = cds(2)/cds(1)
    p = (gm-1)*(cds(3) - 0.5*cds(2)*u)
    ComputeSoundSpeed=sqrt(gm*p/cds(1))  
  End function


  !%%%%%%%%%%%%%%%%% minmod slope limiter %%%%%%%%%%%%%%%%%
  Real function minmod(x,y)
   Real, intent (in) :: x,y
   minmod = 0.5*(sign(1.0,x)+sign(1.0,y))*min(abs(x),abs(y))
  End function minmod

 ! Compute the numerical flux
 Subroutine ComputeFlux(Rkstage)
  Integer i,k,Rkstage
  Real CDL(3),CDR(3),LocalFlux(3) 

  ! Loop over the spatial index i
  do I=0,N
   
   ! call reconstruction procedure at RK stage Rkstage to compute left CDL and right CDR
   ! values of the conserved vector between cells i and i+1
   CALL Reconstruction(CSV(Rkstage,:,i-2:i+3),CDL,CDR)
  
   ! calculate the numerical flux using reconstructed conserved vectors CDL, CDR
   ! Left   initial data  for the local Riemann problem is given by CDL
   ! Right  initial data  for the local Riemann problem is given by CDR

   Select Case(FluxType)
   Case(1)
      CALL LxF(CDL,CDR,LocalFlux)
      InterCellFlux(Rkstage,:,i) = LocalFlux

   Case(2)
      CALL Rusanov(CDL,CDR,LocalFlux)
      InterCellFlux(Rkstage,:,i) = LocalFlux

   Case(3)
      CALL HLL(CDL,CDR,LocalFlux)
      InterCellFlux(Rkstage,:,i) = LocalFlux

   Case(4)
      CALL HLLC(CDL,CDR,LocalFlux)
      InterCellFlux(Rkstage,:,i) = LocalFlux

   Case(5)
      CALL RoeSolver(CDL,CDR,LocalFlux)
      InterCellFlux(Rkstage,:,i) = LocalFlux
      
   Case default
    print*,' the flux is not defined. stop the program'
	read*
	stop
  End select	 

  Enddo
 end subroutine


  !%%%%%%%%% Lax Friedrich flux %%%%%%%%
   Subroutine LxF(CDL,CDR,Flux)
      Real  FL(3), FR(3),CDL(3),CDR(3),Flux(3)
      
        CALL FLUEVAL(CDL,FL)
        CALL FLUEVAL(CDR,FR)
        
        Flux = 0.5*(FL+FR) - 0.5*(h/dt)*(CDR-CDL) 
  End subroutine
  
  
 !%%%%%%%% Rusanov flux %%%%%%%%%%%%%%%%%%
   Subroutine Rusanov(CDL,CDR,Flux)
      Real  FL(3), FR(3),CDL(3),CDR(3),Flux(3),Speed, sa1, sa2
      CALL FLUEVAL(CDL,FL)
      CALL FLUEVAL(CDR,FR)
      sa1=ComputeSoundSpeed(CDL)
      sa2=ComputeSoundSpeed(CDR)
      speed=max(abs(CDL(2)/CDL(1))+sa1,abs(CDR(2)/CDR(1))+sa2)
      Flux = 0.5*(FL+FR) - 0.5*speed*(CDR-CDL)
  End subroutine 

 !%%%%%%%%%%%%%% HLL flux %%%%%%%%%%%%%%%%%%%%%
   Subroutine HLL(CDL,CDR,Flux)
      Real  FL(3), FR(3),CDL(3),CDR(3),Flux(3),sa1,sa2,WSL,WSR
      CALL FLUEVAL(CDL,FL)
      CALL FLUEVAL(CDR,FR)
      sa1=ComputeSoundSpeed(CDL)
      sa2=ComputeSoundSpeed(CDR)
      WSL=(CDL(2)/CDL(1))-sa1
      WSR=(CDR(2)/CDR(1))+sa2
      !WSL=min((),())
      !WSR=min((),())
      if (WSL >= 0.0) then
	    Flux = FL
	  else if (WSL <= 0.0 .and. WSR >= 0.0) then
	    Flux = (WSR*FL-WSL*FR+(WSL*WSR*(CDR-CDL)))/(WSR-WSL)
	  else if (WSR <= 0.0) then
	    Flux = FR
      end if  
  End subroutine 


 !%%%%%%%%%%%% HLLC flux %%%%%%%%%%%%%%%%%%%%%
   Subroutine HLLC(CDL,CDR,Flux)
      Real  FL(3), FR(3),CDL(3),CDR(3),Flux(3),sa1,sa2,WSL,WSR,WSS,CSL(3),CSR(3),PL,PR,UL,UR,RL,RR,EL,ER
      CALL FLUEVAL(CDL,FL)
      CALL FLUEVAL(CDR,FR)

      ul = cdl(2)/cdl(1)
      el=cdl(3)
      rl=cdl(1)
      pl = (gM-1)*(el - 0.5*rl*ul**2)

      ur = cdr(2)/cdr(1)
      er=cdr(3)
      rr=cdr(1)
      pr = (gM-1)*(er - 0.5*rr*ur**2)

      sa1=ComputeSoundSpeed(CDL)
      sa2=ComputeSoundSpeed(CDR)
      !CSL=rl*((WSL-CDL)/(WSL-WSS))
      !CSR=rr*((WSR-CDR)/(WSR-WSS))
      WSL=(CDL(2)/CDL(1))-sa1
      WSR=(CDR(2)/CDR(1))+sa2
      WSS=(pr-pl+(rl*UL*(WSL-UL))-(rr*UR*(WSR-UR)))/((rl*(WSL-UL))-(rr*(WSR-UR)))

	CSL(1)= rl*((WSL-UL)/(WSL-WSS))*1
	CSL(2)= rl*((WSL-UL)/(WSL-WSS))*WSS
	CSL(3)= rl*((WSL-UL)/(WSL-WSS))*((EL/RL)+(WSS-UL)*(WSS+(PL/(RL*(WSL-UL)))))

	CSR(1)= rr*((WSR-UR)/(WSR-WSS))*1
	CSR(2)= rr*((WSR-UR)/(WSR-WSS))*WSS
	CSR(3)= rr*((WSR-UR)/(WSR-WSS))*((ER/RR)+(WSS-UR)*(WSS+(PR/(RR*(WSR-UR)))))
       
	if (WSL >= 0.0) then
	    Flux = FL
	  else if (WSL <= 0.0 .and. WSS >= 0.0) then
	    Flux = FL+(WSL*(CSL-CDL))
	  else if (WSS <= 0.0 .and. WSR >= 0.0) then
	    Flux = FR+(WSR*(CSR-CDR))
	  else if (WSR <= 0.0) then
	    Flux = FR
	end if
  End subroutine 

 !%%%%%%%%%%%% Roe solver %%%%%%%%%%%%%%%%%%%%%
  Subroutine RoeSolver(CDL,CDR,Flux)
      Real  FL(3), FR(3),CDL(3),CDR(3),Flux(3),eigenvectors(3,3)
      Real  A1,A2,A3,L1,L2,L3, R1(3),R2(3),R3(3)
      Real  alpha1,alpha2,alpha3,deltaq1,deltaq2,deltaq3
      Real  eigenvalues1,eigenvalues2,eigenvalues3
      Real  RoeAverage1,roeaverage2,roeaverage3, u, a, H
      ! Declare variables
	  real(kind=8) :: sqrtrhoL, sqrtrhoR
	  real(kind=8) :: HroeL, HroeR, sa1, sa2,PL,PR,UL,UR,RL,RR,EL,ER

      
      !loop index
    
      integer :: k

      CALL FLUEVAL(CDL,FL)
      CALL FLUEVAL(CDR,FR)

   !%%%%%%%%compute roe averages%%%%%%%%!
	   ! Calculate conserved variables for left and right states
	  ul = cdl(2)/cdl(1)
	  el=cdl(3)
	  rl=cdl(1)
	  pl = (gM-1)*(el - 0.5*rl*ul**2)

	  ur= cdr(2)/cdr(1)
	  er=cdr(3)
	  rr= cdr(1)
	  pr = (gM-1)*(er - 0.5*rr*ur**2)

	
	  ! Calculate enthalpies for left and right states
	  HroeL = (el+pl)/rl
	  sa1 = ComputeSoundSpeed(CDL)
	  HroeR = (er+pr)/rr
	  sa2 = ComputeSoundSpeed(CDR)
	  
	   ! Calculate square root of density for left and right states
	  sqrtrhoL = sqrt(rl)
	  sqrtrhoR = sqrt(rr)

	  ! Compute Roe-averaged values
	  RoeAverage1 = (sqrtrhoL*uL + sqrtrhoR*uR)/(sqrtrhoL+sqrtrhoR)
	  
	  RoeAverage2 = (sqrtrhoL*HroeL + sqrtrhoR*HroeR)/(sqrtrhoL+sqrtrhoR)

	  RoeAverage3 = sqrt((gM-1)*(RoeAverage2-0.5*RoeAverage1**2))


    !%%%%%%%%compute eigenvalues%%%%%%%%!
	  eigenvalues1 = RoeAverage1 - RoeAverage3
	  eigenvalues2 = RoeAverage1
	  eigenvalues3 = RoeAverage1 + RoeAverage3

    !%%%%%%%%compute right-eigenvectors%%%%%%%%!
	! Extract Roe-averaged values
	  u = RoeAverage1
	  H = Roeaverage2
	  a = RoeAverage3

	  ! Compute averaged right eigenvectors
	  eigenvectors(1, 1) = 1.0
	  eigenvectors(2, 1) = u - a
	  eigenvectors(3, 1) = H -(u* a)

	  eigenvectors(1, 2) = 1.0
	  eigenvectors(2, 2) = u
	  eigenvectors(3, 2) = 0.5*u*u

	  eigenvectors(1, 3) = 1.0
	  eigenvectors(2, 3) = u + a
	  eigenvectors(3, 3) = H + (u*a)

    !%%%%%%%%compute right-eigenvectors%%%%%%%%!
	 ! Compute deltaQ vector (jump in conserved variables)
	  deltaQ1 = CDR(1) - CDL(1)
	  deltaQ2 = CDR(2) - CDL(2)
	  deltaQ3 = CDR(3) - CDL(3)

	   
	  !compute wave strengths 
	  alpha2 = (gM-1)/(a*a)*(deltaQ1*(H-u**2)+u*deltaQ2-deltaQ3)
	  alpha1 = 1/(2*a)*(deltaQ1*(u+a)-deltaQ2-a*alpha2)
	  alpha3 = deltaQ1-(alpha1+alpha2)

   !%%%%%%%% assemble scalars and matrices for flux calculation %%%%%%%%!
      A1=alpha1
      A2= alpha2
      A3=alpha3

      L1= abs(eigenvalues1)
      L2= abs(eigenvalues2)
      L3= abs(eigenvalues3)

      R1=eigenvectors(:,1)
      R2=eigenvectors(:,2)
      R3=eigenvectors(:,3)

      Flux= 0.5 * (FL + FR)-0.5*((A1*L1*R1)+(A2*L2*R2)+(A3*L3*R3))
      
  End Subroutine RoeSolver

 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 !%%%%%%%%%%%%%%%% Reconstruction procedure%%%%%%%%%%%%%%
 ! Input: one-dimensional array U1D of flow quantities near cell interface i+1/2
 ! Output: left CDL and right CDR values at interface
 Subroutine Reconstruction(U1D,CDL,CDR)
   Integer I, F
   Real U1d(3,-2:3),CDL(3),CDR(3)
   Real r, theta
  Real, Parameter :: epsilon = 1.0e-10  ! singularity prevention measure
   
   select Case(SpatialOrder)
  
	 Case(1)
	 
	 ! First order 
	  Do f=1,3
	   CDL(f) = U1D(f,0)  
	   CDR(f) = U1D(f,1)  	
	  Enddo 
	 
	 ! second order TVD 
	 Case(2)
         Do f = 1, 3
         r = (U1D(f, 1) - U1D(f, 0)) / (U1D(f, 2) - U1D(f, 1) + epsilon)
        theta = minmod(1.0, r)
        !theta =max(0,min(1.0,r))
        CDL(f) = U1D(f, 0) + 0.5 * theta * (U1D(f, 1) - U1D(f, 0))

        r = (U1D(f, 2) - U1D(f, 1)) / (U1D(f, 1) - U1D(f, 0) + epsilon)
        theta = minmod(1.0, r)
        CDR(f) = U1D(f, 1) - 0.5 * theta * (U1D(f, 2) - U1D(f, 1))
         Enddo
        !print*,'no reconstruction found. stop the code!'
        !stop  

	 Case default
	 print*,' Wrong spatial accuracy. Stop the code'
	 stop  
  end select

 end subroutine

 Subroutine OutPutTecplot
    Integer i,j
	real(8) x
    
   55 Format (4(2x,e11.4))
    WRITE(111,*)'ZONE ',',I=',n, ',F="POINT"'
    WRITE(111,*) ', SOLUTIONTIME=',T_
    DO  I = 1,n
	  x = lb + i*h-h/2
      WRITE(111,55)X,pv(1,1,i),pv(1,2,i),pv(1,3,i)
    Enddo

 End subroutine

 END 
