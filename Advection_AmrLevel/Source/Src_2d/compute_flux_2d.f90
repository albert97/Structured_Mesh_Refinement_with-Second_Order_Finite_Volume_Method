
module compute_flux_module

  implicit none

  private :: HLLC_FLUX, MSCLE_HANCOCK, &
             DATA_RECONSTRUCTION, SLOPE_LIMITER

  public :: compute_flux_2d

contains

  subroutine compute_flux_2d(DIRPARAM, lo, hi, dt, dx, &
                             phi,ph_lo,ph_hi, &
                             flxx, fx_lo, fx_hi, &
                             flxy, fy_lo, fy_hi)

    use global_var_module
    use global_function_module
    
	implicit none
    
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: ph_lo(2), ph_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    
    integer, intent(in) :: DIRPARAM	
    
    double precision, intent(in) :: dt, dx(2)
    double precision, intent(in) :: phi (ph_lo(1):ph_hi(1), &
                                         ph_lo(2):ph_hi(2),nVar)    

    double precision, intent(out) :: flxx(fx_lo(1):fx_hi(1), &
                                          fx_lo(2):fx_hi(2),nVar)
    double precision, intent(out) :: flxy(fy_lo(1):fy_hi(1), &
                                          fy_lo(2):fy_hi(2),nVar)
    
    double precision :: QL( nVarQ ), QR( nVarQ ), QLL( nVarQ ), &
                        QRR( nVarQ ), NUM_FLUX( nVarQ )      
    integer :: i, j
    
    
    flxx(:,:,:) = 0.d0 
    flxy(:,:,:)= 0.d0 
    

    
    if ( DIRPARAM == 0 ) then 

    !--------------------- FLUX IN THE X DIR ---------------------------
    	
		do j = lo(2),hi(2)
			do i = lo(1), hi(1) +1
			    
			    QLL = phi(i-2,j,1:nVarQ)
				QL  = phi(i-1,j,1:nVarQ)
			    QR  = phi(i,j,1:nVarQ)
			    QRR = phi(i+1,j,1:nVarQ)		    
		
			!	call HLLC_FLUX ( DIRPARAM, QL, QR, NUM_FLUX ) 
	            		   
				call MSCLE_HANCOCK( DIRPARAM, QLL, QL, QR, QRR, &
					                    & dt, dx( 2 ), NUM_FLUX )  
					                    
					                    
				flxx(i,j,1:nVarQ) = NUM_FLUX( : )
				
			end do
			
		end do
	
	else if ( DIRPARAM == 1 ) then 		
	
!	write(*,*)"Density",phi(:,:,1)
!	write(*,*)"ph(32,0)",phi(32,0,  1)
	!--------------------- FLUX IN THE Y DIR ---------------------------

		do j = lo(2),hi(2) +1
		  
			do i = lo(1), hi(1) 

		        QLL = phi(i,j-2,1:nVarQ)
				QL  = phi(i,j-1,1:nVarQ)
			    QR  = phi(i,j,1:nVarQ)
			    QRR = phi(i,j+1,1:nVarQ)
	!	    	 write(*,*)"QLL", QLL, "QL",QL ,"QR", QR, "QRR", QRR, i ,j
  
  
			!	call HLLC_FLUX ( DIRPARAM, QL, QR, NUM_FLUX ) 
			    
					call MSCLE_HANCOCK( DIRPARAM, QLL, QL, QR, QRR, &
					                    & dt, dx( 2 ), NUM_FLUX )  
		
				flxy(i,j,1:nVarQ) = NUM_FLUX( : )
						
			end do
		end do
			
	end if	
	   

  end subroutine compute_flux_2d
  
  
   !---------------------    WAVESPEED_ESTIMATE  ----------------------!
  
  subroutine P_BASED_WAVESPEED( PL, PR, VNormL, VNormR, RHOL, RHOR, WAVESPEED)
  
  use global_var_module
  use auxiliary_module
  
  implicit none
  
  double precision, intent(in) :: PL, PR, VNormL, VNormR, RHOL, RHOR
  double precision, intent(inout) :: WAVESPEED(2)
  double precision :: AL, AR, A_BAR
  double precision :: RHO_BAR
  double precision :: PPVRS, PSTAR
  double precision :: Q_L, Q_R
  
  
  !	write(*,*) "PL/RHOL", PL/RHOL, "PL", PL , "RHOL",  RHOL
  	
  !	write(*,*) "PR/RHOR", PR/RHOR, "PL", PR , "RHOL",  RHOR
  	
  AL = DSQRT( Gamma * (PL/RHOL))
  AR = DSQRT(Gamma * (PR/RHOR))
  A_BAR = 0.5d0 * (AL+AR)
  RHO_BAR = 0.5d0 * (RHOL + RHOR)
  
  PPVRS = 0.5d0 * (PL + PR) - 0.5d0 * (VNormR - VNormL) * RHO_BAR * A_BAR
  PSTAR = MAX( 0.0, PPVRS )
   
    if ( PSTAR .LE. PL ) then 
    
		Q_L = 1 
		
	else 
	
	!write(*,*) "PSTAR / PL", PSTAR / PL, "PSTAR", PSTAR , "PL",  PL

		Q_L =  DSQRT( 1 + ( Gamma + 1 )/( 2*Gamma )* &
		    &  ( PSTAR / PL - 1 ) )    
    
    end if 
    
    if (PSTAR .LE. PR ) then 
    
		Q_R = 1 
		
	else 

	
	
		Q_R =  DSQRT( 1 + ( Gamma + 1 )/( 2*Gamma )* &
		    &  ( PSTAR/PR ) - 1 )   
    
    end if 
    
    WAVESPEED(1) = VNormL - AL * Q_L
    WAVESPEED(2) = VNormR + AR * Q_R    
              
  
  end subroutine P_BASED_WAVESPEED
  
  
  
  
     !---------------------  STAR_STATE  ----------------------!
     
     
   subroutine Q_STAR(DIRPARAM, SK, SSTAR, VNormK, PK, RHOK, VTanK, EK, QKSTAR)
   
   use global_var_module
   	use auxiliary_module
   implicit none 
   integer, intent(in) :: DIRPARAM
   double precision, intent(in) :: SK, SSTAR, VNormK, PK, RHOK, VTanK, EK
   double precision, intent(inout) :: QKSTAR(nVarQ)
   double precision :: ft
   
   ft = RHOK * ((SK - VNormK)/(SK - SSTAR))
   QKSTAR (1) = ft
   QKSTAR (4) = ft * ( EK / RHOK + ( SStar - VNormK )*(  &
                 & SStar + PK/( RHOK*( SK - VNormK) ) ) ) 
                 
  
    if ( DIRPARAM == 0 ) then 

		QKSTAR ( 2 ) = ft*SSTAR 
		QKSTAR ( 3 ) = ft*VTanK

    else 
 
		QKSTAR ( 2 ) = ft*VTanK
		QKSTAR ( 3 ) = ft*SSTAR 

    end if 
    
   end subroutine Q_STAR
   
  
  !----------------------    HLLC FLUX    -----------------------------!
  
  
  subroutine HLLC_FLUX( DIRPARAM, QL, QR, FLUX )

    use global_var_module
    use auxiliary_module
    
	implicit none
	integer, intent(in)::DIRPARAM
	double precision, intent(in)::QL(nVarQ)
	double precision, intent(in)::QR(nVarQ)
	double precision, intent(inout):: FLUX(nVarQ)
	double precision :: WL(nVarQ), WR(nVarQ)
	double precision :: PL, PR, RHOL, RHOR, VNormL, VNormR, VTanL, VTanR, EL, ER
	double precision :: SL, SR, SSTAR
	double precision :: FL(nVarQ), FR(nVArQ)
	double precision :: FSTARL(nVarQ), FSTARR(nVarQ)
	integer :: i
	double precision :: WAVESPEED(2)
	double precision ::  QLSTAR(nVarQ),  QRSTAR(nVarQ)
	!-----Wavespeed EStimate--------!
	
	
    call C2P( (/1,1/),(/1,1/),QL, WL ) ;
    call C2P( (/1,1/),(/1,1/),QR, WR ) ;
    
    RHOL = WL(1)
    RHOR = WR(1)
    PL = WL(4)
    PR = WR(4)
        
    if ( DIRPARAM == 0 ) then
    
		VNormL = WL( 2 )  
		VNormR = WR( 2 ) 
		VTanL  = WL( 3 )
		VTanR  = WR( 3 )    
    
    else 
    
		VNormL = WL( 3 )  
		VNormR = WR( 3 ) 
		VTanL  = WL( 2 )
		VTanR  = WR( 2 )    
    
    end if 

  
  
	call P_BASED_WAVESPEED(PL, PR, VNormL, VNormR, RHOL, RHOR, WAVESPEED)

	SL = WAVESPEED(1)
	SR = WAVESPEED(2)
	

	
   !----------------------------Star State-----------------------------!
    
    SSTAR= ( PR - PL + RHOL * VNormL*( SL - VNormL ) - &
          &   RHOR * VNormR*( SR - VNormR ) )/( RHOL * ( SL - &
          &   VNormL ) - RHOR*( SR - VNormR ) ) 
              

	    
    call FLUX_FUNCT(DIRPARAM, (/1,1/), (/1,1/), QL, FL)  
    call FLUX_FUNCT(DIRPARAM, (/1,1/), (/1,1/), QR, FR)  


 
    call Q_STAR( DIRPARAM, SL, SSTAR, VNormL, PL, RHOL, VTanL, QL(4), QLSTAR)
    call Q_STAR( DIRPARAM, SR, SSTAR, VNormR, PR, RHOR, VTanR, QR(4), QRSTAR)


    do i = 1, nVarQ
    
		FSTARL( i ) = FL( i ) + SL*( QLSTAR( i ) - QL( i ) ) 	
		FSTARR( i ) = FR( i ) + SR*( QRSTAR( i ) - QR( i ) ) 		
		    
    end do 
    

    
  !------------------------------ HLLC FLUX ---------------------------!
 

 
   if ( SL .GE. 0.d0 ) then 
    
		FLUX( : ) = FL( : ) 
    
    end if 
    
    if ( SL .LT. 0.d0 .AND. SSTAR .GE. 0.d0   ) then 
    
		FLUX( : ) = FSTARL( : ) 
    
    end if 
    
    if ( SSTAR .LT. 0.d0 .AND. SR .GE. 0.d0   ) then 
    
		FLUX( : ) = FSTARR( : ) 
    
    end if 
    
    if ( SR .LT. 0.d0 ) then 
    
		FLUX( : )  = FR( : ) 

    end if 	
	
  end subroutine  HLLC_FLUX
  
  
   
   !-----------------------------MUSCL HANCOCK-------------------------!
   
   subroutine MSCLE_HANCOCK( DIRPARAM,QLL,QL,QR,QRR,dt,dxy,FLUX)
   
	use global_var_module
	use auxiliary_module
	
	implicit none
	
	integer, intent( in ) :: DIRPARAM
	double precision, intent ( in ) :: dt, dxy 
	double precision, intent ( in ) :: QLL( nVarQ ), QL( nVarQ ), &
	                                   QR( nVarQ ), QRR( nVarQ )
	double precision,intent ( inout ) ::FLUX( nVarQ ) 
	
	double precision :: QUPDATEL(2, nVarQ), QUPDATER(2, nVarQ)
	double precision :: QExLL( nVarQ ), QExLR( nVarQ ), &
	                    QExRL( nVarQ ), QExRR( nVarQ ), &
	                    QEvL( nVarQ ), QEvR( nVarQ ) 
	                    
	double precision :: FLL( nVarQ ), FLR( nVarQ ), FRL( nVarQ ), & 
	                    FRR( nVarQ )
	                    
	integer :: i 
	
	

	call DATA_RECONSTRUCTION( QLL, QL, QR,  QUPDATEL)           
	
	QExLL ( : ) = QUPDATEL( 1, : )       
	QExLR ( : ) = QUPDATEL( 2, : )    
	

		
	call DATA_RECONSTRUCTION( QL, QR, QRR,  QUPDATER)           
	QExRL ( : ) = QUPDATER( 1, : )       
	QExRR ( : ) = QUPDATER( 2, : )    
		

		
	call FLUX_FUNCT(DIRPARAM, (/1,1/), (/1,1/), QExLL, FLL)
	call FLUX_FUNCT(DIRPARAM, (/1,1/), (/1,1/), QExLR, FLR)
	call FLUX_FUNCT(DIRPARAM, (/1,1/), (/1,1/), QExRL, FRL)
	call FLUX_FUNCT(DIRPARAM, (/1,1/), (/1,1/), QExRR, FRR)    
	

	
	QEvL( : ) = QExLR( : ) + 0.5d0*dt/dxy*( FLL( : ) - FLR( : ) ) 
	QEvR( : ) = QExRL( : ) + 0.5d0*dt/dxy*( FRL( : ) - FRR( : ) ) 
	

	
	call HLLC_FLUX( DIRPARAM, QEvL, QEvR, FLUX) 
	    	
   end subroutine  MSCLE_HANCOCK
   

  
 subroutine DATA_RECONSTRUCTION( QL, QC, QR, QUPDATE ) 
  
	use global_var_module
  
	implicit none
	
	double precision, intent ( in ) :: QL( nVarQ ), QC( nVarQ ), &
	                                   QR( nVarQ )  
	double precision,intent ( inout ) :: QUPDATE( 2, nVarQ )
	
	double precision :: di1, di2, r,  SlopeLimiter
	double precision :: SlopeLR, TOL
	double precision ::di,sdi, Omega, BetaCR
	integer :: i 
	Omega = 0.d0  
	BetaCR = 1.d0 
	TOL = 1E-10 
	
	
	do i = 1, nVarQ
	
		di1 = QC( i ) - QL( i ) 
		di2 = QR( i ) - QC( i ) 

		
		if ( ABS( di2 ) .LT. TOL ) then 
					
			di2 = TOL 
			
			
				
			if ( ABS( di1 ) .LT. TOL ) then 
				
				r = 1.d0 
				
			!	write(*,*)"di1", di1
				
			else 
				
				r = di1/di2 
				
				if ( r .NE. r .OR. r .GE. HUGE(TOL) ) then 
			
			       r = 0.d0 ;
			
				end if  
				
			end if 
			
		else
			
			r = di1/di2  			
			
		end if 		

		SlopeLR= 2.d0*BetaCR/(( 1.d0 - Omega )+ &
		          & ( 1.d0 + Omega )*r) 
		           		
		call SLOPE_LIMITER(r, SlopeLR, SlopeLimiter)

	    di = 0.5d0*( 1.d0 + Omega )*di1+ &
		      & 0.5d0*( 1.d0 - Omega )*di2		
		sdi = SlopeLimiter * di

	
		QUPDATE(1,  i ) =  QC(i) - 0.5 * sdi
		QUPDATE(2,  i ) =  QC(i) + 0.5 * sdi
			
			
	end do             
	                    
	             

end subroutine DATA_RECONSTRUCTION
  !----------------------   SLOPE LIMITER   ----------------------------
  
  
  subroutine SLOPE_LIMITER(r, SlopLR, SlopeLimiter ) 
  
	use global_var_module
  
	implicit none
	
	double precision, intent ( in ) :: r, SlopLR

	double precision,intent ( inout ) :: SlopeLimiter 
	

	if ( r .LT. 0.d0 ) then 
					
		SlopeLimiter  = 0.d0 
					
	end if 
				
	if ( r .GE. 0.d0 .AND. r .LT. 0.5d0 ) then 
					
		SlopeLimiter  = 2.d0*r
					
	end if 
				
	if ( r .GE. 0.5d0 .AND. r .LT. 1.d0 ) then 
					
		SlopeLimiter  = 1.d0 
					
	end if
				
	if ( r .GE. 1.d0 ) then 
					
		SlopeLimiter  = MIN( 2.d0, MIN( r, SlopLR) )
					
	end if
		
	     
  end subroutine SLOPE_LIMITER
  
end module compute_flux_module

