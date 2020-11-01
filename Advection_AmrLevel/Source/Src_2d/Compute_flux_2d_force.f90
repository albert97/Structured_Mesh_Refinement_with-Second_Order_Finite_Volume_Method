

module compute_flux_module

  implicit none

  private :: FORCE_FLUX, HLLC_FLUX, MUSCL_TVD_FLUX, &
             DATA_RECONSTRUCTION, SLOPE_LIMITER

  public :: compute_flux_2d

contains

  subroutine compute_flux_2d(DIRPARAM, lo, hi, dt, dx, &
                             phi,ph_lo,ph_hi, &
                             flxx, fx_lo, fx_hi, &
                             flxy, fy_lo, fy_hi)

    use global_var_module
    
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
    
    !ULL AND URR FOR SECOND ORDER METHODS 
    !NUM_FLUX IS TO STORE THE FLUXES CALCUALTED FORMEACH FLUX SCHEMES 
    double precision :: UL( nVarQ ), UR( nVarQ ), ULL( nVarQ ), &
                        URR( nVarQ ), NUM_FLUX( nVarQ )      
    integer :: i, j
    
    
    !FluxType DETERMINES THE METHOD USED TO DETERMINE THE FLUX 
    integer :: FluxType = 2
    
    
    flxx(:,:,:) = 0.d0 
    flxy(:,:,:)= 0.d0 
    
    write(*,*)" DIRPARAM",  DIRPARAM	
    
    
    if ( DIRPARAM == 0 ) then 
    
    !--------------------- FLUX IN THE X DIR ---------------------------
    	write(*,*)"X sweep"
    	write(*,*)"Density",phi( : , : ,1)
		do j = lo(2),hi(2)
			do i = lo(1), hi(1) +1
			    
			    ULL = phi(i-2,j,1:nVarQ)
			    
!			    WRITE (*,*) "ULL",   phi(i-2,j,1:nVarQ)
				UL  = phi(i-1,j,1:nVarQ)
			    UR  = phi(i,j,1:nVarQ)
			    URR = phi(i+1,j,1:nVarQ)		    
			    
			    
			    if ( FluxType == 0 ) then 
			    
					call FORCE_FLUX( DIRPARAM, UL, UR, dx( 1 ), dt, &
			                         NUM_FLUX ) 
			    
			    else if ( FluxType == 1 ) then 	
			    
					call HLLC_FLUX ( DIRPARAM, UL, UR, NUM_FLUX ) 
					
			    else if ( FluxType == 2 ) then
			    
					call MUSCL_TVD_FLUX ( DIRPARAM, ULL, UL, UR, URR, &
					                    & dt, dx( 1 ), NUM_FLUX ) 
					                    		    
			    end if		         
				
				flxx(i,j,1:nVarQ) = NUM_FLUX( : )
				
				! WRITE (*,*) "flxx(i,j,1:nVarQ)", flxx(i,j,1:nVarQ)
			end do
			
		end do
	
	else if ( DIRPARAM == 1 ) then 		
	
	!--------------------- FLUX IN THE Y DIR ---------------------------
	
	write(*,*)"Y sweep"
			write(*,*)"Density",phi( : , : ,1)
		do j = lo(2),hi(2) +1
			do i = lo(1), hi(1) 
		
		        ULL = phi(i,j-2,1:nVarQ)
				UL  = phi(i,j-1,1:nVarQ)
			    UR  = phi(i,j,1:nVarQ)
			    URR = phi(i,j+1,1:nVarQ)
			    
			    if ( FluxType == 0 ) then 
			    
					call FORCE_FLUX( DIRPARAM, UL, UR, dx( 2 ), dt, &
			                         NUM_FLUX ) 
			    
			    else if ( FluxType == 1 ) then 	
			    
					call HLLC_FLUX ( DIRPARAM, UL, UR, NUM_FLUX )
					
				else if ( FluxType == 2 ) then
			    
					call MUSCL_TVD_FLUX ( DIRPARAM, ULL, UL, UR, URR, &
					                    & dt, dx( 2 ), NUM_FLUX )  
			    
			    end if	 			         
				
				flxy(i,j,1:nVarQ) = NUM_FLUX( : )
						
			end do
		end do
			
	end if	
	   

  end subroutine compute_flux_2d
  
  
  
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  
  !----------------------    FORCE FLUX    -----------------------------
  
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  
  
  
  subroutine FORCE_FLUX ( DIRPARAM, UL, UR, dxy, dt, FORCE_NF ) 
  

	use global_var_module
	use auxiliary_module
	
	implicit none
	
	integer, intent( in ) :: DIRPARAM	
	double precision, intent ( in ) :: UL( nVarQ ), UR( nVarQ ), dxy, dt	
	double precision,intent (inout ) :: FORCE_NF ( nVarQ )  
	
	double precision :: LFFlux( nVarQ ), RIFlux( nVarQ ), FL( nVarQ ), &
	                    FR( nVarQ ), QInter( nVarQ )   
	
	
	call Lax_Friedrich_flux(DIRPARAM, UL, UR, LFFlux, dxy, dt)

	call Richtmyer_flux(DIRPARAM, UL, UR, RIFlux, dxy, dt)

    FORCE_NF = 0.5d0*( LFFlux( : ) + RIFlux( : ) )
    	
	end subroutine FORCE_FLUX



  !--------------------------------------------------------------------!
  !--------------------------------------------------------------------!
   
  !------------------------ Lax_Friedrich FLUX  -----------------------!
 
  !--------------------------------------------------------------------!
  !--------------------------------------------------------------------!
  
  subroutine Lax_Friedrich_flux(DIRPARAM, UL, UR, LFFlux, dxy, dt)
  
  	use global_var_module
	use auxiliary_module
	
	implicit none
	
	integer, intent( in ) :: DIRPARAM	
	double precision, intent ( in ) :: UL( nVarQ ), UR( nVarQ ), dxy, dt	
	double precision, intent( out ) :: LFFlux( nVArQ ) 
	
	call FLUX_FUNCT(DIRPARAM, (/1,1/), (/1,1/), UL, FL) 
    call FLUX_FUNCT(DIRPARAM, (/1,1/), (/1,1/), UR, FR) 
				    
    LFFlux( : ) = 0.5d0*( FL( : ) + FR( : ) ) - dxy/dt* &
		        & 0.5d0*( UR( : ) - UL( : ) ) 

  end subroutine Lax_Friedrich_flux
  
  
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  
  !----------------------  RICHTMYER FLUX   -------------------------
  
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  
  
  subroutine Richtmyer_flux(DIRPARAM, UL, UR,  RIFlux, dxy, dt)
	
  	use global_var_module
	use auxiliary_module
	
	implicit none
	
	integer, intent( in ) :: DIRPARAM	
	double precision, intent ( in ) :: UL( nVarQ ), UR( nVarQ ), dxy, dt	
	double precision, intent( out ) ::  RIFlux( nVArQ ) 
	double precision:: QInter
	
	call FLUX_FUNCT(DIRPARAM, (/1,1/), (/1,1/), UL, FL) 
    call FLUX_FUNCT(DIRPARAM, (/1,1/), (/1,1/), UR, FR) 
				    
	QInter( : ) = 0.5d0*( UL( : ) + UR( : ) ) - dt/dxy* &
		        & 0.5d0*( FR( : ) - FL( : ) )    
		         
    call FLUX_FUNCT(DIRPARAM, (/1,1/), (/1,1/), QInter, RIFlux)  
  
  
  end subroutine
  
end module compute_flux_module

