


module global_var_module

	implicit none

	double precision :: Gamma = 1.4d0
	integer :: nVar = 8
	integer :: nVarQ = 4
	double precision :: tStop = 0.25d0
	integer :: No_cell = 128

end module global_var_module

module global_function_module

  implicit none
   
  private

  public :: P2C, C2P, FLUX_FUNCT

contains

  subroutine C2P ( Q_lo, Q_hi, Q, W ) 
 
	use global_var_module
	implicit none
	
	
	integer, intent (in) :: Q_lo(2), Q_hi(2) 
    double precision, intent(in)  :: Q (Q_lo(1):Q_hi(1),Q_lo(2):Q_hi(2), nVarQ)
    double precision, intent(out) :: W (Q_lo(1):Q_hi(1),Q_lo(2):Q_hi(2), nVarQ)
			
	W (:,:, 1 ) = Q (:,:, 1 ) 
	W (:,:, 2 ) = Q (:,:, 2 )/Q (:,:, 1 ) 
	W (:,:, 3 ) = Q (:,:, 3 )/Q (:,:, 1 ) 
	W (:,:, 4 ) = ( Q (:,:, 4 ) - 0.5*( Q ( :,:, 2 ) **2.d0 / Q ( :,:, 1 ) + &
	            &   Q (:,:,3 )** 2.d0 / Q (:,:, 1 ) ) )*( Gamma - 1 )	   

  end subroutine C2P
  
  
  !---------------------------------------------------------------------
  
  subroutine P2C ( W_lo, W_hi, W, Q ) 
  
	
	use global_var_module
	implicit none
		
	integer, intent (in) :: W_lo(2), W_hi(2) 
	double precision, intent(in)  :: W (W_lo(1):W_hi(1),W_lo(2):W_hi(2), nVarQ)
	double precision, intent(out) :: Q (W_lo(1):W_hi(1),W_lo(2):W_hi(2), nVarQ)

	Q (:,:, 1 ) = W (:,:, 1 ) 
	Q (:,:, 2 ) = W (:,:, 1 )*W (:,:, 2 ) 
	Q (:,:, 3 ) = W (:,:, 1 )*W (:,:, 3 ) 
	Q (:,:, 4 ) = W (:,:, 4 )/( Gamma - 1.d0 ) + 0.5d0*W (:,:, 1 )*( &
	            & W (:,:, 2 )**2.d0 + W (:,:, 3 )**2.d0 )     
	      
  end subroutine P2C
  
  
  !---------------------------------------------------------------------
  
  
  subroutine FLUX_FUNCT(DIRPARAM, Q_lo, Q_hi, Q, F)
	
	use global_var_module
	implicit none
	
	integer, intent (in) :: DIRPARAM, Q_lo(2), Q_hi(2) 
    double precision, intent(in)  :: Q (Q_lo(1):Q_hi(1),Q_lo(2):Q_hi(2), nVarQ)
	double precision, intent(out)  :: F (Q_lo(1):Q_hi(1),Q_lo(2):Q_hi(2), nVarQ)
	double precision ::  P (Q_lo(1):Q_hi(1),Q_lo(2):Q_hi(2))
	
	P (:,:) = ( Q (:,:, 4 ) - 0.5*( Q ( :,:, 2 ) **2.d0 / Q ( :,:, 1 ) + &
	        &   Q (:,:, 3 )** 2.d0 / Q (:,:, 1 ) ) )*( Gamma - 1 )
	
	if ( DIRPARAM == 0 ) then             
	
		F (:,:,1)  = Q (:,:,2) 	
		F (:,:,2)  = Q (:,:,2) ** 2.d0 / Q (:,:,1) + P(:,:) 		
		F (:,:,3)  = Q (:,:,2) *Q (:,:,3) / Q (:,:,1) 	
		F (:,:,4)  = Q (:,:,2) / Q (:,:,1) * ( Q (:,:,4) + P(:,:) ) 
		
	else if ( DIRPARAM == 1 ) then             
	
		F (:,:,1)  = Q (:,:,3) 	
		F (:,:,2)  = Q (:,:,2) *Q (:,:,3) / Q (:,:,1) 	
		F (:,:,3)  = Q (:,:,3) ** 2.d0 / Q (:,:,1) + P(:,:) 		
		F (:,:,4)  = Q (:,:,3) / Q (:,:,1) * ( Q (:,:,4) + P(:,:) ) 
		
	end if

	
  end subroutine FLUX_FUNCT
  
  
  !---------------------------------------------------------------------
  
  
  

end module auxiliary_module
