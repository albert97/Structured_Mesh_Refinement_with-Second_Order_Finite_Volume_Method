
!
!PHI IS THE 2D VECTOR WITH ALL CONSERVTIVE AND PRIMITIVE VARIABLES
!

!
!VX, VY ARE USED TO STORE THE VELOCITY IN X AND Y DIRECTION 
!

!
!THE VELOCITY OUTPUTS ARE THE SUM OF LOCAL VELOCITY AND SOUND SPEED
!

!
!TO FIND THE SOUND SPEED, WE NEED PRESSURE, VELOCITY, DENSITY AND GAMMA 
!



subroutine dt_function_2d(level, time, &
					phi,ph_lo,ph_hi, &
     dx, dt) bind(C, name="dt_function_2d")



  use global_var_module
  
  implicit none

  integer, intent(in) :: level
  double precision, intent(in) :: time
  integer, intent(in) :: ph_lo(2), ph_hi(2)
  double precision, intent(in) :: phi(ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2), nVar)
  double precision, intent(in) :: dx(2)
  double precision, intent(out) :: dt

  double precision :: soundspeed
  double precision ::cxmax, cymax
  double precision :: ctemp
  
  integer :: i, j
  
  cxmax = -10.0
  cymax = -10.0
  
  
!write(*,*) "dx", dx



	
  !
  !IN X DIRECTION 
  !
  
  
  !$omp parallel do private(i,j, soundspeed, ctemp) collapse(2)
  do j = ph_lo(2), ph_hi(2)
	do i = ph_lo(1), ph_lo(1)
	
		soundspeed = sqrt(Gamma * (phi(i,j,7)/phi(i,j,1)))
		
		ctemp = phi(i,j,5 ) + soundspeed 
		
   !    write (*,*) "ctemp", ctemp," uspeed",  phi(i,j,5 ), "soundspeed", soundspeed
		
		if(ctemp .gt. cxmax)then
		
			cxmax = ctemp
			
		end if
	!	write(*,*)"cxmax", cxmax
	end do
	
  end do
  		
!$omp end parallel do
  !
  !IN Y DIRECTION 
  !
  
  
  !$omp parallel do private(i,j, soundspeed, ctemp) collapse(2)
  do j = ph_lo(2), ph_hi(2)
	do i = ph_lo(1), ph_lo(1)
	
		soundspeed = sqrt(Gamma * (phi(i,j,7)/phi(i,j,1)))
		
		ctemp = phi(i,j,6 ) + soundspeed 
		
!		write (*,*) " vspeed",  phi(i,j,6 ), "soundspeed", soundspeed
	
		if(ctemp .gt. cymax)then
		
			cymax = ctemp
		
		end if
		
	end do
  end do
  		
!$omp end parallel do

!write(*,*)"cxmax", cxmax , "cymax", cymax
  if(cxmax .gt. cymax)then
  
		dt = (dx(1)/cxmax)
		
		
  else
  
		dt = (dx(2)/cymax)
		
  end if  



end subroutine dt_function_2d
