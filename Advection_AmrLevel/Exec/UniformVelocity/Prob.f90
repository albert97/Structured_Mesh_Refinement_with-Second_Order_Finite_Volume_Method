



subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  double precision, intent(in) :: problo(*), probhi(*)

  integer untin,i

  namelist /fortin/adv_vel
!  namelist /fortin/ adv_rhoL
!  namelist /fortin/ adv_preL
!  namelist /fortin/ adv_velR
!  namelist /fortin/ adv_rhoR
!  namelist /fortin/ adv_preR
  
  
  !
  ! Build "probin" filename -- the name of file containing fortin namelist.
  !     
  integer maxlen
  parameter (maxlen=256)
  character probin*(maxlen)
  
  if (namlen .gt. maxlen) then
     write(6,*) 'probin file name too long'
     stop
  end if
  
  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do
  
  ! set the namelist default
  adv_vel(:) = 1.d0

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

end subroutine amrex_probinit


subroutine initdata(level, time, lo, hi, &
     phi, phi_lo, phi_hi, &
     dx, prob_lo) bind(C, name="initdata")

  implicit none
  integer, intent(in) :: level, lo(3), hi(3), phi_lo(3), phi_hi(3)
  double precision, intent(in) :: time
  double precision, intent(inout) :: phi(phi_lo(1):phi_hi(1), &
       &                                 phi_lo(2):phi_hi(2), &
       &                                 phi_lo(3):phi_hi(3), 8)
  double precision, intent(in) :: dx(3), prob_lo(3)
  double precision tempangle
  integer          :: dm
  integer          :: i,j,k, nVar
  double precision :: u, v, rho, p, Gamma
  double precision :: DfC, x, y, z
  
  nVAr = 8
  Gamma =1.4d0



print *, "prob_lo(2)", prob_lo(2), "prob_lo(1)",prob_lo(1)
print *, "ph_lo(1)", phi_lo(1), "ph_hi(1)",phi_hi(1)
print *, "ph_lo(2)", phi_lo(2), "ph_hi(2)",phi_hi(2)
print *, "lo(1)", lo(1), "hi(1)",hi(1)
print *, "lo(2)", lo(2), "hi(2)",hi(2)

	

if (phi_lo(3) .eq. 0 .and. phi_hi(3) .eq. 0) then

	dm = 2
	
	else
	
	dm = 3
	
end if

				
				
!------------------- Initialise Inner Domain  ---------------------!

 
  !$omp parallel do private(i,j,k,x,y,z) collapse(2)
	do k=lo(3),hi(3)
			do j=lo(2), hi(2)
			
				z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
				y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
				
				do i= lo(1), hi(1)
						
					x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
					
				!	write(*,*) "x", x, "y", y
					
					DfC = sqrt((x-1d0)**2 + (y-1d0)**2)
					
				!	write(*,*)"DfC", DfC 
											
					if(DfC .LE. 0.4d0) then
						u = 0.0
						v = 0.0
						rho = 1.0
						p = 1.0
						phi(i,j,k, 1) = rho
						phi(i,j,k, 2) = rho * u
						phi(i,j,k, 3) = rho * v
						phi(i,j,k, 4) = p /( Gamma - 1.d0 )  + 0.5d0 * rho * (u**2.d0 + v**2.d0)
						phi(i,j,k, 5) = u
						phi(i,j,k, 6) = v
						phi(i,j,k, 7) = p
						phi(i,j,k, 8) = p /( Gamma - 1.d0 )
					else
						u = 0.0
						v = 0.0
						rho = 0.125
						p = 0.1
						phi(i,j,k, 1) = rho
						phi(i,j,k, 2) = rho * u
						phi(i,j,k, 3) = rho * v
						phi(i,j,k, 4) = p /( Gamma - 1.d0 )  + 0.5d0 * rho * (u**2.d0 + v**2.d0)
						phi(i,j,k, 5) = u
						phi(i,j,k, 6) = v
						phi(i,j,k, 7) = p
						phi(i,j,k, 8) = p /( Gamma - 1.d0 )
								
					end if
							 
				end do
					
	
			end do
		end do 

	 !$omp end parallel do

	

end subroutine initdata
