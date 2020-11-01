
!----------------------------------------------------------------------!
!																	   !
!		THIS FUNCTION IS USED TO DISPLAY THE STATE VECTOR              !
!																	   !
!																	   !	
!----------------------------------------------------------------------!




subroutine Display_vector(level, time, lo, hi, &
					phi,ph_lo,ph_hi ) bind(C, name="Display_vector")



  use global_var_module
  
  implicit none

  integer, intent(in) :: level
  double precision, intent(in) :: time
  integer, intent(in) :: ph_lo(2), ph_hi(2)
  integer, intent(in) :: lo(2), hi(2)
  double precision, intent(in) :: phi(ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2), nVar)


  integer :: i, j, k

	do j = lo(2), hi(2)
		do i = lo(1), lo(1)	
			write(*,*) phi(i, j,1) ,   i , j 
		end do
		write(*,*)
    end do




end subroutine Display_vector
