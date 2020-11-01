
subroutine Max_Gradient( uin , ui_lo, ui_hi, &
                        tag_lo,tag_hi,  &
                         level, dx, MAXG) bind(C, name="Max_Gradient")
                          
    use global_var_module
    
	implicit none
    
    integer, intent(in) :: tag_lo(2), tag_hi(2)
    integer, intent(in) :: ui_lo(2), ui_hi(2)
    double precision, intent(in) :: dx(2), level
    double precision, intent(in) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),nVar)
	double precision, intent(inout) :: MAXG
	double precision :: ax, ay
	double precision :: log_d_rho_mod
	integer          :: i, j
	
	
!	write(*,*) "un_lo", ui_lo, "ui_hi", ui_hi
!	write(*,*) " tag_lo",  tag_lo, "tag_hi",  tag_hi
	MAXG = 0.0
!	!$omp parallel do private(i,j, ax,ay,log_d_rho_mod)

    do	   j= tag_lo(2), tag_hi(2)
		do i =  tag_lo(1),  tag_hi(1)
		
			!Calculating the maxmium gradient using central difference
			
		    ax = (uin(i-1,j, 1)-(uin(i+1,j, 1)))/2.d0/dx(1)
			ay = (uin(i,j-1, 1)-uin(i,j+1, 1))/2.d0/dx(2)
!			write(*,*)"ax", ax, "ay",ay
			log_d_rho_mod = log(1 + sqrt( ax**2.d0 + ay**2.d0)/uin(i,j,1))
		!	write (*,*) "log_d_rho_mod", log_d_rho_mod
			if (log_d_rho_mod > MAXG) then
				MAXG = log_d_rho_mod
			end if 	
		
		end do
    end do                 
 !   !$omp end parallel do   
 ! write(*,*)"Max_Grad", MAXG                    
end subroutine  Max_Gradient

             
