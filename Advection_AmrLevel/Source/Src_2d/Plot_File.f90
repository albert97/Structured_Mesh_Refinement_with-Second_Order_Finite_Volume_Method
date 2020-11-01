
subroutine Plot_File(level, time,&
			    phi, ph_lo, ph_hi,&
			   dx) bind(C, name="Plot_File")
				

   use global_var_module
   
   implicit none
		
   integer, intent(in) :: ph_lo(2), ph_hi(2)
   double precision, intent(in) :: phi(ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2), nVar)
   double precision, intent(in) :: dx(2)
   integer, intent(in) :: level, time 
   double precision :: xlocation
  double precision :: x, y, z
	integer :: i, j 
	
	
!print *, "ph_lo(1)", ph_lo(1), "ph_hi(1)",ph_hi(1)
!print *, "ph_lo(2)", ph_lo(2), "ph_hi(2)",ph_hi(2)
!print*, "no_cell", no_cell
!print *, "level", level


!Box one data output
if(ph_lo(1) == 0 .and. ph_hi(1)==(No_cell/4)-1) then
	if(ph_lo(2) == 0 .and. ph_hi(2)==( No_cell/4)-1) then
!	print*, "I am working here"
		if (time .lt. tStop) then
!	print*, "I am working here"
			if (level == 0) then
	!		print*, "I am working here"
				!Output data into a file
					
				open(1, file = "ArmRex_data-densityBox_1_Level0.dat", action="write",status="replace")
		
				do i = ph_lo(1), ph_hi(1)
					xlocation = i * dx(1)
					write (1, *) xlocation, 0,  phi(i, 0, 1)
					
				end do
				close(1)	
	
			end if 
	
		end if
	end if 
end if 

!Box two data output
if((ph_lo(1) == ( No_cell/4))  .and. ph_hi(1)== ( No_cell/2 -1)) then
	if(ph_lo(2) == 0 .and. ph_hi(2)==( No_cell/4)-1) then
		if (time .lt. tStop) then

			if (level== 0) then
		
			!Output data into a file
			open(1, file = "ArmRex_data-densityBox_2_Level0d.dat", action="write",status="replace")
	
			!	do j= ph_lo(2), ph_hi(2)
			do i = ph_lo(1), ph_hi(1)
				xlocation = i * dx(1)
				write (1, *) xlocation, 0,  phi(i, 0, 1)
			end do
			close(1)	
	
			end if 
	
		end if
	end if 
end if 


!Box three data output
if(ph_lo(1) == (No_cell/2 ).and. ph_hi(1)== (3* No_cell/4 -1)) then
	if(ph_lo(2) == 0 .and. ph_hi(2)==( No_cell/4)-1) then
	
		if (time .lt. tStop) then

				if (level == 0) then
		
					!Output data into a file
					open(1, file = "ArmRex_data-densityBox_3_Level0.dat", action="write",status="replace")
	
					do i = ph_lo(1), ph_hi(1)
						xlocation = i * dx(1)
						write (1, *) xlocation, 0,  phi(i, 0, 1)
					end do

					close(1)	
				
				end if 
		end if
	end if 
end if 	

!Box four data output
if(ph_lo(1) == (3* No_cell/4 ) .and. ph_hi(1)== ( No_cell-1)) then
	if(ph_lo(2) == 0 .and. ph_hi(2)==( No_cell/4)-1) then
	
		if (time .lt. tStop) then

			if (level == 0) then
		
				!Output data into a file
				open(1, file = "ArmRex_data-densityBox_4_Level0.dat", action="write",status="replace")
	
				do i = ph_lo(1), ph_hi(1)			
						xlocation = i * dx(1)
						write (1, *) xlocation, 0,  phi(i, 0, 1)
				end do

				close(1)	
	
			end if 
	
		end if
	end if 
end if

	
end subroutine Plot_File


