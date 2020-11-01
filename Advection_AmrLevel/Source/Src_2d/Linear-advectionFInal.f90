

subroutine advect(time, lo, hi, &
     &            uin , ui_lo, ui_hi, &
     &            uout, uo_lo, uo_hi, &
     &            vx  , vx_lo, vx_hi, &
     &            vy  , vy_lo, vy_hi, &
     &            flxx, fx_lo, fx_hi, &
     &            flxy, fy_lo, fy_hi, &
     &            dx,dt, xsweep) bind(C, name="advect")
  
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use compute_flux_module, only : compute_flux_2d

  implicit none

  integer, intent(in) :: lo(2), hi(2)
  double precision, intent(in) :: dx(2), dt, time
  integer, intent(in) :: ui_lo(2), ui_hi(2)
  integer, intent(in) :: uo_lo(2), uo_hi(2)
  integer, intent(in) :: vx_lo(2), vx_hi(2)
  integer, intent(in) :: vy_lo(2), vy_hi(2)
  integer, intent(in) :: fx_lo(2), fx_hi(2)
  integer, intent(in) :: fy_lo(2), fy_hi(2)
  integer, intent(in) :: xsweep
  double precision, intent(in   ) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2))
  double precision, intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2))
  double precision, intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2))
  double precision, intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2))
  double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
  double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))

  integer  :: i, j
  integer :: glo(2), ghi(2)
  double precision :: dtdx(2), umax, vmax

  !Dummy array to store results from x direction sweep
 
  dtdx = dt/dx
  
  glo = lo - 1
  ghi = hi + 1
 

  ! check if CFL condition is violated.
  umax = maxval(abs(vx))
  vmax = maxval(abs(vy))
  if ( umax*dt .ge. dx(1) .or. &
       vmax*dt .ge. dx(2) ) then
     print *, "umax = ", umax, ", vmax = ", vmax, ", dt = ", dt, ", dx = ", dx
     call bl_error("CFL violation. Use smaller adv.cfl.")
  end if
  
 

!PRINT *, "lox= ", lo(1), "loy= ", lo(2), "hix= ", hi(1), "hiy= ", hi(2)
!Print result shows the dimension are not consistent and does not neccessary start from zero	

!PRINT *, "ui_loX= ",ui_lo(1), "ui_loY= ", ui_lo(2), "ui_hix= ",ui_hi(1), "ui_hi= ", ui_hi(2)

!PRINT *, "vx_loX= ",vx_lo(1), "vx_loY= ", vx_lo(2), "vx_hix= ",vx_hi(1), "vx_hi= ", vx_hi(2)

!PRINT *, "vy_loX= ",vy_lo(1), "vy_loY= ", vy_lo(2), "vy_hix= ",vy_hi(1), "vy_hi= ", vy_hi(2)
	
!PRINT *, "vx(i-1,j) = ", vx(lo(1)-1,j), "uin(i-1, j)", uin(lo(1)-1, j), " vx(i+1,j)",  vx(lo(1)+1,j), " uin(i+1, j) ",  uin(lo(1)+1, j)
		
 ! print the value of vx, vy and  dx 
 
   
!  PRINT *, "PRINT THE X VELOCITY "
!  do    j = vx_lo(2), vx_hi(2)
!     do i = vx_lo(1), vx_hi(1)      
!		PRINT *, vx(i,j)
!     end do
!     PRINT *
!  end do
!    PRINT *, "PRINT THE Y VELOCITY "
!   do    j = vy_lo(2), vy_hi(2)
!     do i = vy_lo(1), vy_hi(1)      
!		PRINT *, vy(i,j)
!     end do
!     PRINT *
!  end do
 
  
 
  ! call a function to compute flux
  call compute_flux_2d(lo, hi, dt, dx, &
                       uin, ui_lo, ui_hi, &
                       vx, vx_lo, vx_hi, &
                       vy, vy_lo, vy_hi, &
                       flxx, fx_lo, fx_hi, &
                       flxy, fy_lo, fy_hi, xsweep)


  
    ! Do a conservative update
  do    j = lo(2),hi(2)
     do i = lo(1),hi(1)
        uout(i,j) = uin(i,j) + &
             (flxx(i,j) - flxx(i+1,j)) * dtdx(1) &
             + (flxy(i,j) - flxy(i,j+1)) * dtdx(2) 
             
 !      PRINT *, "uout(i,j)", uout(i,j), "uin(i,j)", uin(i,j),  "(flxx(i,j)", flxx(i,j), "(flxx(i+1,j)", flxx(i+1,j), " dtdx(1)",  dtdx(1), " dtdx(2)",  dtdx(2), "flxy(i,j) ", flxy(i,j), "flxy(i,j+1) ", flxy(i,j+1) 
     enddo
  enddo

  ! Scale by face area in order to correctly reflx
  do    j = lo(2), hi(2)
     do i = lo(1), hi(1)+1
        flxx(i,j) = flxx(i,j) * ( dt * dx(2))
     enddo
  enddo
  
  ! Scale by face area in order to correctly reflx
  do    j = lo(2), hi(2)+1 
     do i = lo(1), hi(1)
        flxy(i,j) = flxy(i,j) * (dt * dx(1))
     enddo
  enddo
  		
    		
  !Scall bl_deallocate(ux_out)
  
end subroutine advect

