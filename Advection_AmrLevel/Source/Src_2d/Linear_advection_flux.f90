module compute_flux_module

  implicit none

  private

  public :: compute_flux_2d

contains

  subroutine compute_flux_2d(lo, hi, dt, dx, &
                             phi,ph_lo,ph_hi, &
                             umac,  u_lo,  u_hi, &
                             vmac,  v_lo,  v_hi, &
                             flxx, fx_lo, fx_hi, &
                             flxy, fy_lo, fy_hi, &
                             phix_1d, phiy_1d, phix, phiy, slope, glo, ghi)

    use slope_module, only: slopex, slopey

    integer, intent(in) :: lo(2), hi(2), glo(2), ghi(2)
    double precision, intent(in) :: dt, dx(2)
    integer, intent(in) :: ph_lo(2), ph_hi(2)
    integer, intent(in) ::  u_lo(2),  u_hi(2)
    integer, intent(in) ::  v_lo(2),  v_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    double precision, intent(in   ) :: phi (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2))
    double precision, intent(in   ) :: umac( u_lo(1): u_hi(1), u_lo(2): u_hi(2))
    double precision, intent(in   ) :: vmac( v_lo(1): v_hi(1), v_lo(2): v_hi(2))
    double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
    double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))
    double precision, dimension(glo(1):ghi(1),glo(2):ghi(2)) :: &
         phix_1d, phiy_1d, phix, phiy, slope
         
    integer :: i, j
    double precision :: hdtdx(2)

    hdtdx = 0.5*(dt/dx)
    
PRINT *, "dx1", dx(1), "dx2", dx(2)
PRINT *, "lox= ", lo(1), "loy= ", lo(2), "hix= ", hi(1), "hiy= ", hi(2)
!Print result shows the dimension are not consistent and does not neccessary start from zero	

PRINT *, "vx_loX= ",u_lo(1), "vx_loY= ", u_lo(2), "vx_hix= ",u_hi(1), "vx_hi= ", u_hi(2)

PRINT *, "vy_loX= ",v_lo(1), "vy_loY= ", v_lo(2), "vy_hix= ",v_hi(1), "vy_hi= ", v_hi(2)
	
PRINT *, " fx_lo(1) = ", fx_lo(1), "fx_hi(1)", fx_hi(1), "fx_lo(2)",fx_lo(2), " x_hi(2) ", fx_hi(2)

PRINT *, " fy_lo(1) = ", fy_lo(1), "fx_hi(1)", fy_hi(1), "fx_lo(2)",fy_lo(2), " x_hi(2) ", fy_hi(2)
	
    call slopex(glo, ghi, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi)

    ! compute phi on x faces using umac to upwind; ignore transverse terms
    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)  , hi(1)+1

          if (umac(i,j) .lt. 0.d0) then
             phix_1d(i,j) = phi(i  ,j) - (0.5d0 + hdtdx(1)*umac(i,j))*slope(i  ,j)
          else
             phix_1d(i,j) = phi(i-1,j) + (0.5d0 - hdtdx(1)*umac(i,j))*slope(i-1,j)
          end if

       end do
    end do

    call slopey(glo, ghi, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi)

    ! compute phi on y faces using umac to upwind; ignore transverse terms
    do    j = lo(2)  , hi(2)+1
       do i = lo(1)-1, hi(1)+1

          if (vmac(i,j) .lt. 0.d0) then
             phiy_1d(i,j) = phi(i,j  ) - (0.5d0 + hdtdx(2)*vmac(i,j))*slope(i,j  )
          else
             phiy_1d(i,j) = phi(i,j-1) + (0.5d0 - hdtdx(2)*vmac(i,j))*slope(i,j-1)
          end if

       end do
    end do

    ! update phi on x faces by adding in y-transverse terms
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)+1

          if (umac(i,j) .lt. 0.d0) then
             phix(i,j) = phix_1d(i,j) &
                  - hdtdx(2)*( 0.5d0*(vmac(i  ,j+1)+vmac(i  ,j)) * (phiy_1d(i  ,j+1)-phiy_1d(i  ,j)) )
          else
             phix(i,j) = phix_1d(i,j) &
                  - hdtdx(2)*( 0.5d0*(vmac(i-1,j+1)+vmac(i-1,j)) * (phiy_1d(i-1,j+1)-phiy_1d(i-1,j)) )
          end if

          ! compute final x-fluxes
          flxx(i,j) = phix(i,j)*umac(i,j)

       end do
    end do

    ! update phi on y faces by adding in x-transverse terms
    do    j = lo(2), hi(2)+1
       do i = lo(1), hi(1)

          if (vmac(i,j) .lt. 0.d0) then
             phiy(i,j) = phiy_1d(i,j) &
                  - hdtdx(1)*( 0.5d0*(umac(i+1,j  )+umac(i,j  )) * (phix_1d(i+1,j  )-phix_1d(i,j  )) )
          else
             phiy(i,j) = phiy_1d(i,j) &
                  - hdtdx(1)*( 0.5d0*(umac(i+1,j-1)+umac(i,j-1)) * (phix_1d(i+1,j-1)-phix_1d(i,j-1)) )
          end if

          ! compute final y-fluxes
          flxy(i,j) = phiy(i,j)*vmac(i,j)

       end do
    end do

  end subroutine compute_flux_2d
  
end module compute_flux_module
module compute_flux_module

  implicit none

  private

  public :: compute_flux_2d

contains

  subroutine compute_flux_2d(lo, hi, dt, dx, &
                             uin,uin_lo,uin_hi, &
                             umac,  u_lo,  u_hi, &
                             vmac,  v_lo,  v_hi, &
                             flxx, fx_lo, fx_hi, &
                             flxy, fy_lo, fy_hi, xsweep)
   ! use slope_module, only: slopex, slopey

    integer, intent(in) :: lo(2), hi(2)
    double precision, intent(in) :: dt, dx(2)
    integer, intent(in) :: uin_lo(2), uin_hi(2)
    integer, intent(in) ::  u_lo(2),  u_hi(2)
    integer, intent(in) ::  v_lo(2),  v_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    logical, intent(in) :: xsweep
    !dummy array to store values for extended uin in x direction 
    !dummy array to store values for extended uin in y direction 
    double precision :: uiny_extend(uin_lo(1):uin_hi(1)+1,uin_lo(2)-1:uin_hi(2)+1)
    double precision, intent(in   ) :: uin(uin_lo(1):uin_hi(1), uin_lo(2): uin_hi(2))
    double precision, intent(in   ) :: umac( u_lo(1): u_hi(1), u_lo(2): u_hi(2))
    double precision, intent(in   ) :: vmac( v_lo(1): v_hi(1), v_lo(2): v_hi(2))
    double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
    double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))
!    double precision, intent(  out) :: fluxR(fluxR_lo(1):fluxR_hi(1),fluxR_lo(2):fluxR_hi(2))
   
    integer :: i, j
    double precision :: dtdx(2)
   ! double precision :: CFL 
    double precision :: c1,c2
    dtdx = (dt/dx)
  
	
	! Fil lthe value into the extended array of Uin

	! Compute the left flux 
	
	if ( xsweep .eqv. .true.) then


	   do  j = lo(2), hi(2)
			do i = lo(1), hi(1)+1
		
				c1 = dtdx(1) * umac(i,j)
				c2 = dtdx(1) * umac(i+1,j)
				flxx(i,j) = (((1+c1)/(2.d0*c1)) * (umac(i,j) * uin(i, j ))) + (((c2-1)/(2.d0*c2)) * (umac(i+1,j) * uin(i+1, j )))
				flxy(i,j) = 0.d0
			end do
		end do
	else

		
	   do  j = lo(2), hi(2)+1
			do i = lo(1), hi(1)
			! compute y-fluxes
				c1 = dtdx(2) * vmac(i,j)
				c2 = dtdx(2) * vmac(i,j+1)
				flxy(i,j) =  (((1+c1)/(2.d0*c1)) * (vmac(i,j) * uin(i, j ))) + (((c2-1)/(2.d0*c2)) * (umac(i,j+1) * uin(i, j+1)))
		
			end do
	  end do
  end if
  end subroutine compute_flux_2d
  
end module compute_flux_module
module compute_flux_module

  implicit none

  private

  public :: compute_flux_2d

contains

  subroutine compute_flux_2d(lo, hi, dt, dx, &
                             uin,uin_lo,uin_hi, &
                             umac,  u_lo,  u_hi, &
                             vmac,  v_lo,  v_hi, &
                             flxx, fx_lo, fx_hi, &
                             flxy, fy_lo, fy_hi, xsweep)
  !as  use slope_module, only: slopex, slopey

    integer, intent(in) :: lo(2), hi(2)
    double precision, intent(in) :: dt, dx(2)
    integer, intent(in) :: uin_lo(2), uin_hi(2)
    integer, intent(in) ::  u_lo(2),  u_hi(2)
    integer, intent(in) ::  v_lo(2),  v_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    logical, intent(in) :: xsweep
    !dummy array to store values for extended uin in x direction 
    !dummy array to store values for extended uin in y direction 
    double precision :: uiny_extend(uin_lo(1):uin_hi(1)+1,uin_lo(2)-1:uin_hi(2)+1)
    double precision, intent(in   ) :: uin(uin_lo(1):uin_hi(1), uin_lo(2): uin_hi(2))
    double precision, intent(in   ) :: umac( u_lo(1): u_hi(1), u_lo(2): u_hi(2))
    double precision, intent(in   ) :: vmac( v_lo(1): v_hi(1), v_lo(2): v_hi(2))
    double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
    double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))
!    double precision, intent(  out) :: fluxR(fluxR_lo(1):fluxR_hi(1),fluxR_lo(2):fluxR_hi(2))
   
    integer :: i, j
    double precision :: dtdx(2)
    double precision :: CFL 
    double precision :: c1,c2
    dtdx = (dt/dx)
  
	
	! Fil lthe value into the extended array of Uin

	! Compute the left flux 
	
	if ( xsweep .eqv. .true.) then


	   do  j = lo(2), hi(2)
			do i = lo(1), hi(1)+1
		
				c1 = dtdx(1) * umac(i,j)
				c2 = dtdx(1) * umac(i+1,j)
				flxx(i,j) = (((1+c1)/(2.d0*c1)) * (umac(i,j) * uin(i, j ))) + (((c2-1)/(2.d0*c2)) * (umac(i+1,j) * uin(i+1, j )))
				flxy(i,j) = 0.d0
			end do
		end do
	else

		
	   do  j = lo(2), hi(2)+1
			do i = lo(1), hi(1)
			! compute y-fluxes
				c1 = dtdx(2) * vmac(i,j)
				c2 = dtdx(2) * vmac(i,j+1)
				flxy(i,j) =  (((1+c1)/(2.d0*c1)) * (vmac(i,j) * uin(i, j ))) + (((c2-1)/(2.d0*c2)) * (umac(i,j+1) * uin(i, j+1)))
				flxx(i,j) = 0.d0
			end do
	  end do
  end if
  end subroutine compute_flux_2d
  
end module compute_flux_module

module compute_flux_module

  implicit none

  private

  public :: compute_flux_2d

contains

  subroutine compute_flux_2d(lo, hi, dt, dx, &
                             uin,uin_lo,uin_hi, &
                             umac,  u_lo,  u_hi, &
                             vmac,  v_lo,  v_hi, &
                             flxx, fx_lo, fx_hi, &
                             flxy, fy_lo, fy_hi, xsweep)
   ! use slope_module, only: slopex, slopey

    integer, intent(in) :: lo(2), hi(2)
    double precision, intent(in) :: dt, dx(2)
    integer, intent(in) :: uin_lo(2), uin_hi(2)
    integer, intent(in) ::  u_lo(2),  u_hi(2)
    integer, intent(in) ::  v_lo(2),  v_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
	integer, intent(in) :: xsweep
    !dummy array to store values for extended uin in x direction 
    !dummy array to store values for extended uin in y direction 
    double precision, intent(in   ) :: uin(uin_lo(1):uin_hi(1), uin_lo(2): uin_hi(2))
    double precision, intent(in   ) :: umac( u_lo(1): u_hi(1), u_lo(2): u_hi(2))
    double precision, intent(in   ) :: vmac( v_lo(1): v_hi(1), v_lo(2): v_hi(2))
    double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
    double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))
!    double precision, intent(  out) :: fluxR(fluxR_lo(1):fluxR_hi(1),fluxR_lo(2):fluxR_hi(2))
   
    integer :: i, j
    double precision :: dtdx(2)
   ! double precision :: CFL 
    double precision :: c1,c2
    dtdx = (dt/dx)
  
!PRINT *, "lox= ", lo(1), "loy= ", lo(2), "hix= ", hi(1), "hiy= ", hi(2)
!Print result shows the dimension are not consistent and does not neccessary start from zero	

!PRINT *, "uin_loX= ",uin_lo(1), "uin_loY= ", uin_lo(2), "ui_hix= ",uin_hi(1), "ui_hi= ", uin_hi(2)

!PRINT *, "fx_loX= ",fx_lo(1), "fx_loY= ", fx_lo(2), "fx_hix= ",fx_hi(1), "fx_hi= ", fx_hi(2)

!PRINT *, "fy_loX= ",fy_lo(1), "fy_loY= ", fy_lo(2), "fy_hix= ",fy_hi(1), "fy_hi= ", fy_hi(2)
	
	
!	print*, "I am working here1"
	  do  j = lo(2), hi(2)
		do i = lo(1), hi(1)+1
			flxx( i, j ) = 0.0
		end do
	 end do
	 
!	print*, "I am working here2"
	  do  j = lo(2), hi(2)+1
		do i = lo(1), hi(1)
			flxy( i,j ) = 0.0
		end do
	  end do
!print*, "I am working here3"
	!print*, xsweep
	! Fil lthe value into the extended array of Uin

	! Compute the left flux 
	
	if ( xsweep .eq. 0) then
! 	print*, "I am working here5"
	   do  j = lo(2), hi(2)
			do i = lo(1), hi(1)+1
			
		!		print*, "dtdx(1)", dtdx(1), "umac(i,j)", umac(i,j), "umac(i+1,j)", umac(i+1,j), "uin(i, j )", uin(i, j ), "uin(i+1, j )", uin(i+1, j)
				c1 = dtdx(1) * umac(i,j)
				c2 = dtdx(1) * umac(i+1,j)
				flxx(i,j) = (((1+c1)/(2.d0*c1)) * (umac(i,j) * uin(i, j ))) + (((c2-1)/(2.d0*c2)) * (umac(i+1,j) * uin(i+1, j )))
		!		print *, "flxx(i,j)", flxx(i,j), "((1+c1)/(2.d0*c1)) ", ((1+c1)/(2.d0*c1)) , "((c2-1)/(2.d0*c2))", ((c2-1)/(2.d0*c2))
				!flxy(i,j) = 0.0
			!	flxy(i,j) =  (((1+c1)/(2.d0*c1)) * (vmac(i,j) * uin(i, j ))) + (((c2-1)/(2.d0*c2)) * (vmac(i,j+1) * uin(i, j+1)))
			end do
		end do
		
	else if (xsweep .eq. 1) then
! 	print*, "I am working here5"
	   do  j = lo(2), hi(2)+1
			do i = lo(1), hi(1)
			! compute y-fluxes
				c1 = dtdx(2) * vmac(i,j)
				c2 = dtdx(2) * vmac(i,j+1)
			!	flxx(i,j) = (((1+c1)/(2.d0*c1)) * (umac(i,j) * uin(i, j ))) + (((c2-1)/(2.d0*c2)) * (umac(i+1,j) * uin(i+1, j )))
				flxy(i,j) =  (((1+c1)/(2.d0*c1)) * (vmac(i,j) * uin(i, j ))) + (((c2-1)/(2.d0*c2)) * (vmac(i,j+1) * uin(i, j+1)))
			!flxx(i,j) = 0.0
				!print *, "flxy(i,j)", flxy(i,j)
			end do
	  end do
  end if
  end subroutine compute_flux_2d
  
end module compute_flux_module

