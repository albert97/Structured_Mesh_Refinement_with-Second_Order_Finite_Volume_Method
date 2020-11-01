


subroutine FluxANDUpdate(DIRPARAM, time, lo, hi, &
     &                   uin , ui_lo, ui_hi, &
     &          	     uout, uo_lo, uo_hi, &
     &           	     flxx, fx_lo, fx_hi, &
     &           	     flxy, fy_lo, fy_hi, &
     &            	     dx,dt) bind(C, name="FluxANDUpdate")
  
  use compute_flux_module, only : compute_flux_2d
  use global_var_module
  use global_function_module

  implicit none

  integer, intent(in) :: lo(2), hi(2)
  double precision, intent(in) :: dx(2), dt, time
  integer, intent(in) :: ui_lo(2), ui_hi(2)
  integer, intent(in) :: uo_lo(2), uo_hi(2)
  integer, intent(in) :: fx_lo(2), fx_hi(2)
  integer, intent(in) :: fy_lo(2), fy_hi(2)
  
  integer , intent(in) ::  DIRPARAM 
  
  double precision, intent(in) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),nVar)
  double precision, intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),nVar)
  double precision, intent(out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),nVar)
  double precision, intent(out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),nVar)
  double precision :: W(nVarQ), QPASS(nVarQ)


  integer :: i, j
  double precision :: dtdx(2)

  dtdx = dt/dx
  
  ! call a function to compute flux
  
  call compute_flux_2d(DIRPARAM, lo, hi, dt, dx, &
                       uin, ui_lo, ui_hi, &
                       flxx, fx_lo, fx_hi, &
                       flxy, fy_lo, fy_hi)

  ! Do a conservative update
  
  !write(*,*)"statein density", uin(:,:,1)
 !if(DIRPARAM == 1)then
 ! write(*,*)"statein density", uin(:,:,1)
 !  write(*,*)"flxy", flxy(:,:,1)
 !  write(*,*)"flxx", flxx(:,:,1)
 !write(*,*)"dtdx",dtdx

 !end if 
  !  write(*,*)"flxx", flxx(:,:,1)
 ! write(*,*)"dtdx",dtdx
 
  !$omp parallel do private(i,j) collapse(2)
  do    j = lo(2),hi(2)
     do i = lo(1),hi(1)
     
        uout(i,j,1 : nVarQ) = uin(i,j,1 : nVarQ) + &
             ( (flxx(i,j,1 : nVarQ) - flxx(i+1,j,1 : nVarQ)) * dtdx(1) &
             + (flxy(i,j,1 : nVarQ) - flxy(i,j+1,1 : nVarQ)) * dtdx(2) )
              
     enddo
  enddo
 !$omp end parallel do
 ! if(DIRPARAM == 1)then
   
  ! write(*,*) "stateout density", uout(:,:,1)
 !end if 
  
  ! write(*,*)"stateout density", uout(:,:,1)
  ! Scale by face area in order to correctly reflx
   !  !$omp parallel do private(i,j) collapse(2)
  do    j = lo(2), hi(2)
     do i = lo(1), hi(1)+1
        flxx(i,j,:) = flxx(i,j,:) * ( dt * dx(2))
     enddo
  enddo
 !  !$omp end parallel do
   !Scale by face area in order to correctly reflx
!  !$omp parallel do private(i,j) collapse(2)
  do    j = lo(2), hi(2)+1 
     do i = lo(1), hi(1)
        flxy(i,j,:) = flxy(i,j,:) * (dt * dx(1))
     enddo
  enddo
 !  !$omp end parallel do
!write(*,*)"Conservative varibales",uin(:,:,1:nVarQ)
!write(*,*)"Conservative varibales",uout(:,:,1:nVarQ)
end subroutine FluxANDUpdate
