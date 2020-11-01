

subroutine get_face_velocity(level, time, &
     vx, vx_l1, vx_l2, vx_h1, vx_h2, &
     dx, prob_lo) bind(C, name="get_face_velocity")

  use probdata_module, only : adv_velL,   adv_velR
    
    
    
    
  implicit none

  integer, intent(in) :: level
  double precision, intent(in) :: time
  integer, intent(in) :: vx_l1, vx_l2, vx_h1, vx_h2
  double precision, intent(out) :: vx(vx_l1:vx_h1,vx_l2:vx_h2)
  double precision, intent(in) :: dx(2), prob_lo(2)

  vx = adv_vel(1)


end subroutine get_face_velocity
