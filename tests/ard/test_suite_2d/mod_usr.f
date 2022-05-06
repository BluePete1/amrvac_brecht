module mod_usr

  use mod_ard
  implicit none

  integer, parameter :: dp = kind(0.0d0)

contains

  subroutine usr_init()
     integer :: i

     call set_coordinate_system('Cartesian')
     call ard_activate()

     usr_init_one_grid => gs_init
  end subroutine usr_init

  subroutine gs_init(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,&
     ixmax2,w,x)
     integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
         ixmin1,ixmin2,ixmax1,ixmax2
     double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
        1:ndim)
     double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
        1:nw)
     double precision :: x1, x2, l1, l2

     ! Utility variables
     x1 = xprobmin1 + 0.5d0 * (xprobmax1 - xprobmin1)
     l1 = xprobmax1 - xprobmin1
     x2 = xprobmin2 + 0.5d0 * (xprobmax2 - xprobmin2)
     l2 = xprobmax2 - xprobmin2

     ! Center interval
     w(ixmin1:ixmax1,ixmin2:ixmax2,u_) = 1.0d0
     w(ixmin1:ixmax1,ixmin2:ixmax2,v_) = 0.0d0
     where ( (abs(x(ixmin1:ixmax1,ixmin2:ixmax2,&
         1) - x1) < 0.1d0 * l1) .and. (abs(x(ixmin1:ixmax1,ixmin2:ixmax2,&
         2) - x2) < 0.1d0 * l2) )
        w(ixmin1:ixmax1,ixmin2:ixmax2,u_) = 0.5d0
        w(ixmin1:ixmax1,ixmin2:ixmax2,v_) = 0.25d0
     endwhere
  end subroutine gs_init

end module mod_usr
