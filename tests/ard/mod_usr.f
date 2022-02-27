!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_rd, mod_mhd ...
  use mod_mhd

  implicit none

  ! Custom variables can be defined here
  ! ...

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_2.5D")

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see mod_usr_methods.t
    ! ...

    ! Choose independent normalization units if using dimensionless variables.
    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm-3

    ! Active the physics module: rho_activate(), rd_activate(), mhd_activate()
    call mrd_activate()

  end subroutine usr_init

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImax1, ixOmin1,ixOmax1, w, x)

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)

    ! index for do loop and cell faces
    integer :: idir, ixCmin1,ixCmax1

    ! Set initial values for w
    ! 1.d0 and 0.d0 are just examples should be adjusted to user's problem

    ! cell-center density
    w(ixOmin1:ixOmax1, rho_) = 1.d0
    ! cell-center velocity 
    w(ixOmin1:ixOmax1, mom(1)) = 0.d0
    w(ixOmin1:ixOmax1, mom(2)) = 0.d0
    ! cell-center pressure
    w(ixOmin1:ixOmax1, e_) = 1.d0
    if(stagger_grid) then
      ! set cell-face magnetic field B using CT method for divB control

      ! set B from vector potential (zero divB guaranteed) given 
      ! usr_init_vector_potential pointing to  a subroutine to set vector potential
      call b_from_vector_potential(ixGslo1,ixGshi1,ixImin1,ixImax1,ixOmin1,&
         ixOmax1,block%ws,x)

      ! or directly set cell-face B (divB maybe non-zero) as following:
      do idir=1,ndim
        ixCmin1=ixImin1;
        ixCmax1=ixImax1-kr(idir,1);
        ! cell-face B_idir
        block%ws(ixCmin1:ixCmax1,idir)=1.d0
      end do
      ! update cell-center B from cell-face B
      call mrd_face_to_center(ixOmin1,ixOmax1,block)
    else
      ! cell-center magnetic field
      w(ixOmin1:ixOmax1,mag(1)) = 1.d0
      w(ixOmin1:ixOmax1,mag(2)) = 1.d0
    end if

    ! convert primitive variables to conservative variables
    call mrd_to_conserved(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)

  end subroutine initial_conditions

  ! Extra routines can be placed here
  ! ...

end module mod_usr
