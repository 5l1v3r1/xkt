module interface_f_mod

  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env
  use :: flcl_mod
  
  implicit none

  public

  interface
    subroutine f_kokkos_initialize() &
      & bind(c, name='c_kokkos_initialize')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine f_kokkos_initialize
  end interface

  interface
    subroutine f_kokkos_finalize() &
      & bind(c, name='c_kokkos_finalize')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine f_kokkos_finalize
  end interface

  interface
    subroutine f_refine_mesh_status( nd_new_mesh_status, nd_mesh_status, dimension, refine_count ) &
      & bind(c, name="c_refine_mesh_status")
      use, intrinsic :: iso_c_binding
      use :: flcl_mod
      implicit none
      type(nd_array_t), intent(inout) :: nd_mesh_status
      type(nd_array_t), intent(inout) :: nd_new_mesh_status
      integer(c_int32_t), intent(inout) :: dimension
      integer(c_int32_t), intent(inout) :: refine_count
    end subroutine f_refine_mesh_status
  end interface

  interface
    subroutine f_axpy_full( v_x, v_y, s_val, c_levs ) &
      & bind(c, name="c_axpy_full")
      use, intrinsic :: iso_c_binding
      use :: flcl_mod
      use :: levels_f_mod
      implicit none
      type(c_ptr), intent(inout) :: v_x
      type(c_ptr), intent(inout) :: v_y
      real(c_double), intent(inout) :: s_val
      type(c_levels_t), intent(inout) :: c_levs
    end subroutine f_axpy_full
  end interface

  interface
    subroutine f_axpy_ltop( v_x, v_y, s_val, c_levs ) &
      & bind(c, name="c_axpy_ltop")
      use, intrinsic :: iso_c_binding
      use :: flcl_mod
      use :: levels_f_mod
      implicit none
      type(c_ptr), intent(inout) :: v_x
      type(c_ptr), intent(inout) :: v_y
      real(c_double), intent(inout) :: s_val
      type(c_levels_t), intent(inout) :: c_levs
    end subroutine f_axpy_ltop
  end interface

  interface
    subroutine f_axpy_ltop_chunks( v_x, v_y, s_val, c_levs ) &
      & bind(c, name="c_axpy_ltop_chunks")
      use, intrinsic :: iso_c_binding
      use :: flcl_mod
      use :: levels_f_mod
      implicit none
      type(c_ptr), intent(inout) :: v_x
      type(c_ptr), intent(inout) :: v_y
      real(c_double), intent(inout) :: s_val
      type(c_levels_t), intent(inout) :: c_levs
    end subroutine f_axpy_ltop_chunks
  end interface

  contains

    subroutine kokkos_initialize()
      use, intrinsic :: iso_c_binding
      implicit none
      call f_kokkos_initialize()
    end subroutine kokkos_initialize
    
    subroutine kokkos_finalize()
      use, intrinsic :: iso_c_binding
      implicit none
      call f_kokkos_finalize()
    end subroutine kokkos_finalize

    subroutine refine_mesh_status( p_mesh_status, p_new_mesh_status, dimension, refine_count )

      use, intrinsic :: iso_c_binding
      use :: flcl_mod
      implicit none
      integer(c_int32_t), dimension(:), intent(inout), pointer :: p_mesh_status
      integer(c_int32_t), dimension(:), intent(inout), pointer :: p_new_mesh_status
      integer(c_int32_t), intent(inout) :: dimension
      integer(c_int32_t), intent(inout) :: refine_count


      ! integer(c_int32_t), dimension(:), pointer :: p_temp_mesh_status
      type(nd_array_t) :: nd_mesh_status
      type(nd_array_t) :: nd_new_mesh_status

      nd_mesh_status = to_nd_array( p_mesh_status )
      nd_new_mesh_status = to_nd_array( p_new_mesh_status )
      ! write(*,*)'refine_mesh_status: mesh_status_dims(1) ', mesh_status_dims(1)
      call f_refine_mesh_status( nd_mesh_status, nd_new_mesh_status, dimension, refine_count )

      ! p_temp_mesh_status => p_mesh_status
      ! p_mesh_status => p_new_mesh_status
      ! nullify(p_temp_mesh_status)

    end subroutine refine_mesh_status

    subroutine axpy_full( v_x, v_y, s_val, c_levs )
      use, intrinsic :: iso_c_binding
      use :: flcl_mod
      use :: levels_f_mod
      implicit none
      real(c_double), intent(inout) :: s_val
      type(c_levels_t), intent(inout) :: c_levs
      type(c_ptr), intent(inout) :: v_x
      type(c_ptr), intent(inout) :: v_y

      call f_axpy_full( v_x, v_y, s_val, c_levs )
    end subroutine axpy_full

    subroutine axpy_ltop( v_x, v_y, s_val, c_levs )
      use, intrinsic :: iso_c_binding
      use :: flcl_mod
      use :: levels_f_mod
      implicit none
      real(c_double), intent(inout) :: s_val
      type(c_levels_t), intent(inout) :: c_levs

      type(c_ptr), intent(inout) :: v_x
      type(c_ptr), intent(inout) :: v_y

      call f_axpy_ltop( v_x, v_y, s_val, c_levs )
    end subroutine axpy_ltop

    subroutine axpy_ltop_chunks( v_x, v_y, s_val, c_levs )
      use, intrinsic :: iso_c_binding
      use :: flcl_mod
      use :: levels_f_mod
      implicit none
      real(c_double), intent(inout) :: s_val
      type(c_levels_t), intent(inout) :: c_levs

      type(c_ptr), intent(inout) :: v_x
      type(c_ptr), intent(inout) :: v_y

      call f_axpy_ltop_chunks( v_x, v_y, s_val, c_levs )
    end subroutine axpy_ltop_chunks


end module interface_f_mod
