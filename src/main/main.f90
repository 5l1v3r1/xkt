program main
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env

  use :: flcl_mod
  use :: levels_f_mod

  use :: interface_f_mod
  use :: amr_f_mod

  implicit none

  type(levels_t) :: levs
  integer :: m
  integer :: ii
  integer :: refine_count
  integer :: refine_steps
  real(c_double) :: percentage
  integer(c_int32_t) ::dimension
  integer :: block_size
  integer(c_int32_t) :: ntop
  integer(c_int32_t) :: nchunks
  integer(c_int32_t), dimension(:), pointer :: mesh_ids
  integer(c_int32_t), dimension(:), allocatable, target :: mesh_status
  integer(c_int32_t), dimension(:), allocatable, target :: new_mesh_status
  integer(c_int32_t), dimension(:), pointer :: p_mesh_status
  integer(c_int32_t), dimension(:), pointer :: p_new_mesh_status
  ! integer(c_int32_t), dimension(:), allocatable, target :: ltop
  integer(c_int32_t), dimension(:), pointer :: p_ltop
  real(c_double), dimension(:), pointer :: array_x
  real(c_double), dimension(:), pointer :: array_y
  real(c_double) :: s_val
  type(c_ptr) :: v_x
  type(c_ptr) :: v_y

  type(c_levels_t) :: c_levs
  
  ! sensible values for dimension might be 1,2,3
  dimension = 2
  block_size = 2**dimension

  ! sensible values for m might be 1e4 through 1e6 (or so)
  m = 2000


  ! set up values for "refinement"
  percentage = 0.1
  refine_steps = 1

  ! alpha in axpy
  s_val = 2.0


  ! allocate( mesh_ids(m) )
  ! allocate( mesh_status(m) )



  call kokkos_initialize()

  call initialize_mesh_status( m, mesh_status )
  refine_count = get_refine_count ( mesh_status, percentage )
  call initialize_mesh_status( m+(refine_count*block_size), new_mesh_status )
  p_mesh_status => mesh_status
  p_new_mesh_status => new_mesh_status
  ! write(*,*)'main: refine_count ',refine_count
  ! write(*,*)'main: dimension ',dimension
  ! call list_mesh_status( p_mesh_status )
  ! call list_mesh_status( p_new_mesh_status )

  call refine_mesh_status( p_mesh_status, p_new_mesh_status, dimension, refine_count )

  ! call list_mesh_status( p_mesh_status )
  ! call list_mesh_status( p_new_mesh_status )
  
  do ii = 1, refine_steps
    refine_count = get_refine_count ( new_mesh_status, percentage )
    nullify(p_mesh_status)
    deallocate(mesh_status)
    mesh_status = new_mesh_status
    p_mesh_status => mesh_status
    deallocate(new_mesh_status)
    call initialize_mesh_status( size(p_mesh_status)+(refine_count*block_size), new_mesh_status )
    p_new_mesh_status => new_mesh_status
    call refine_mesh_status( p_mesh_status, p_new_mesh_status, dimension, refine_count )
  end do

  ! unroll
  nullify(p_mesh_status)
  deallocate(mesh_status)
  mesh_status = new_mesh_status
  p_mesh_status => mesh_status
  deallocate(new_mesh_status)
  levs%n = size(mesh_status)
  write (*,*)'levs%n ',levs%n



  ! enable print to screen for debugging
  ! call list_mesh_status( p_mesh_status )

  ! create ltop and chunks
  ntop = count_top_level( p_mesh_status )
  levs%ntop = ntop
  allocate(levs%ltop(ntop))
  p_ltop => levs%ltop
  call populate_ltop( p_mesh_status, p_ltop )

  ! enable for debugging
  ! call list_mesh_status( p_ltop )

  nchunks = count_chunks( p_mesh_status )
  allocate(levs%chunk_start(nchunks))
  allocate(levs%chunk_end(nchunks))
  allocate(levs%chunk_ids(nchunks))
  call populate_chunks( p_mesh_status, levs%chunk_start, levs%chunk_end )
  levs%chunk_count = nchunks
  write (*,*)'levs%chunk_count ',levs%chunk_count

  ! enable for debugging
  ! call list_mesh_status( levs%chunk_start )
  ! call list_mesh_status( levs%chunk_end )

  c_levs%nd_ltop = to_nd_array( levs%ltop )
  c_levs%nd_chunk_start = to_nd_array( levs%chunk_start )
  c_levs%nd_chunk_end = to_nd_array( levs%chunk_end )
  c_levs%nd_chunk_ids = to_nd_array( levs%chunk_ids )
  c_levs%n = levs%n
  
  ! allocate 'physics arrays'
  call initialize_mesh_status_r64( size(p_mesh_status), array_x, v_x)
  call initialize_mesh_status_r64( size(p_mesh_status), array_y, v_y)

  ! full axpy
  call axpy_full( v_x, v_y, s_val, c_levs )

  ! ltop axpy
  call axpy_ltop( v_x, v_y, s_val, c_levs )

  ! chunk axpy
  call axpy_ltop_chunks( v_x, v_y, s_val, c_levs )

  ! full axpy w scatterview
  call axpy_full_scatter( v_x, v_y, s_val, c_levs )

  call kokkos_deallocate_dualview( array_x, v_x )
  call kokkos_deallocate_dualview( array_y, v_y )

  call kokkos_finalize()


end program
