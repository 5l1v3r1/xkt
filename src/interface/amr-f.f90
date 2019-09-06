module amr_f_mod

  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env
  use :: flcl_mod
  
  implicit none

  public

  contains

    subroutine initialize_mesh_status( m, mesh_status )
      
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in) :: m
      integer(c_int32_t), dimension(:), intent(inout), target, allocatable :: mesh_status

      integer :: ii

      allocate( mesh_status(m) )

      do ii = 1, m
        mesh_status(ii) = 1
      end do 

    end subroutine initialize_mesh_status

    subroutine initialize_mesh_status_r64( m, mesh_status, view_ptr )
      
      use, intrinsic :: iso_c_binding
      use :: flcl_mod, only: kokkos_allocate_dualview
      integer(c_int32_t), intent(in) :: m
      real(c_double), dimension(:), intent(inout), pointer :: mesh_status
      type(c_ptr) :: view_ptr
      integer :: ii
      ! type(c_ptr) :: fptr

      ! character(len=255) ::fword1
      ! fword1 = C_CHAR_"foo"//C_NULL_CHAR
 

      ! fptr = c_loc(fword1)

      !allocate( mesh_status(m) )
      call kokkos_allocate_dualview(mesh_status, view_ptr, "foo", m)

      do ii = 1, m
        mesh_status(ii) = 1.0
      end do 

    end subroutine initialize_mesh_status_r64

    subroutine list_mesh_status( mesh_status )
      
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), dimension(:), intent(inout) :: mesh_status

      integer :: ii

      write(*,*)'list_mesh_status: ',mesh_status(:)

    end subroutine list_mesh_status

    ! subroutine refine_mesh_status( mesh_status, mesh_levels )
      
    !   use, intrinsic :: iso_c_binding
    !   implicit none
    !   integer(c_int32_t), dimension(:), intent(inout), target :: mesh_status
    !   integer(c_int32_t), dimension(:), intent(inout), target :: mesh_levels
    !   integer :: ii

    !   ! mark cells as refined by one level
    !   ! calculate length of new mesh_status
    !   ! allocate temporary new mesh_status
    !   ! deallocate old mesh_status
    !   ! re-allocate mesh_status
    !   ! copy new mesh_status into re-allocated mesh_status

    ! end subroutine refine_mesh_status

    function get_refine_count ( mesh_status, percentage ) result( refine_count )
      
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int32_t),intent(in), dimension(:) :: mesh_status
      real(c_double) :: percentage
      integer(c_int32_t) :: refine_count

      integer(c_int32_t) :: length
      real(c_double) :: d_length
      real(c_double) :: d_refine_count

      length = size( mesh_status, 1, c_int32_t )
      d_length = real( length, c_double )
      d_refine_count = floor( d_length * percentage )
      refine_count = int( d_refine_count )
      
      ! write(*,*)'get_refine_count: amr vector refine count: ',refine_count

    end function get_refine_count

    function count_top_level( p_mesh_status ) result( ntop )
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int32_t), intent(in), dimension(:) :: p_mesh_status
      integer(c_int32_t) :: ntop

      integer :: ii 

      ntop = 0

      do ii = 1, size(p_mesh_status)
        if (p_mesh_status(ii) > 0) then
          ntop = ntop + 1
        end if
      end do 

      write(*,*)'ntop: ',ntop
    end function count_top_level

    subroutine populate_ltop( p_mesh_status, p_ltop )
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int32_t), intent(in), dimension(:) :: p_mesh_status
      integer(c_int32_t), intent(in), dimension(:), pointer :: p_ltop
      
      integer :: ii
      integer :: ntop

      ntop = 0

      do ii = 1, size(p_mesh_status)
        if (p_mesh_status(ii) > 0) then
          ntop = ntop + 1
          p_ltop(ntop) = ii
        end if
      end do 

    end subroutine

    function count_chunks( p_mesh_status ) result( nchunks )
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int32_t), intent(in), dimension(:) :: p_mesh_status
      integer(c_int32_t) :: nchunks

      integer :: ii 

      nchunks = 1

      do ii = 2, size(p_mesh_status)
        if ( (p_mesh_status(ii-1) < 0 ) .and. (p_mesh_status(ii) > 0 ) ) then
          nchunks = nchunks + 1
        end if
      end do 

      write(*,*)'nchunks: ',nchunks
    end function count_chunks

    subroutine  populate_chunks( p_mesh_status, chunk_start, chunk_end )
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int32_t), intent(in), dimension(:) :: p_mesh_status
      integer(c_int32_t), intent(inout), dimension(:) :: chunk_start
      integer(c_int32_t), intent(inout), dimension(:) :: chunk_end
      integer(c_int32_t) :: nchunks

      integer :: ii 

      nchunks = 1
      chunk_start(1) = 1

      do ii = 2, size(p_mesh_status)
        if ( (    p_mesh_status(ii-1) < 0 ) .and. (p_mesh_status(ii) > 0 ) ) then
          nchunks = nchunks + 1
          chunk_start(nchunks) = ii
        elseif ( (p_mesh_status(ii-1) > 0 ) .and. (p_mesh_status(ii) < 0 ) ) then
          chunk_end(nchunks) = ii-1
        elseif ( (p_mesh_status(ii-1) > 0 ) .and. (p_mesh_status(ii) > 0 ) ) then
          chunk_end(nchunks) = ii
        end if
      end do

      write(*,*)'nchunks: ',nchunks      
    end subroutine populate_chunks

end module amr_f_mod
