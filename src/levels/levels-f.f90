module levels_f_mod

  use, intrinsic :: iso_c_binding
  use :: flcl_mod

  type, bind(C) :: c_levels_t
    integer(c_int32_t) :: n
    type(nd_array_t) :: nd_ltop
    type(nd_array_t) :: nd_chunk_start
    type(nd_array_t) :: nd_chunk_end
    type(nd_array_t) :: nd_chunk_ids
  end type c_levels_t

  type :: levels_t
    integer(c_int32_t) :: n
    integer(c_int32_t) :: ntop
    integer(c_int32_t) :: chunk_count
    integer(c_int32_t), dimension(:), pointer :: ltop
    integer(c_int32_t), dimension(:), pointer :: chunk_start
    integer(c_int32_t), dimension(:), pointer :: chunk_end
    integer(c_int32_t), dimension(:), pointer :: chunk_ids
  end type levels_t

end module levels_f_mod