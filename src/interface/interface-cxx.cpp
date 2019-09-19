#include "flcl-cxx.hpp"
#include "levels-cxx.hpp"
#include <vector>
#include <iostream>
#include <cstdlib>
#include <Kokkos_ScatterView.hpp>

using std::vector;
using std::cout;
using std::endl;


template<class T>
void list_amr_vector( vector<T>& amr_vector) {

  cout << "list_amr_vector: amr vector contains:";
  for (auto it=amr_vector.begin(); it<amr_vector.end(); it++) {
    cout << ' ' << *it;
  }
  cout << endl;

}

enum cell_loop_t { clt_all, clt_ltop, clt_chunks, clt_all_scatter };

template< class T>
void cell_loop( const c_levels_t& c_levs, const cell_loop_t loop_type, const T& t ) {
  if (loop_type == clt_all) {

    Kokkos::parallel_for("clt_all",c_levs.n, t);
  
  } else if(loop_type == clt_ltop) {

    auto ltop_host = flcl::view_from_ndarray<int32_t*>(c_levs.nd_ltop);
    auto ltop_dv = Kokkos::DualView<int32_t*>("ltop_dv", ltop_host.extent(0));
    for(int i = 0; i < ltop_host.extent(0);i++){
	    ltop_dv.h_view(i) = ltop_host(i);
    }

    auto ltop = ltop_dv.d_view;
    Kokkos::parallel_for("clt_ltop",ltop.extent(0),
      KOKKOS_LAMBDA (const size_t ii) {
        int32_t idx = ltop(ii)-1;
        t(idx);
      }
    );

  } else if (clt_chunks) {
    
    auto chunk_start_host = flcl::view_from_ndarray<int32_t*>(c_levs.nd_chunk_start);
    auto chunk_end_host = flcl::view_from_ndarray<int32_t*>(c_levs.nd_chunk_end);    
    auto chunk_start_dv = Kokkos::DualView<int32_t*>("chunk_start_dv", chunk_start_host.extent(0));
    auto chunk_end_dv = Kokkos::DualView<int32_t*>("chunk_end_dv", chunk_end_host.extent(0));
    for(int i = 0; i < chunk_start_host.extent(0);i++){
	    chunk_start_dv.h_view(i) = chunk_start_host(i);
	    chunk_end_dv.h_view(i) = chunk_end_host(i);
    }
    auto chunk_start = chunk_start_dv.d_view;
    auto chunk_end = chunk_end_dv.d_view;

    Kokkos::parallel_for("clt_chunks",chunk_start_host.extent(0),
      KOKKOS_LAMBDA (const size_t ii) {
        auto chunk_size = chunk_end(ii)-chunk_start(ii);
        for (int32_t jj = 0; jj < chunk_size; jj++) {
          auto idx = chunk_start(ii) + jj - 1;
          t(idx);
        }
      }      
    );

  } else if (loop_type == clt_all_scatter) {
    
    Kokkos::parallel_for("clt_all_scatter",c_levs.n, t);
  
  }
}

extern "C" {

void c_kokkos_initialize() {

  Kokkos::initialize();

}

void c_kokkos_finalize( void ) {

  Kokkos::finalize();

}

void c_refine_mesh( flcl_ndarray_t &nd_mesh_ids, flcl_ndarray_t &nd_mesh_status )
{

    auto mesh_ids = flcl::view_from_ndarray<int32_t*>(nd_mesh_ids);
    auto mesh_status = flcl::view_from_ndarray<int32_t*>(nd_mesh_status);

    
}

void c_refine_mesh_status( flcl_ndarray_t &nd_mesh_status, flcl_ndarray_t &nd_new_mesh_status, int32_t &dimension, int32_t &refine_count ) {

  srand (time(NULL)); // comment if want determinism

  vector<int32_t> amr_vector (nd_mesh_status.dims[0]);
  // cout << "refine_amr_vector: nd_mesh_status.dims[0] " << nd_mesh_status.dims[0] << endl;
  // vector<int32_t> new_amr_vector (nd_new_mesh_status.dims[0]);
  auto v_mesh_status = flcl::view_from_ndarray<int32_t*>(nd_mesh_status);
  auto v_new_mesh_status = flcl::view_from_ndarray<int32_t*>(nd_new_mesh_status);

  for (int ii = 0; ii < amr_vector.size(); ii++) {
    amr_vector[ii] = v_mesh_status(ii);
  }
  // list_amr_vector<int32_t>( amr_vector );

  int actually_refined = 0;
  auto length = amr_vector.size();
  int32_t block_size = pow(2,dimension);
  // cout << "refine_amr_vector: length " << length << endl;
  // cout << "refine_amr_vector: refine_count " << refine_count << endl;
  // cout << "refine_amr_vector: dimension " << dimension << endl;
  // cout << "refine_amr_vector: block_size " << block_size << endl;

  while ( actually_refined < refine_count ) {
    auto it = amr_vector.begin();
    auto length = amr_vector.size();
    // cout << "refine_amr_vector: current amr vector length " << length << endl;

    int32_t min_val = floor(.1 * length);
    int32_t range = floor(.8 * length);
    int32_t element_to_refine = ( rand() % range ) + min_val;
    // cout << "refine_amr_vector: refining element " << element_to_refine << endl;

    if ( amr_vector.at(element_to_refine) < 0 ) {
      // cout << "refine_amr_vector: which is already refined, so skipping" << endl;
      continue;
    }

    int32_t block_to_refine = floor( element_to_refine / block_size )+1;
    // cout << "refine_amr_vector: which is in block " << block_to_refine << endl;
    
    int32_t insert_after = ( ( block_to_refine ) * block_size ) - 1;
    // cout << "refine_amr_vector: refined cells inserted after cell " << insert_after << endl;

    auto value_to_refine = amr_vector.at(element_to_refine);
    // cout << "refine_amr_vector: value to be refined from " << value_to_refine << endl;
    // it = amr_vector.insert( it+(element_to_refine+1), block_size, value_to_refine+1 );
    for (int32_t ii = 0; ii < block_size; ii++ ){
      amr_vector.emplace(it+(insert_after+1), value_to_refine+1 );    
      it = amr_vector.begin();
    }
    amr_vector.at(element_to_refine) *= -1;

    // list_amr_vector<int32_t>( amr_vector );
    actually_refined++;

  }

  for (int ii = 0; ii < nd_new_mesh_status.dims[0]; ii++) {
    v_new_mesh_status(ii) = amr_vector[ii];
  }

  // amr_vector.insert(it, length, 1);
  // cout << "amr vector initialized to length " << length << endl;
}

void c_axpy_full( flcl::dualview_r64_1d_t**  v_x, flcl::dualview_r64_1d_t** v_y, const double &s_val, const c_levels_t &c_levs ) {
  flcl::dualview_r64_1d_t::t_dev_const array_x = (*v_x)->d_view;
  flcl::dualview_r64_1d_t::t_dev array_y = (*v_y)->d_view;
  
  cell_loop( c_levs, clt_all,  KOKKOS_LAMBDA(const int& ii) {
    array_y(ii) += s_val * array_x(ii);
  } );

  cout << "c_axpy_full: finished" << endl;
}

void c_axpy_ltop( flcl::dualview_r64_1d_t** v_x, flcl::dualview_r64_1d_t** v_y, const double &s_val, const c_levels_t &c_levs ) {

  flcl::dualview_r64_1d_t::t_dev_const array_x = (*v_x)->d_view;
  flcl::dualview_r64_1d_t::t_dev array_y = (*v_y)->d_view;

  cell_loop( c_levs, clt_ltop, KOKKOS_LAMBDA(const int& ii) {
    array_y(ii) += s_val * array_x(ii);
  } );
  cout << "c_axpy_ltop: finished" << endl;
}

void c_axpy_ltop_chunks( flcl::dualview_r64_1d_t** v_x, flcl::dualview_r64_1d_t** v_y, const double &s_val, const c_levels_t &c_levs ) {

  flcl::dualview_r64_1d_t::t_dev_const array_x = (*v_x)->d_view;
  flcl::dualview_r64_1d_t::t_dev array_y = (*v_y)->d_view;
  
  cell_loop( c_levs, clt_chunks,
    KOKKOS_LAMBDA(const int& ii) {
      array_y(ii) += s_val * array_x(ii);
    }
  );
  cout << "c_axpy_ltop_chunks: finished" << endl;
}

void c_axpy_full_scatter( flcl::dualview_r64_1d_t**  v_x, flcl::dualview_r64_1d_t** v_y, const double &s_val, const c_levels_t &c_levs ) {
  flcl::dualview_r64_1d_t::t_dev_const array_x = (*v_x)->d_view;
  flcl::dualview_r64_1d_t::t_dev array_y = (*v_y)->d_view;
  using view_type = typename flcl::dualview_r64_1d_t::t_dev;
  using data_type = typename view_type::data_type;
  using array_layout = typename view_type::array_layout;
  using device_type = typename view_type::device_type;
  using memory_traits = typename view_type::memory_traits;
  using scatterview_type = Kokkos::Experimental::ScatterView<data_type, array_layout, Kokkos::DefaultExecutionSpace, Kokkos::Experimental::ScatterSum>;

  scatterview_type scatter(array_y);
  cell_loop( c_levs, clt_all_scatter,  KOKKOS_LAMBDA(const int& ii) {
    auto access = scatter.access();
    // array_y(ii) += s_val * array_x(ii);
    access(ii).update(s_val * array_x(ii));
  } );
  contribute(array_y, scatter);
  cout << "c_axpy_full_scatter: finished" << endl;
}

}
