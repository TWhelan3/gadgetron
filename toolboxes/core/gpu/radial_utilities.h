#pragma once
#include "gpucore_export.h"

#include "cuNDArray.h"
#include "vector_td.h"

#include <boost/smart_ptr.hpp>

namespace Gadgetron{
// Compute fixed angle radial trajectory in the normalized range [-1/2;1/2]
template<class REAL> EXPORTGPUCORE boost::shared_ptr< cuNDArray< typename reald<REAL,2>::Type > >
compute_radial_trajectory_fixed_angle_2d( unsigned int num_samples_per_profile, unsigned int num_profiles_per_frame, unsigned int num_frames, REAL angular_offset = REAL(0) );

// Compute golden ratio radial trajectory in the normalized range [-1/2;1/2]
template<class REAL> EXPORTGPUCORE boost::shared_ptr< cuNDArray< typename reald<REAL,2>::Type > >
compute_radial_trajectory_golden_ratio_2d( unsigned int num_samples_per_profile, unsigned int num_profiles_per_frame, unsigned int num_frames, unsigned int profile_offset = 0 );

// Compute fixed angle radial density compensation weights (a function of the chose reconstruction settings: matrix_size and oversampling factor)
template<class REAL> EXPORTGPUCORE boost::shared_ptr< cuNDArray<REAL> >
compute_radial_dcw_fixed_angle_2d( unsigned int num_samples_per_profile, unsigned int num_profiles, REAL alpha, REAL one_over_radial_oversampling_factor);

// Compute golden ratio radial density compensation weights (a function of the chose reconstruction settings: matrix_size and oversampling factor)
template<class REAL> EXPORTGPUCORE boost::shared_ptr< cuNDArray<REAL> >
compute_radial_dcw_golden_ratio_2d( unsigned int num_samples_per_profile, unsigned int num_profiles, REAL alpha, REAL one_over_radial_oversampling_factor, unsigned int profile_offset = 0 );
}
