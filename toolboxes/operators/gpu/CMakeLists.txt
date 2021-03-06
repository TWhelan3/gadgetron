if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_GPUOPERATORS__)
endif ()

if(WIN32)
link_directories(${Boost_LIBRARY_DIRS})
endif()

include_directories(
  ${CUDA_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CUDA_INCLUDE_DIRS}
  )

cuda_add_library(gadgetron_toolbox_gpuoperators SHARED 
  ../generalOperator.h
  ../linearOperator.h
  cuPartialDerivativeOperator.h
  cuLaplaceOperator.h
  cuTvOperator.h
  cuTv1dOperator.h
  cuConvolutionOperator.h
  cuPartialDerivativeOperator.cu
  cuPartialDerivativeOperator2.cu
  cuLaplaceOperator.cu
  cuTvOperator.cu
  cuTv1dOperator.cu
  cuConvolutionOperator.cu
  )

set_target_properties(gadgetron_toolbox_gpuoperators PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_gpuoperators 
  gadgetron_toolbox_gpucore 
  gadgetron_toolbox_gpunfft
  ${CUDA_LIBRARIES}
  ${CUDA_CUBLAS_LIBRARIES} 
  )

install(TARGETS gadgetron_toolbox_gpuoperators DESTINATION lib COMPONENT main)

install(FILES 
  cuImageOperator.h
  cuDiagonalOperator.h
  cuDiagonalSumOperator.h
  cuPartialDerivativeOperator.h
  cuPartialDerivativeOperator2.h
  cuConvolutionOperator.h
  cuLaplaceOperator.h
  cuTvOperator.h
  cuTvPicsOperator.h
  cuTv1dOperator.h
  cuDownsampleOperator.h
  cuFFTOperator.h
  cuUpsampleOperator.h
  hoCuIdentityOperator.h
  hoCuPartialDerivativeOperator.h
  hoCuTvOperator.h
  hoCuTvPicsOperator.h
  hoCuEncodingOperatorContainer.h
  gpuoperators_export.h
  hoCuOperator.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
