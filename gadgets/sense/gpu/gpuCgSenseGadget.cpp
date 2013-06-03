#include "gpuCgSenseGadget.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_blas.h"
#include "cuNDArray_utils.h"
#include "Gadgetron.h"
#include "GadgetMRIHeaders.h"
#include "b1_map.h"
#include "GPUTimer.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "vector_td_utilities.h"

#include "hoNDArray_fileio.h"

namespace Gadgetron{

  gpuCgSenseGadget::gpuCgSenseGadget()
  : channels_(0)
  , device_number_(0)
  , number_of_iterations_(5)
  , cg_limit_(1e-6)
  , oversampling_(1.25)
  , kernel_width_(5.5)
  , kappa_(0.1)
  , is_configured_(false)
  , image_series_(0)
  , image_counter_(0)
  {
    matrix_size_ = uintd2(0,0);
    matrix_size_os_ = uintd2(0,0);
  }

  gpuCgSenseGadget::~gpuCgSenseGadget() {}

  int gpuCgSenseGadget::process_config( ACE_Message_Block* mb )
  {
    GADGET_DEBUG1("gpuCgSenseGadget::process_config\n");
    GPUTimer timer("gpuCgSenseGadget::process_config");

    device_number_ = get_int_value(std::string("deviceno").c_str());

    int number_of_devices = 0;
    if (cudaGetDeviceCount(&number_of_devices)!= cudaSuccess) {
      GADGET_DEBUG1( "Error: unable to query number of CUDA devices.\n" );
      return GADGET_FAIL;
    }

    if (number_of_devices == 0) {
      GADGET_DEBUG1( "Error: No available CUDA devices.\n" );
      return GADGET_FAIL;
    }

    if (device_number_ >= number_of_devices) {
      GADGET_DEBUG2("Adjusting device number from %d to %d\n", device_number_,  (device_number_%number_of_devices));
      device_number_ = (device_number_%number_of_devices);
    }

    if (cudaSetDevice(device_number_)!= cudaSuccess) {
      GADGET_DEBUG1( "Error: unable to set CUDA device.\n" );
      return GADGET_FAIL;
    }

    number_of_iterations_ = get_int_value(std::string("number_of_iterations").c_str());
    cg_limit_ = get_double_value(std::string("cg_limit").c_str());
    oversampling_ = get_double_value(std::string("oversampling").c_str());
    kernel_width_ = get_double_value(std::string("kernel_width").c_str());
    kappa_ = get_double_value(std::string("kappa").c_str());
    pass_on_undesired_data_ = get_bool_value(std::string("pass_on_undesired_data").c_str());
    image_series_ = this->get_int_value("image_series");

    boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));

    std::vector<long> dims;
    ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();
    if (e_seq.size() != 1) {
      GADGET_DEBUG2("Number of encoding spaces: %d\n", e_seq.size());
      GADGET_DEBUG1("This Gadget only supports one encoding space\n");
      return GADGET_FAIL;
    }

    ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
    ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
    ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

    if (!is_configured_) {

      cudaDeviceProp deviceProp;
      if( cudaGetDeviceProperties( &deviceProp, device_number_ ) != cudaSuccess) {
	GADGET_DEBUG1( "\nError: unable to query device properties.\n" );
	return GADGET_FAIL;
      }

      unsigned int warp_size = deviceProp.warpSize;

      channels_ = cfg->acquisitionSystemInformation().present() ?
	(cfg->acquisitionSystemInformation().get().receiverChannels().present() ? cfg->acquisitionSystemInformation().get().receiverChannels().get() : 1) : 1;

      matrix_size_ = uintd2(e_space.matrixSize().x(), e_space.matrixSize().y());

      GADGET_DEBUG2("Matrix size  : [%d,%d] \n", matrix_size_.vec[0], matrix_size_.vec[1]);

      matrix_size_os_ =
	uintd2(static_cast<unsigned int>(std::ceil((matrix_size_.vec[0]*oversampling_)/warp_size)*warp_size),
	       static_cast<unsigned int>(std::ceil((matrix_size_.vec[1]*oversampling_)/warp_size)*warp_size));

      GADGET_DEBUG2("Matrix size OS: [%d,%d] \n", matrix_size_os_.vec[0], matrix_size_os_.vec[1]);

      // Allocate encoding operator for non-Cartesian Sense
      E_ = boost::shared_ptr< cuNonCartesianSenseOperator<float,2> >( new cuNonCartesianSenseOperator<float,2>() );

      // Allocate preconditioner
      D_ = boost::shared_ptr< cuCgPreconditioner<float_complext> >( new cuCgPreconditioner<float_complext>() );

      // Allocate regularization image operator
      R_ = boost::shared_ptr< cuImageOperator<float_complext> >( new cuImageOperator<float_complext>() );
      R_->set_weight( kappa_ );

      // Setup solver
      cg_.set_encoding_operator( E_ );        // encoding matrix
      cg_.add_regularization_operator( R_ );  // regularization matrix
      cg_.set_preconditioner( D_ );           // preconditioning matrix
      cg_.set_max_iterations( number_of_iterations_ );
      cg_.set_tc_tolerance( cg_limit_ );
      cg_.set_output_mode( cuCgSolver<float_complext>::OUTPUT_WARNINGS);

      is_configured_ = true;
    }

    return GADGET_OK;
  }

  int gpuCgSenseGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage<SenseJob> * m2)
  {
    GADGET_DEBUG1("gpuCgSenseGadget::process");
    GPUTimer timer("gpuCgSenseGadget::process");

    if (!is_configured_) {
      GADGET_DEBUG1("\nData received before configuration complete\n");
      return GADGET_FAIL;
    }

    SenseJob* j = m2->getObjectPtr();

    //Let's first check that this job has the required stuff...
    if (!j->csm_host_.get() || !j->dat_host_.get() || !j->tra_host_.get() || !j->dcw_host_.get()) {
      GADGET_DEBUG1("Received an incomplete CGSense JOB\n");
      m1->release();
      return GADGET_FAIL;
    }

    unsigned int samples = j->dat_host_->get_size(0);
    unsigned int channels = j->dat_host_->get_size(1);
    unsigned int rotations = samples / j->tra_host_->get_number_of_elements();
    unsigned int frames = j->tra_host_->get_size(1)*rotations;

    if( samples%j->tra_host_->get_number_of_elements() ) {
      GADGET_DEBUG2("Mismatch between number of samples (%d) and number of k-space coordinates (%d).\nThe first should be a multiplum of the latter.\n", 
		    samples, j->tra_host_->get_number_of_elements());
      m1->release();
      return GADGET_FAIL;
    }

    boost::shared_ptr< cuNDArray<floatd2> > traj(new cuNDArray<floatd2> (j->tra_host_.get()));
    boost::shared_ptr< cuNDArray<float> > dcw(new cuNDArray<float> (j->dcw_host_.get()));
    boost::shared_ptr< cuNDArray<float_complext> > csm(new cuNDArray<float_complext> (j->csm_host_.get()));
    boost::shared_ptr< cuNDArray<float_complext> > device_samples(new cuNDArray<float_complext> (j->dat_host_.get()));

    std::vector<unsigned int> image_dims = to_std_vector(matrix_size_);
    image_dims.push_back(frames);
    
    E_->set_domain_dimensions(&image_dims);
    E_->set_codomain_dimensions(device_samples->get_dimensions().get());
    E_->set_dcw(dcw);
    E_->set_csm(csm);

    try{ E_->setup( matrix_size_, matrix_size_os_, static_cast<float>(kernel_width_) ); }
    catch (runtime_error& err){
      GADGET_DEBUG_EXCEPTION(err, "\nError: unable to setup encoding operator.\n" );
      return GADGET_FAIL;
    }
    
    try{ E_->preprocess(traj.get());}
    catch (runtime_error& err){
      GADGET_DEBUG_EXCEPTION(err,"\nError during cuOperatorNonCartesianSense::preprocess()\n");
      return GADGET_FAIL;
    }

    boost::shared_ptr< cuNDArray<float_complext> > reg_image(new cuNDArray<float_complext> (j->reg_host_.get()));
    R_->compute(reg_image.get());

    // Define preconditioning weights
    boost::shared_ptr< cuNDArray<float> > _precon_weights = sum(abs_square(csm.get()).get(), 2);
    boost::shared_ptr<cuNDArray<float> > R_diag = R_->get();
    *R_diag *= float(kappa_);
    *_precon_weights += *R_diag;
    R_diag.reset();
    reciprocal_sqrt_inplace(_precon_weights.get());	
    boost::shared_ptr< cuNDArray<float_complext> > precon_weights = real_to_complex<float_complext>( _precon_weights.get() );
    _precon_weights.reset();
    D_->set_weights( precon_weights );
	
    // Invoke solver
    // 

    boost::shared_ptr< cuNDArray<float_complext> > cgresult;

    try{
      cgresult = cg_.solve(device_samples.get());
    }
    catch (runtime_error& err){
      GADGET_DEBUG_EXCEPTION(err,"\nError during solve()\n");
      return GADGET_FAIL;
    }

    if (!cgresult.get()) {
      GADGET_DEBUG1("\nIterative_sense_compute failed\n");
      return GADGET_FAIL;
    }

    static int counter = 0;
    char filename[256];
    sprintf((char*)filename, "recon_%2d.real", counter);
    write_nd_array<float>( abs(cgresult.get())->to_host().get(), filename );
    counter++;

    m2->release();

    // Now pass on the reconstructed images

    for( unsigned int frame=0; frame<frames; frame++ ){
      
      GadgetContainerMessage< hoNDArray< std::complex<float> > > *cm2 =
	new GadgetContainerMessage< hoNDArray< std::complex<float> > >();
      
      GadgetContainerMessage<ISMRMRD::ImageHeader> *m = 
	(frame==frames-1) ? m1 : new GadgetContainerMessage<ISMRMRD::ImageHeader>();

      if( frame<frames-1 )
	*m->getObjectPtr() = *m1->getObjectPtr();

      m->cont(cm2);

      std::vector<unsigned int> img_dims(2);
      img_dims[0] = matrix_size_.vec[0];
      img_dims[1] = matrix_size_.vec[1];

      try {cm2->getObjectPtr()->create(&img_dims);}
      catch (runtime_error &err){
	GADGET_DEBUG_EXCEPTION(err,"\nUnable to allocate host image array");
	m->release();
	return GADGET_FAIL;
      }

      size_t data_length = prod(matrix_size_);

      cudaMemcpy(cm2->getObjectPtr()->get_data_ptr(),
		 cgresult->get_data_ptr()+frame*data_length,
		 data_length*sizeof(std::complex<float>),
		 cudaMemcpyDeviceToHost);

      cudaError_t err = cudaGetLastError();
      if( err != cudaSuccess ){
	GADGET_DEBUG2("\nUnable to copy result from device to host: %s", cudaGetErrorString(err));
	m->release();
	return GADGET_FAIL;
      }

      m->getObjectPtr()->matrix_size[0] = img_dims[0];
      m->getObjectPtr()->matrix_size[1] = img_dims[1];
      m->getObjectPtr()->matrix_size[2] = 1;
      m->getObjectPtr()->channels       = 1;

      if (this->next()->putq(m) < 0) {
	GADGET_DEBUG1("\nFailed to result image on to Q\n");
	m->release();
	return GADGET_FAIL;
      }
    }
    
    return GADGET_OK;
  }

  GADGET_FACTORY_DECLARE(gpuCgSenseGadget)
}
