#include "BipolarCorrectMatlabGadget.h"

namespace Gadgetron{

int BipolarCorrectMatlabGadget::process(GadgetContainerMessage<IsmrmrdImageArray>* m1)
{
	// Initialize a string for matlab commands
	std::string cmd;

        // Create a mxArray of bytes for the ISMRMRD::ImageHeader
	// Not really used for Bipolar Correction, but this way I don't have to redo BaseGadget
	// Could be useful in future
	ISMRMRD::ImageHeader *img_hdr = m1->getObjectPtr()->headers_.get_data_ptr();
        mwSize img_hdr_dims[2] = {sizeof(ISMRMRD::ImageHeader), 1};
        mxArray *img_hdr_bytes = mxCreateNumericArray(2, img_hdr_dims, mxUINT8_CLASS, mxREAL);
        memcpy(mxGetData(img_hdr_bytes), img_hdr, sizeof(ISMRMRD::ImageHeader));

	/*IMAGE DATA*/
	hoNDArray< std::complex<float> >* im_arr_ptr = &(m1->getObjectPtr()->data_);

	if (!im_arr_ptr) {
		GDEBUG("Broken raw_data pointer\n");
		return GADGET_FAIL;
	}

	/*DIMENSIONS*/
	std::vector<size_t> imSize;
	im_arr_ptr->get_dimensions(imSize);

	int xres = imSize[0];
	int yres = imSize[1];
	int numSlices = imSize[2];
	int numChannels = imSize[3];
	int numEchos = imSize[4];

	GINFO("%d %d %d %d %d\n", imSize[0], imSize[1], imSize[2], imSize[3], imSize[4]);

	if (xres == 0) xres = 1;
	if (yres == 0) yres = 1;
	if (numSlices == 0) numSlices = 1;
	if (numChannels == 0) numChannels = 1;
	if (numEchos == 0) numEchos = 1;

	mwSize ndim = 5;

	// Create a mxArray for the Image data
	mxArray *img_data = mxCreateNumericArray(5,imSize.data(), mxSINGLE_CLASS, mxCOMPLEX);

	float *real_data = (float *)mxGetData(img_data);
	float *imag_data = (float *)mxGetImagData(img_data);
	unsigned long num_elements = im_arr_ptr->get_number_of_elements();
	std::complex<float>* data_ptr = im_arr_ptr->get_data_ptr();
	for (int i = 0; i < num_elements; i++) {
		real_data[i] = data_ptr[i].real();
		imag_data[i] = data_ptr[i].imag();
	}
	engPutVariable(engine_, "hdr_bytes", img_hdr_bytes);
	engPutVariable(engine_, "data", img_data);
	cmd = "Q = matgadget.run_process(2, hdr_bytes, data); matgadget.emptyQ();";
	GDEBUG("command to send: %s\n", cmd.c_str());
	send_matlab_command(cmd);

	// Get the size of the gadget's queue
	// Right now (Bipolar correction) it should be one
	mxArray *Q = engGetVariable(engine_, "Q");

	if (Q == NULL) {
		GDEBUG("Failed to get the Queue from matgadget\n");
		return GADGET_FAIL;
	}
	size_t qlen = mxGetNumberOfElements(Q);

	// Loop over the elements of the Q, reading one entry at a time
	// to get a structure with type, headerbytes, and data
	mwIndex idx;
	for (idx = 0; idx < qlen; idx++)
	{ //for now all we have to do is overwrite data

		mxArray *res_data = mxGetField(Q, idx, "data");

		float *real_data = (float *)mxGetData(res_data);
		float *imag_data = (float *)mxGetImagData(res_data);
		for (int i = 0; i < im_arr_ptr->get_number_of_elements(); i++) {
			data_ptr[i] = std::complex<float>(real_data[i],imag_data[i]);
		}

		if (this->next()->putq(m1) < 0) {
			GDEBUG("Failed to put Image message on queue\n");
			return GADGET_FAIL;
		}

		break;

	}

	// Match all mxCreate___s with mxDestroy___s
	mxDestroyArray(img_hdr_bytes);
	mxDestroyArray(img_data);

	// Match engGetVariable with mxDestroy___s
	mxDestroyArray(Q);

	return GADGET_OK;
}


GADGET_FACTORY_DECLARE(BipolarCorrectMatlabGadget)
}
