/** \file   FlippingCartesianReconReferencePrepGadget.h
    \brief  This is the class gadget for both 2DT and 3DT cartesian reconstruction to prepare the reference data, working on the IsmrmrdReconData.
    \author Hui Xue
*/

#pragma once

#include <complex>
#include "gadgetron_mricore_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/meta.h"
#include "GadgetronTimer.h"

#include "hoNDArray_utils.h"

// #include "gtPlusIOAnalyze.h"

#include "GadgetStreamController.h"

#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"

#include "mri_core_data.h"
#include "mri_core_utility.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE FlippingCartesianReconReferencePrepGadget : public Gadget1<IsmrmrdReconData>
    {
    public:
        GADGET_DECLARE(FlippingCartesianReconReferencePrepGadget);

        typedef Gadget1<IsmrmrdReconData> BaseClass;

        FlippingCartesianReconReferencePrepGadget();
        ~FlippingCartesianReconReferencePrepGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------

        /// ------------------------------------------------------------------------------------
        /// debug and timing
        GADGET_PROPERTY(verbose, bool, "Whether to print more information", false);
        GADGET_PROPERTY(debug_folder, std::string, "If set, the debug output will be written out", "");
        GADGET_PROPERTY(perform_timing, bool, "Whether to perform timing on some computational steps", false);

        /// ref preparation
        /// whether to average all N for ref generation
        /// for the interleaved mode, the sampling times will be counted and used for averaging
        /// it is recommended to set N as the interleaved dimension
        GADGET_PROPERTY(average_all_ref_N, bool, "Whether to average all N for ref generation", false);
        /// whether to average all S for ref generation
        GADGET_PROPERTY(average_all_ref_S, bool, "Whether to average all S for ref generation", false);

	//Whether to flip in Readout Direction
	GADGET_PROPERTY(flip_RO, int, "Whether to flip in RO direction", 0);
	//Whether to flip in Encoding(1) Direction
	GADGET_PROPERTY(flip_E1, int, "Whether to flip in E1 direction", 0);
	//Whether to flip in Encoding(2) Direction
	//--NOT IMPLEMENTED YET--
	//GADGET_PROPERTY(flip_E2, bool, "Whether to flip in E2 direction", false);
	

    protected:

        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

        // number of encoding spaces in the protocol
        size_t num_encoding_spaces_;

        // for every encoding space
        // calibration mode
        std::vector<Gadgetron::ismrmrdCALIBMODE> calib_mode_;

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------

        // number of times the process function is called
        size_t process_called_times_;

        // --------------------------------------------------
        // variables for debug and timing
        // --------------------------------------------------

        // in verbose mode, more info is printed out
        bool verbose_;

        //// debug folder
        //std::string debug_folder_full_path_;

        // clock for timing
        Gadgetron::GadgetronTimer gt_timer_local_;
        Gadgetron::GadgetronTimer gt_timer_;

        //// exporter
        //Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;

        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1);
    };
}
