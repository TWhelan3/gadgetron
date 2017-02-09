
#include "MultiChannelCartesianGrappaReconGadget.h"
#include <iomanip>

#include "mri_core_kspace_filter.h"
#include "hoNDArray_reductions.h"
#include "mri_core_def.h"
#include "mri_core_coil_map_estimation.h"
#include "mri_core_grappa.h"
#include "mri_core_utility.h"

/*
    The input is IsmrmrdReconData and output is single 2D or 3D ISMRMRD images

    If required, the gfactor map can be sent out

    If the  number of required destination channel is 1, the GrappaONE recon will be performed

    The image number computation logic is implemented in compute_image_number function, which can be overloaded
*/

namespace Gadgetron {

    MultiChannelCartesianGrappaReconGadget::MultiChannelCartesianGrappaReconGadget() : num_encoding_spaces_(1), process_called_times_(0)
    {
        gt_timer_.set_timing_in_destruction(false);
        gt_timer_local_.set_timing_in_destruction(false);
    }

    MultiChannelCartesianGrappaReconGadget::~MultiChannelCartesianGrappaReconGadget()
    {
    }

    int MultiChannelCartesianGrappaReconGadget::process_config(ACE_Message_Block* mb)
    {
        ISMRMRD::IsmrmrdHeader h;
        try
        {
            deserialize(mb->rd_ptr(), h);
        }
        catch (...)
        {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        if (!h.acquisitionSystemInformation)
        {
            GDEBUG("acquisitionSystemInformation not found in header. Bailing out");
            return GADGET_FAIL;
        }

        // -------------------------------------------------

        size_t NE = h.encoding.size();
        num_encoding_spaces_ = NE;
        GDEBUG_CONDITION_STREAM(verbose.value(), "Number of encoding spaces: " << NE);

        meas_max_idx_.resize(NE);
        acceFactorE1_.resize(NE, 1);
        acceFactorE2_.resize(NE, 1);
        calib_mode_.resize(NE, ISMRMRD_noacceleration);

        recon_obj_.resize(NE);

        size_t e;
        for (e = 0; e < h.encoding.size(); e++)
        {
            ISMRMRD::EncodingSpace e_space = h.encoding[e].encodedSpace;
            ISMRMRD::EncodingSpace r_space = h.encoding[e].reconSpace;
            ISMRMRD::EncodingLimits e_limits = h.encoding[e].encodingLimits;

            GDEBUG_CONDITION_STREAM(verbose.value(), "---> Encoding space : " << e << " <---");
            GDEBUG_CONDITION_STREAM(verbose.value(), "Encoding matrix size: " << e_space.matrixSize.x << " " << e_space.matrixSize.y << " " << e_space.matrixSize.z);
            GDEBUG_CONDITION_STREAM(verbose.value(), "Encoding field_of_view : " << e_space.fieldOfView_mm.x << " " << e_space.fieldOfView_mm.y << " " << e_space.fieldOfView_mm.z);
            GDEBUG_CONDITION_STREAM(verbose.value(), "Recon matrix size : " << r_space.matrixSize.x << " " << r_space.matrixSize.y << " " << r_space.matrixSize.z);
            GDEBUG_CONDITION_STREAM(verbose.value(), "Recon field_of_view :  " << r_space.fieldOfView_mm.x << " " << r_space.fieldOfView_mm.y << " " << r_space.fieldOfView_mm.z);

            meas_max_idx_[e].kspace_encode_step_1 = (uint16_t)e_space.matrixSize.y - 1;
            meas_max_idx_[e].set = (e_limits.set && (e_limits.set->maximum > 0)) ? e_limits.set->maximum : 0;
            meas_max_idx_[e].phase = (e_limits.phase && (e_limits.phase->maximum > 0)) ? e_limits.phase->maximum : 0;

            meas_max_idx_[e].kspace_encode_step_2 = (uint16_t)e_space.matrixSize.z - 1;

            meas_max_idx_[e].contrast = (e_limits.contrast && (e_limits.contrast->maximum > 0)) ? e_limits.contrast->maximum : 0;
            meas_max_idx_[e].slice = (e_limits.slice && (e_limits.slice->maximum > 0)) ? e_limits.slice->maximum : 0;
            meas_max_idx_[e].repetition = e_limits.repetition ? e_limits.repetition->maximum : 0;
            meas_max_idx_[e].slice = (e_limits.slice && (e_limits.slice->maximum > 0)) ? e_limits.slice->maximum : 0;
            meas_max_idx_[e].average = e_limits.average ? e_limits.average->maximum : 0;
            meas_max_idx_[e].segment = 0;

            if (!h.encoding[e].parallelImaging)
            {
                GDEBUG_STREAM("Parallel Imaging section not found in header");
                calib_mode_[e] = ISMRMRD_noacceleration;
                acceFactorE1_[e] = 1;
                acceFactorE2_[e] = 1;
            }
            else
            {
                ISMRMRD::ParallelImaging p_imaging = *h.encoding[0].parallelImaging;

                acceFactorE1_[e] = p_imaging.accelerationFactor.kspace_encoding_step_1;
                acceFactorE2_[e] = p_imaging.accelerationFactor.kspace_encoding_step_2;
                GDEBUG_CONDITION_STREAM(verbose.value(), "acceFactorE1 is " << acceFactorE1_[e]);
                GDEBUG_CONDITION_STREAM(verbose.value(), "acceFactorE2 is " << acceFactorE2_[e]);

                std::string calib = *p_imaging.calibrationMode;

                bool separate = (calib.compare("separate") == 0);
                bool embedded = (calib.compare("embedded") == 0);
                bool external = (calib.compare("external") == 0);
                bool interleaved = (calib.compare("interleaved") == 0);
                bool other = (calib.compare("other") == 0);

                calib_mode_[e] = Gadgetron::ISMRMRD_noacceleration;
                if (acceFactorE1_[e] > 1 || acceFactorE2_[e] > 1)
                {
                    if (interleaved)
                        calib_mode_[e] = Gadgetron::ISMRMRD_interleaved;
                    else if (embedded)
                        calib_mode_[e] = Gadgetron::ISMRMRD_embedded;
                    else if (separate)
                        calib_mode_[e] = Gadgetron::ISMRMRD_separate;
                    else if (external)
                        calib_mode_[e] = Gadgetron::ISMRMRD_external;
                    else if (other)
                        calib_mode_[e] = Gadgetron::ISMRMRD_other;
                }
            }
        }

        // ---------------------------------------------------------------------------------------------------------

        /*if (!debug_folder.value().empty())
        {
            Gadgetron::get_debug_folder_path(debug_folder.value(), debug_folder_full_path_);
            GDEBUG_CONDITION_STREAM(verbose.value(), "Debug folder is " << debug_folder_full_path_);
        }
        else
        {
            GDEBUG_CONDITION_STREAM(verbose.value(), "Debug folder is not set ... ");
        }*/

        return GADGET_OK;
    }

    int MultiChannelCartesianGrappaReconGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1)
    {
        process_called_times_++;

        IsmrmrdReconData* recon_bit_ = m1->getObjectPtr();
        if (recon_bit_->rbit_.size() > num_encoding_spaces_)
        {
            GWARN_STREAM("Incoming recon_bit has more encoding spaces than the protocol : " << recon_bit_->rbit_.size() << " instead of " << num_encoding_spaces_);
        }

        // for every encoding space
        for (size_t e = 0; e < recon_bit_->rbit_.size(); e++)
        {
            std::stringstream os;
            os << "_encoding_" << e;

            GDEBUG_CONDITION_STREAM(verbose.value(), "Calling " << process_called_times_ << " , encoding space : " << e);
            GDEBUG_CONDITION_STREAM(verbose.value(), "======================================================================");

            // ---------------------------------------------------------------

            if (recon_bit_->rbit_[e].ref_)
            {
                // ---------------------------------------------------------------
		
		
		
		        // after this step, the recon_obj_[e].ref_calib_ and recon_obj_[e].ref_coil_map_ are set

	        this->make_ref_coil_map(*recon_bit_->rbit_[e].ref_,*recon_bit_->rbit_[e].data_.data_.get_dimensions(), recon_obj_[e], e);

	        // after this step, coil map is computed and stored in recon_obj_[e].coil_map_
	        if(!flatmap.value())
			this->perform_coil_map_estimation(recon_bit_->rbit_[e], recon_obj_[e], e);
	        else
		{
			recon_obj_[e].coil_map_=recon_obj_[e].ref_coil_map_;
        		recon_obj_[e].coil_map_.fill(1);//create(*recon_bit_->rbit_[e].data_.data_.get_dimensions());
		}
                // after this step, recon_obj_[e].kernel_, recon_obj_[e].kernelIm_, recon_obj_[e].unmixing_coeff_ are filled
                // gfactor is computed too

                this->perform_calib(recon_bit_->rbit_[e], recon_obj_[e], e);


                // ---------------------------------------------------------------

                recon_bit_->rbit_[e].ref_ = boost::none;
            }

            if (recon_bit_->rbit_[e].data_.data_.get_number_of_elements() > 0)
            {
                // ---------------------------------------------------------------


                this->perform_unwrapping(recon_bit_->rbit_[e], recon_obj_[e], e);

                // ---------------------------------------------------------------

                this->compute_image_header(recon_bit_->rbit_[e], recon_obj_[e], e);

                // ---------------------------------------------------------------

                this->send_out_image_array(recon_bit_->rbit_[e], recon_obj_[e].recon_res_, e, image_series.value() + ((int)e + 1), GADGETRON_IMAGE_REGULAR);
            }

            recon_obj_[e].recon_res_.data_.clear();
            recon_obj_[e].gfactor_.clear();
            recon_obj_[e].recon_res_.headers_.clear();
            recon_obj_[e].recon_res_.meta_.clear();
        }

        m1->release();
        return GADGET_OK;
    }

    size_t MultiChannelCartesianGrappaReconGadget::compute_image_number(ISMRMRD::ImageHeader& imheader, size_t encoding, size_t CHA, size_t cha, size_t E2)
    {
        if (encoding >= meas_max_idx_.size())
        {
            GWARN_STREAM("encoding >= meas_max_idx_.size()");
            encoding = 0;
        }

        size_t SET = meas_max_idx_[encoding].set + 1;
        size_t REP = meas_max_idx_[encoding].repetition + 1;
        size_t PHS = meas_max_idx_[encoding].phase + 1;
        size_t SLC = meas_max_idx_[encoding].slice + 1;
        size_t CON = meas_max_idx_[encoding].contrast + 1;
        if (E2 == 0) E2 = 1;

        size_t imageNum = imheader.average*REP*SET*PHS*CON*SLC*E2*CHA + imheader.repetition*SET*PHS*CON*SLC*E2*CHA 
            + imheader.set*PHS*CON*SLC*E2*CHA + imheader.phase*CON*SLC*E2*CHA + imheader.contrast*SLC*E2*CHA + imheader.slice*E2*CHA + cha + 1;

        return imageNum;
    }

    int MultiChannelCartesianGrappaReconGadget::send_out_image_array(IsmrmrdReconBit& recon_bit, IsmrmrdImageArray& res, size_t encoding, int series_num, const std::string& data_role)
    {
        try
        {
            size_t RO = res.data_.get_size(0);
            size_t E1 = res.data_.get_size(1);
            size_t E2 = res.data_.get_size(2);
            size_t CHA = res.data_.get_size(3);
            size_t N = res.data_.get_size(4);
            size_t S = res.data_.get_size(5);
            size_t SLC = res.data_.get_size(6);

            GDEBUG_CONDITION_STREAM(true, "sending out image array, acquisition boundary [RO E1 E2 CHA N S SLC] = [" << RO << " " << E1 << " " << E2 << " " << CHA << " " << N << " " << S << " " << SLC << "] ");

            // compute image numbers and fill the image meta
            size_t n, s, slc;
            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (n = 0; n < N; n++)
                    {
                        ISMRMRD::ImageHeader header = res.headers_(n, s, slc);

                        if (header.measurement_uid == 0) continue;

                        res.headers_(n, s, slc).image_index = (uint16_t)this->compute_image_number(res.headers_(n, s, slc), encoding, CHA, 0, E2);
                        res.headers_(n, s, slc).image_series_index = series_num;

                        size_t offset = n + s*N + slc*N*S;
                        res.meta_[offset].set(GADGETRON_IMAGENUMBER, (long)res.headers_(n, s, slc).image_index);
                        res.meta_[offset].set(GADGETRON_IMAGEPROCESSINGHISTORY, "GT");

                        if (data_role == GADGETRON_IMAGE_REGULAR)
                        {
                            res.headers_(n, s, slc).image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

                            res.meta_[offset].append(GADGETRON_IMAGECOMMENT, "GT");

                            res.meta_[offset].append(GADGETRON_SEQUENCEDESCRIPTION, "_GT");
                            res.meta_[offset].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_REGULAR);
                        }
                        else if (data_role == GADGETRON_IMAGE_GFACTOR)
                        {
                            res.headers_(n, s, slc).image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

                            res.meta_[offset].append(GADGETRON_IMAGECOMMENT, GADGETRON_IMAGE_GFACTOR);
                            res.meta_[offset].append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_GFACTOR);
                            res.meta_[offset].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_GFACTOR);

                            // set the skip processing flag, so gfactor map will not be processed during e.g. partial fourier handling or kspace filter gadgets
                            res.meta_[offset].set(GADGETRON_SKIP_PROCESSING_AFTER_RECON, (long)1);
                            // set the flag to use dedicated scaling factor
                            res.meta_[offset].set(GADGETRON_USE_DEDICATED_SCALING_FACTOR, (long)1);
                        }

                        if (verbose.value())
                        {
                            for (size_t cha = 0; cha < CHA; cha++)
                            {
                                GDEBUG_STREAM("sending out " << data_role << " image [CHA SLC CON PHS REP SET AVE] = [" << cha << " "<< res.headers_(n, s, slc).slice << " " << res.headers_(n, s, slc).contrast << " "<< res.headers_(n, s, slc).phase << " " << res.headers_(n, s, slc).repetition << " " << res.headers_(n, s, slc).set << " " << res.headers_(n, s, slc).average << " " << "] "<< " -- Image number -- " << res.headers_(n, s, slc).image_index); 
                            }
                        }
                    }
                }
            }

            // send out the images
            Gadgetron::GadgetContainerMessage<IsmrmrdImageArray>* cm1 = new Gadgetron::GadgetContainerMessage<IsmrmrdImageArray>();
            *(cm1->getObjectPtr()) = res;

            if (this->next()->putq(cm1) < 0)
            {
                GERROR_STREAM("Put image array to Q failed ... ");
                return GADGET_FAIL;
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors in MultiChannelCartesianGrappaReconGadget::send_out_image_array(...) ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    // ----------------------------------------------------------------------------------------

    void MultiChannelCartesianGrappaReconGadget::make_ref_coil_map(IsmrmrdDataBuffered& ref_, std::vector<size_t>  recon_dims, ReconObjType& recon_obj, size_t encoding)
    {

        try
        {
            hoNDArray< std::complex<float> >& ref_data = ref_.data_;
            hoNDArray< std::complex<float> >& ref_calib = recon_obj.ref_calib_;
            hoNDArray< std::complex<float> >& ref_coil_map = recon_obj.ref_coil_map_;

            // sampling limits
            size_t sRO = ref_.sampling_.sampling_limits_[0].min_;
            size_t eRO = ref_.sampling_.sampling_limits_[0].max_;
            size_t cRO = ref_.sampling_.sampling_limits_[0].center_;

            size_t sE1 = ref_.sampling_.sampling_limits_[1].min_;
            size_t eE1 = ref_.sampling_.sampling_limits_[1].max_;
            size_t cE1 = ref_.sampling_.sampling_limits_[1].center_;

            size_t sE2 = ref_.sampling_.sampling_limits_[2].min_;
            size_t eE2 = ref_.sampling_.sampling_limits_[2].max_;
            size_t cE2 = ref_.sampling_.sampling_limits_[2].center_;

            // recon size
            size_t recon_RO = recon_dims[0];
            size_t recon_E1 = recon_dims[1];
            size_t recon_E2 = recon_dims[2];

            // ref array size
            size_t CHA = ref_data.get_size(3);
            size_t N = ref_data.get_size(4);
            size_t S = ref_data.get_size(5);
            size_t SLC = ref_data.get_size(6);

            // determine the ref_coil_map size
            size_t RO = 2 * cRO;
            if (sRO>0 || eRO<RO - 1)
            {
                RO = 2 * std::max(cRO - sRO, eRO - cRO+1);
                if (RO>recon_RO) RO = recon_RO;
            }

            size_t E1 = eE1 - sE1 + 1;
            size_t E2 = eE2 - sE2 + 1;

            if ((calib_mode_[encoding] == Gadgetron::ISMRMRD_interleaved) || (calib_mode_[encoding] == Gadgetron::ISMRMRD_noacceleration))
            {
                E1 = 2 * std::max(cE1 - sE1, eE1 - cE1+1);
                if (E1>recon_E1) E1 = recon_E1;

                if (E2 > 1)
                {
                    E2 = 2 * std::max(cE2 - sE2, eE2 - cE2 + 1);
                    if (E2 > recon_E2) E2 = recon_E2;
                }
            }

            ref_coil_map.create(RO, E1, E2, CHA, N, S, SLC);
            Gadgetron::clear(ref_coil_map);

            size_t slc, s, n, cha, e2, e1;
            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (n = 0; n < N; n++)
                    {
                        for (cha = 0; cha < CHA; cha++)
                        {
                            for (e2 = sE2; e2 <= eE2; e2++)
                            {
                                for (e1 = sE1; e1 <= eE1; e1++)
                                {
                                    std::complex<float>* pSrc = &(ref_data(0, e1-sE1, e2-sE2, cha, n, s, slc));
                                    std::complex<float>* pDst = &(ref_coil_map(0, e1, e2, cha, n, s, slc));

                                    memcpy(pDst + sRO, pSrc, sizeof(std::complex<float>)*(eRO - sRO + 1));
                                }
                            }
                        }
                    }
                }
            }

            // filter the ref_coil_map
            if (filter_RO_ref_coi_map_.get_size(0) != RO)
            {
                Gadgetron::generate_symmetric_filter_ref(ref_coil_map.get_size(0), ref_.sampling_.sampling_limits_[0].min_, ref_.sampling_.sampling_limits_[0].max_, filter_RO_ref_coi_map_);
            }

            if (filter_E1_ref_coi_map_.get_size(0) != E1)
            {
                Gadgetron::generate_symmetric_filter_ref(ref_coil_map.get_size(1), ref_.sampling_.sampling_limits_[1].min_, ref_.sampling_.sampling_limits_[1].max_, filter_E1_ref_coi_map_);
            }

            if ( (E2 > 1) && (filter_E2_ref_coi_map_.get_size(0) != E2) )
            {
                Gadgetron::generate_symmetric_filter_ref(ref_coil_map.get_size(2), ref_.sampling_.sampling_limits_[2].min_, ref_.sampling_.sampling_limits_[2].max_, filter_E2_ref_coi_map_);
            }

            hoNDArray< std::complex<float> > ref_recon_buf;

            if (E2 > 1)
            {
                Gadgetron::apply_kspace_filter_ROE1E2(ref_coil_map, filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, filter_E2_ref_coi_map_, ref_recon_buf);
            }
            else
            {
                Gadgetron::apply_kspace_filter_ROE1(ref_coil_map, filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, ref_recon_buf);
            }

            // pad the ref_coil_map into the data array
            Gadgetron::pad(recon_RO, recon_E1, recon_E2, &ref_recon_buf, &ref_coil_map);

            std::vector<size_t> dim = *ref_data.get_dimensions();
            ref_calib.create(dim, ref_data.begin());

        }
        catch (...)
        {
            GADGET_THROW("Errors happened in MultiChannelCartesianGrappaReconGadget::make_ref_coil_map(...) ... ");
        }
    }

    void MultiChannelCartesianGrappaReconGadget::perform_coil_map_estimation(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t e)
    {
        try
        {
            recon_obj.coil_map_ = recon_obj.ref_coil_map_;
            Gadgetron::clear(recon_obj.coil_map_);

            size_t E2 = recon_obj.ref_coil_map_.get_size(2);
            if (E2 > 1)
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft3c(recon_obj.ref_coil_map_, complex_im_recon_buf_);
            }
            else
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft2c(recon_obj.ref_coil_map_, complex_im_recon_buf_);
            }

            size_t ks = 7;
            size_t kz = 5;
            size_t power = 3;

            Gadgetron::coil_map_Inati(complex_im_recon_buf_, recon_obj.coil_map_, ks, kz, power);

        }
        catch (...)
        {
            GADGET_THROW("Errors happened in MultiChannelCartesianGrappaReconGadget::perform_coil_map_estimation(...) ... ");
        }
    }

    void MultiChannelCartesianGrappaReconGadget::perform_calib(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t e)
    {
        try
        {
            size_t RO = recon_bit.data_.data_.get_size(0);
            size_t E1 = recon_bit.data_.data_.get_size(1);
            size_t E2 = recon_bit.data_.data_.get_size(2);

            hoNDArray< std::complex<float> >& src = recon_obj.ref_calib_;
            hoNDArray< std::complex<float> >& dst = recon_obj.ref_calib_;

            size_t ref_RO = src.get_size(0);
            size_t ref_E1 = src.get_size(1);
            size_t ref_E2 = src.get_size(2);
            size_t srcCHA = src.get_size(3);
            size_t ref_N = src.get_size(4);
            size_t ref_S = src.get_size(5);
            size_t ref_SLC = src.get_size(6);

            size_t dstCHA = dst.get_size(3);

            recon_obj.unmixing_coeff_.create(RO, E1, E2, srcCHA, ref_N, ref_S, ref_SLC);
            recon_obj.gfactor_.create(RO, E1, E2, 1, ref_N, ref_S, ref_SLC);

            Gadgetron::clear(recon_obj.unmixing_coeff_);
            Gadgetron::clear(recon_obj.gfactor_);

            if (acceFactorE1_[e] <= 1 && acceFactorE2_[e] <= 1)
            {
                Gadgetron::conjugate(recon_obj.coil_map_, recon_obj.unmixing_coeff_);
            }
            else
            {
                // allocate buffer for kernels
                size_t kRO = grappa_kSize_RO.value();
                size_t kNE1 = grappa_kSize_E1.value();
                size_t kNE2 = grappa_kSize_E2.value();

                size_t convKRO(1), convKE1(1), convKE2(1);

                if (E2 > 1)
                {
                    std::vector<int> kE1, oE1;
                    std::vector<int> kE2, oE2;
                    bool fitItself = true;
                    grappa3d_kerPattern(kE1, oE1, kE2, oE2, convKRO, convKE1, convKE2, (size_t)acceFactorE1_[e], (size_t)acceFactorE2_[e], kRO, kNE1, kNE2, fitItself);
                }
                else
                {
                    std::vector<int> kE1, oE1;
                    bool fitItself = true;

                    Gadgetron::grappa2d_kerPattern(kE1, oE1, convKRO, convKE1, (size_t)acceFactorE1_[e], kRO, kNE1, fitItself);
                    recon_obj.kernelIm_.create(RO, E1, 1, srcCHA, dstCHA, ref_N, ref_S, ref_SLC);
                }

                recon_obj.kernel_.create(convKRO, convKE1, convKE2, srcCHA, dstCHA, ref_N, ref_S, ref_SLC);
                Gadgetron::clear(recon_obj.kernel_);
                Gadgetron::clear(recon_obj.kernelIm_);

                long long num = ref_N*ref_S*ref_SLC;

                long long ii;

#pragma omp parallel for default(none) private(ii) shared(src, dst, recon_obj, e, num, ref_N, ref_S, ref_RO, ref_E1, ref_E2, RO, E1, E2, dstCHA, srcCHA, convKRO, convKE1, convKE2, kRO, kNE1, kNE2) if(num>1)
                for (ii = 0; ii < num; ii++)
                {
                    size_t slc = ii / (ref_N*ref_S);
                    size_t s = (ii - slc*ref_N*ref_S) / (ref_N);
                    size_t n = ii - slc*ref_N*ref_S - s*ref_N;

                    std::stringstream os;
                    os << "n" << n << "_s" << s << "_slc" << slc << "_encoding_" << e;
                    std::string suffix = os.str();

                    std::complex<float>* pSrc = &(src(0, 0, 0, 0, n, s, slc));
                    hoNDArray< std::complex<float> > ref_src(ref_RO, ref_E1, ref_E2, srcCHA, pSrc);

                    std::complex<float>* pDst = &(dst(0, 0, 0, 0, n, s, slc));
                    hoNDArray< std::complex<float> > ref_dst(ref_RO, ref_E1, ref_E2, dstCHA, pDst);

                    // -----------------------------------

                    if (E2 > 1)
                    {
                        hoNDArray< std::complex<float> > ker(convKRO, convKE1, convKE2, srcCHA, dstCHA, &(recon_obj.kernel_(0, 0, 0, 0, 0, n, s, slc)));
                        Gadgetron::grappa3d_calib_convolution_kernel(ref_src, ref_dst, (size_t)acceFactorE1_[e], (size_t)acceFactorE2_[e], grappa_reg_lamda.value(), grappa_calib_over_determine_ratio.value(), kRO, kNE1, kNE2, ker);

                        hoNDArray< std::complex<float> > coilMap(RO, E1, E2, dstCHA, &(recon_obj.coil_map_(0, 0, 0, 0, n, s, slc)));
                        hoNDArray< std::complex<float> > unmixC(RO, E1, E2, srcCHA, &(recon_obj.unmixing_coeff_(0, 0, 0, 0, n, s, slc)));
                        hoNDArray<float> gFactor(RO, E1, E2, 1, &(recon_obj.gfactor_(0, 0, 0, 0, n, s, slc)));
                        Gadgetron::grappa3d_unmixing_coeff(ker, coilMap, (size_t)acceFactorE1_[e], (size_t)acceFactorE2_[e], unmixC, gFactor);

                    }
                    else
                    {
                        hoNDArray< std::complex<float> > acsSrc(ref_RO, ref_E1, srcCHA, const_cast< std::complex<float>*>(ref_src.begin()));
                        hoNDArray< std::complex<float> > acsDst(ref_RO, ref_E1, dstCHA, const_cast< std::complex<float>*>(ref_dst.begin()));

                        hoNDArray< std::complex<float> > convKer(convKRO, convKE1, srcCHA, dstCHA, &(recon_obj.kernel_(0, 0, 0, 0, 0, n, s, slc)));
                        hoNDArray< std::complex<float> > kIm(RO, E1, srcCHA, dstCHA, &(recon_obj.kernelIm_(0, 0, 0, 0, 0, n, s, slc)));

                        Gadgetron::grappa2d_calib_convolution_kernel(acsSrc, acsDst, (size_t)acceFactorE1_[e], grappa_reg_lamda.value(), kRO, kNE1, convKer);
                        Gadgetron::grappa2d_image_domain_kernel(convKer, RO, E1, kIm);

                        hoNDArray< std::complex<float> > coilMap(RO, E1, dstCHA, &(recon_obj.coil_map_(0, 0, 0, 0, n, s, slc)));
                        hoNDArray< std::complex<float> > unmixC(RO, E1, srcCHA, &(recon_obj.unmixing_coeff_(0, 0, 0, 0, n, s, slc)));
                        hoNDArray<float> gFactor;

                        Gadgetron::grappa2d_unmixing_coeff(kIm, coilMap, (size_t)acceFactorE1_[e], unmixC, gFactor);
                        memcpy(&(recon_obj.gfactor_(0, 0, 0, 0, n, s, slc)), gFactor.begin(), gFactor.get_number_of_bytes());
                    }

                    // -----------------------------------
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in MultiChannelCartesianGrappaReconGadget::perform_calib(...) ... ");
        }
    }

    void MultiChannelCartesianGrappaReconGadget::perform_unwrapping(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t e)
    {
        try
        {
            typedef std::complex<float> T;

            size_t RO = recon_bit.data_.data_.get_size(0);
            size_t E1 = recon_bit.data_.data_.get_size(1);
            size_t E2 = recon_bit.data_.data_.get_size(2);
            size_t dstCHA = recon_bit.data_.data_.get_size(3);
            size_t N = recon_bit.data_.data_.get_size(4);
            size_t S = recon_bit.data_.data_.get_size(5);
            size_t SLC = recon_bit.data_.data_.get_size(6);

            hoNDArray< std::complex<float> >& src = recon_obj.ref_calib_;
            hoNDArray< std::complex<float> >& dst = recon_obj.ref_calib_;

            size_t ref_RO = src.get_size(0);
            size_t ref_E1 = src.get_size(1);
            size_t ref_E2 = src.get_size(2);
            size_t srcCHA = src.get_size(3);
            size_t ref_N = src.get_size(4);
            size_t ref_S = src.get_size(5);
            size_t ref_SLC = src.get_size(6);

            size_t convkRO = recon_obj.kernel_.get_size(0);
            size_t convkE1 = recon_obj.kernel_.get_size(1);
            size_t convkE2 = recon_obj.kernel_.get_size(2);

            recon_obj.recon_res_.data_.create(RO, E1, E2, dstCHA, N, S, SLC);

            // compute aliased images
            data_recon_buf_.create(RO, E1, E2, dstCHA, N, S, SLC);

            if (E2>1)
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft3c(recon_bit.data_.data_, complex_im_recon_buf_, data_recon_buf_);
            }
            else
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft2c(recon_bit.data_.data_, complex_im_recon_buf_, data_recon_buf_);
            }

            // SNR unit scaling
            float effectiveAcceFactor = acceFactorE1_[e] * acceFactorE2_[e];
            if (effectiveAcceFactor > 1)
            {
                float fftCompensationRatio = (float)(1.0 / std::sqrt(effectiveAcceFactor));
                Gadgetron::scal(fftCompensationRatio, complex_im_recon_buf_);
            }
            // unwrapping

            long long num = N*S*SLC;

            long long ii;

#pragma omp parallel default(none) private(ii) shared(num, N, S, RO, E1, E2, srcCHA, convkRO, convkE1, convkE2, ref_N, ref_S, recon_obj, dstCHA, e) if(num>1)
            {
#pragma omp for 
                for (ii = 0; ii < num; ii++)
                {
                    size_t slc = ii / (N*S);
                    size_t s = (ii - slc*N*S) / N;
                    size_t n = ii - slc*N*S - s*N;

                    // combined channels
                    T* pIm = &(complex_im_recon_buf_(0, 0, 0, 0, n, s, slc));
                    hoNDArray< std::complex<float> > aliasedIm(RO, E1, E2, srcCHA, 1, pIm);

                    size_t usedN = n;
                    if (n >= ref_N) usedN = ref_N - 1;

                    size_t usedS = s;
                    if (s >= ref_S) usedS = ref_S - 1;

                    T* pUnmix = &(recon_obj.unmixing_coeff_(0, 0, 0, 0, usedN, usedS, slc));
                    hoNDArray< std::complex<float> > unmixing(RO, E1, E2, srcCHA, pUnmix);

                    T* pRes = &(recon_obj.recon_res_.data_(0, 0, 0, 0, n, s, slc));
                    hoNDArray< std::complex<float> > res(RO, E1, E2, dstCHA, pRes);
		
                    Gadgetron::apply_unmix_coeff_aliased_image_3D(aliasedIm, unmixing, res);
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in MultiChannelCartesianGrappaReconGadget::perform_unwrapping(...) ... ");
        }
    }

    void MultiChannelCartesianGrappaReconGadget::compute_image_header(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t e)
    {
        try
        {
            size_t RO = recon_obj.recon_res_.data_.get_size(0);
            size_t E1 = recon_obj.recon_res_.data_.get_size(1);
            size_t E2 = recon_obj.recon_res_.data_.get_size(2);
            size_t CHA = recon_obj.recon_res_.data_.get_size(3);
            size_t N = recon_obj.recon_res_.data_.get_size(4);
            size_t S = recon_obj.recon_res_.data_.get_size(5);
            size_t SLC = recon_obj.recon_res_.data_.get_size(6);

            GADGET_CHECK_THROW(N == recon_bit.data_.headers_.get_size(2));
            GADGET_CHECK_THROW(S == recon_bit.data_.headers_.get_size(3));

            recon_obj.recon_res_.headers_.create(N, S, SLC);
            recon_obj.recon_res_.meta_.resize(N*S*SLC);

            size_t n, s, slc;

            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (n = 0; n < N; n++)
                    {
                        size_t header_E1 = recon_bit.data_.headers_.get_size(0);
                        size_t header_E2 = recon_bit.data_.headers_.get_size(1);

                        // for every kspace, find the recorded header which is closest to the kspace center [E1/2 E2/2]
                        ISMRMRD::AcquisitionHeader acq_header;

                        long long bestE1 = E1 + 1;
                        long long bestE2 = E2 + 1;

                        size_t e1, e2;
                        for (e2 = 0; e2 < header_E2; e2++)
                        {
                            for (e1 = 0; e1 < header_E1; e1++)
                            {
                                ISMRMRD::AcquisitionHeader& curr_header = recon_bit.data_.headers_(e1, e2, n, s, slc);

                                // if (curr_header.measurement_uid != 0) // a valid header
                                {
                                    if (E2 > 1)
                                    {
                                        if (std::abs((long long)curr_header.idx.kspace_encode_step_1 - (long long)(E1 / 2)) < bestE1
                                            && std::abs((long long)curr_header.idx.kspace_encode_step_2 - (long long)(E2 / 2)) < bestE2)
                                        {
                                            bestE1 = std::abs((long long)curr_header.idx.kspace_encode_step_1 - (long long)E1 / 2);
                                            bestE2 = std::abs((long long)curr_header.idx.kspace_encode_step_2 - (long long)E2 / 2);

                                            acq_header = curr_header;
                                        }
                                    }
                                    else
                                    {
                                        if (std::abs((long long)curr_header.idx.kspace_encode_step_1 - (long long)(E1 / 2)) < bestE1)
                                        {
                                            bestE1 = std::abs((long long)curr_header.idx.kspace_encode_step_1 - (long long)E1 / 2);

                                            acq_header = curr_header;
                                        }
                                    }
                                }
                            }
                        }

                        //if (acq_header.measurement_uid == 0)
                        //{
                        //    std::ostringstream ostr;
                        //    ostr << "Cannot create valid image header : n = " << n << ", s = " << s << ", slc = " << slc;
                        //    GADGET_THROW(ostr.str());
                        //}
                        //else
                        //{
                            ISMRMRD::ImageHeader& im_header = recon_obj.recon_res_.headers_(n, s, slc);
                            ISMRMRD::MetaContainer& meta = recon_obj.recon_res_.meta_[n + s*N + slc*N*S];

                            im_header.version = acq_header.version;
                            im_header.data_type = ISMRMRD::ISMRMRD_CXFLOAT;
                            im_header.flags = acq_header.flags;
                            im_header.measurement_uid = acq_header.measurement_uid;

                            im_header.matrix_size[0] = (uint16_t)RO;
                            im_header.matrix_size[1] = (uint16_t)E1;
                            im_header.matrix_size[2] = (uint16_t)E2;

                            im_header.field_of_view[0] = recon_bit.data_.sampling_.recon_FOV_[0];
                            im_header.field_of_view[1] = recon_bit.data_.sampling_.recon_FOV_[1];
                            im_header.field_of_view[2] = recon_bit.data_.sampling_.recon_FOV_[2];

                            im_header.channels = (uint16_t)CHA;

                            im_header.position[0] = acq_header.position[0];
                            im_header.position[1] = acq_header.position[1];
                            im_header.position[2] = acq_header.position[2];

                            im_header.read_dir[0] = acq_header.read_dir[0];
                            im_header.read_dir[1] = acq_header.read_dir[1];
                            im_header.read_dir[2] = acq_header.read_dir[2];

                            im_header.phase_dir[0] = acq_header.phase_dir[0];
                            im_header.phase_dir[1] = acq_header.phase_dir[1];
                            im_header.phase_dir[2] = acq_header.phase_dir[2];

                            im_header.slice_dir[0] = acq_header.slice_dir[0];
                            im_header.slice_dir[1] = acq_header.slice_dir[1];
                            im_header.slice_dir[2] = acq_header.slice_dir[2];

                            im_header.patient_table_position[0] = acq_header.patient_table_position[0];
                            im_header.patient_table_position[1] = acq_header.patient_table_position[1];
                            im_header.patient_table_position[2] = acq_header.patient_table_position[2];

                            im_header.average = acq_header.idx.average;
                            im_header.slice = acq_header.idx.slice;
                            im_header.contrast = acq_header.idx.contrast;
                            im_header.phase = acq_header.idx.phase;
                            im_header.repetition = acq_header.idx.repetition;
                            im_header.set = acq_header.idx.set;

                            im_header.acquisition_time_stamp = acq_header.acquisition_time_stamp;

                            im_header.physiology_time_stamp[0] = acq_header.physiology_time_stamp[0];
                            im_header.physiology_time_stamp[1] = acq_header.physiology_time_stamp[1];
                            im_header.physiology_time_stamp[2] = acq_header.physiology_time_stamp[2];

                            im_header.image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;
                            im_header.image_index = (uint16_t)(n + s*N + slc*N*S);
                            im_header.image_series_index = 0;

                            memcpy(im_header.user_int, acq_header.user_int, sizeof(int32_t)*ISMRMRD::ISMRMRD_USER_INTS);
                            memcpy(im_header.user_float, acq_header.user_float, sizeof(float)*ISMRMRD::ISMRMRD_USER_FLOATS);

                            im_header.attribute_string_len = 0;

                            meta.set("encoding", (long)e);

                            meta.set("encoding_FOV", recon_bit.data_.sampling_.encoded_FOV_[0]);
                            meta.append("encoding_FOV", recon_bit.data_.sampling_.encoded_FOV_[1]);
                            meta.append("encoding_FOV", recon_bit.data_.sampling_.encoded_FOV_[2]);

                            meta.set("recon_FOV", recon_bit.data_.sampling_.recon_FOV_[0]);
                            meta.append("recon_FOV", recon_bit.data_.sampling_.recon_FOV_[1]);
                            meta.append("recon_FOV", recon_bit.data_.sampling_.recon_FOV_[2]);

                            meta.set("encoded_matrix", (long)recon_bit.data_.sampling_.encoded_matrix_[0]);
                            meta.append("encoded_matrix", (long)recon_bit.data_.sampling_.encoded_matrix_[1]);
                            meta.append("encoded_matrix", (long)recon_bit.data_.sampling_.encoded_matrix_[2]);

                            meta.set("recon_matrix", (long)recon_bit.data_.sampling_.recon_matrix_[0]);
                            meta.append("recon_matrix", (long)recon_bit.data_.sampling_.recon_matrix_[1]);
                            meta.append("recon_matrix", (long)recon_bit.data_.sampling_.recon_matrix_[2]);

                            meta.set("sampling_limits_RO", (long)recon_bit.data_.sampling_.sampling_limits_[0].min_);
                            meta.append("sampling_limits_RO", (long)recon_bit.data_.sampling_.sampling_limits_[0].center_);
                            meta.append("sampling_limits_RO", (long)recon_bit.data_.sampling_.sampling_limits_[0].max_);

                            meta.set("sampling_limits_E1", (long)recon_bit.data_.sampling_.sampling_limits_[1].min_);
                            meta.append("sampling_limits_E1", (long)recon_bit.data_.sampling_.sampling_limits_[1].center_);
                            meta.append("sampling_limits_E1", (long)recon_bit.data_.sampling_.sampling_limits_[1].max_);

                            meta.set("sampling_limits_E2", (long)recon_bit.data_.sampling_.sampling_limits_[2].min_);
                            meta.append("sampling_limits_E2", (long)recon_bit.data_.sampling_.sampling_limits_[2].center_);
                            meta.append("sampling_limits_E2", (long)recon_bit.data_.sampling_.sampling_limits_[2].max_);
                        //}
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in MultiChannelCartesianGrappaReconGadget::compute_image_header(...) ... ");
        }
    }

    GADGET_FACTORY_DECLARE(MultiChannelCartesianGrappaReconGadget)
}
