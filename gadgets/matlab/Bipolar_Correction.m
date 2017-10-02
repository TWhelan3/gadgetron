classdef Bipolar_Correction < handle & BaseGadget
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
    end

    methods

        function g = config(g)
        end

        function g = process(g, head, data)
            addpath('/home/tim/Desktop/B0_NICE/Common_code')
            addpath('/home/tim/Desktop/B0_NICE/mfile-specific')
            addpath('/home/tim/Desktop/B0_NICE/PUROR2D_MatlabCode/')
            addpath('/home/tim/Desktop/B0_NICE/PUROR3D_MatlabCode/')
            addpath('/home/tim/Desktop/B0_NICE/PolyfitnTools/PolyfitnTools')
            tic()
            matrix_size = size(data)
            complex_image = data;
            total_ch = matrix_size(4);
            for index_ch = 1:total_ch
                if index_ch == 1
                    phase_diff(:,:,:) = angle((complex_image(:,:,:,index_ch,2).^2).*conj(complex_image(:,:,:,index_ch,1).*complex_image(:,:,:,index_ch,3)));
                    complex_diff(:,:,:) = abs(complex_image(:,:,:,index_ch,2)).*exp(1i*phase_diff);%complex_31.*conj(complex_21);
                else
                    phase_diff(:,:,:) = angle((complex_image(:,:,:,index_ch,2).^2).*conj(complex_image(:,:,:,index_ch,1).*complex_image(:,:,:,index_ch,3)));
                    complex_diff = complex_diff + abs(complex_image(:,:,:,index_ch,2)).*exp(1i*phase_diff);
                end
            %
            end
        %
            complex_diff(isnan(complex_diff)) = 1;
            [LinearPhaseError_3D] = estimate_bpPhaseError(complex_diff);
            LinearPhaseError_3D = LinearPhaseError_3D./2;
            %
            for index_ch = 1:total_ch
                %---------------------------------------
                for index_echo = 1:matrix_size(5)
                    complex_3D(:,:,:) = complex_image(:,:,:,index_ch,index_echo);
                    if mod(index_echo,2) == 1
                    complex_3D = complex_3D.*exp(1i*(LinearPhaseError_3D./2));
                    else
                    complex_3D = complex_3D.*exp(-1i*(LinearPhaseError_3D./2));
                    end
                    complex_image(:,:,:,index_ch,index_echo) = complex_3D;
                end
            end
            %---------------------------------------
            %complex_tmp(:,:,:,1,:) = complex_image;
            %complex_tmp(:,:,:,1,:) = mag_SoS.*exp(1i*angle(complex_image));
            data=complex_image;
            g.putQ(head,data);
            toc()
        end
    end
end
