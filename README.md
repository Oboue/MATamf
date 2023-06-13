# MATamf
MATamf is a Matlab package for  the advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets. The MaATamf package has a variety of applications exploration and earthquake seismology.


Copyright (C) Oboue et al., 2022

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details: http://www.gnu.org/licenses/


test_amf_synlcre.m This script demonstrates the denoising performance on 2D synthetic seismic data containing linear and curved events corrupted by random and strong erratic noise. 

test_amf_rs.m This script demonstrates the denoising performance on pre-stack and post-stack 2D real reflection seismic data contaminated by erratic and random noise. 

test_amf_das1.m This script demonstrates the denoising performance on 2D raw DAS seismic data corrupted by a mixture of strong noise. 

test_amf_das2.m This script demonstrates the adaptability of the proposed method for weak and strong DAS seismic signal denoising. 

test_amf_rf.m This script shows the application of the AMF method in improving the receiver function imaging.

test_amf_ssp.m This script shows the application of the AMF method in enhancing the SS precursor signals arising from the earth’s mantle transition zone discontinuities.

test_amf_sosvmf_somf_mf.m This script is used to conduct a comparison of the denoising performance between the MF, SOMF, SOSVMF, and the proposed AMF methods.

test_amf_bp_sosvmf_ct_drr_syn.m This script demonstrates the denoising performance of different methods on 2D synthetic containing linear and curved events, corrupted by random and strong erratic noise. 

test_amf_bp_sosvmf_ct_drr_rs.m This script demonstrates the denoising performance of different methods on 2D raw reflection seismic data.

test_amf_bp_sosvmf_ct_drr_rf.m This script demonstrates the denoising performance of different methods on 2D receiver function data.

test_amf_bp_sosvmf_ct_drr_ssp.m This script demonstrates the denoising performance of different methods on 2D SS precursor data.

test_amf_bp_sosvmf_fk_ct_drr_das.m This script demonstrates the denoising performance of different methods on 2D raw data recorded on DAS.
























References

Oboue et al., 2022

Wang, H., Chen, Y., Saad, O.M., Chen, W., Oboué, Y.A.S.I., Yang, L., Fomel, S. and Chen, Y., 2022. A Matlab code package for 2D/3D local slope estimation and structural filtering. Geophysics, 87(3), pp.F1–F14.

Huang, G., M. Bai, Q. Zhao, W. Chen, and Y. Chen, 2021, Erratic noise suppression using iterative structure-oriented space-varying median filtering with sparsity constraint, Geophysical Prospecting, 69, 101-121.

Chen, Y., S. Zu, Y. Wang, and X. Chen, 2020, Deblending of simultaneous-source data using a structure-oriented space varying median filter, Geophysical Journal International, 222, 1805-1723.

Zhao, Q., Q. Du, X. Gong, and Y. Chen, 2018, Signal-preserving erratic noise attenuation via iterative robust sparsity-promoting filter, IEEE Transactions on Geoscience and Remote Sensing, 56, 1558-0644.

For any questions regarding the package, please contact Yangkang Chen (chenyk2016@gmail.com) or Innocent Oboue (obouesonofgod1@gmail.com).
