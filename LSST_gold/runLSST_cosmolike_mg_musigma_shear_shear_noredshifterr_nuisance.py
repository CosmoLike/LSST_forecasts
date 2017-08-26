#!/usr/bin/env python

import sys
#sys.path.append('/Users/teifler/Dropbox/cosmolike/top-level/miyatake/')

from cosmolike_libs import * 

file_source_z = os.path.join(dirname, "../../zdistris/zdistribution_LSST")
file_lens_z = os.path.join(dirname, "../../zdistris/zdistribution_LSST_gold")
data_file = os.path.join(dirname, "datav/LSST_shear_shear_fid")
cov_file = os.path.join(dirname, "cov/Multi_Probe_LSST_2.600000e+01_1.800000e+04_Rmin10LSST_cov_Ncl25_Ntomo10_noredshifterr_shear_inv")
chain_file = os.path.join(dirname, "./like/like_LSST_shear_shear_noredshifterr_nuisance")

update_Ntable_Na(20)
initcosmo()
initbins(25,20.0,1000.0,1000.0,10.0,7)
initsurvey("LSST")
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","gold")
initclusters()
initia("none","DEEP2")
initpriors("none","none","none","none")
initprobes("shear_shear")
initdatainv(cov_file ,data_file)
initclusters()
update_nuisance_LSST_source()

#sample_params= sample_cosmology_only(MG = True)
sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear(), MG=True)
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance_IA_marg(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering()) 
sample_main(sample_params,4000,64,24,chain_file, blind=False)
