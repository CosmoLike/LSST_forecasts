#!/usr/bin/env python

import sys
#sys.path.append('/Users/teifler/Dropbox/cosmolike/top-level/miyatake/')

from cosmolike_libs_onlyMG import * 

file_source_z = os.path.join(dirname, "../../zdistris/zdistribution_LSST")
file_lens_z = os.path.join(dirname, "../../zdistris/zdistribution_LSST_gold")
data_file = os.path.join(dirname, "datav/LSST_all_2pt_fid")
cov_file = os.path.join(dirname, "cov/Multi_Probe_LSST_2.600000e+01_1.800000e+04_Rmin10_lmin30_lens-zmin0.15LSST_cov_Ncl20_Ntomo10_2pt_inv")
chain_file = os.path.join(dirname, "./like/like_LSST_all_2pt_noredshifterr_MGonly_zmin0.15")

update_Ntable_Na(20)
initcosmo()
initbins(20,30.0,1000.0,1000.0,10.0,7)
initsurvey("LSST")
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","gold")
initclusters()
initia("none","DEEP2")
initpriors("none","none","none","none")
#initprobes("shear_shear")
initprobes("all_2pt")
initdatainv(cov_file ,data_file)
initclusters()
update_nuisance_LSST_source()

sample_params= sample_cosmology_only(MGonly = True)
#sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear())
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance_IA_marg(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering()) 
sample_main(sample_params,4000,64,24,chain_file, blind=False)
