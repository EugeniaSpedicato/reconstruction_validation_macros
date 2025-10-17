rule all:
	input:
		expand('/home/espedica/modifica_fair_install/modificaFairRoot/share/MUonE/common/jobs/snakemake/try/only_reco_myTRMesmer_{number_params}_{n_hits}hit.yaml',n_hits=['0','1','2'],number_params=['0_5','5_10','10_15','15_20','20_25','25_32']),
		expand('/home/espedica/macros_fairmu/snakemake/plots/umb_mu_bv_{n_hits}hit_index{number_params}_{info}.root',n_hits=['0','1','2'],number_params=['0_5','5_10','10_15','15_20','20_25','25_32'],info=['defaultModificaRealGeom']),
		expand('/home/espedica/macros_fairmu/snakemake/plots/hits_shared_{n_hits}hit_{number_params}_{info}.pdf',n_hits=['0','1','2'],number_params=['0_5','5_10','10_15','15_20','20_25','25_32'],info=['defaultModificaRealGeom']),
		expand('/home/espedica/macros_fairmu/snakemake/plots/reco_el_{n_hits}hit_LO_{number_params}_{info}.root',n_hits=['0','1','2'],number_params=['0_5','5_10','10_15','15_20','20_25','25_32'],info=['defaultModificaRealGeom'])
		#expand('/home/espedica/macros_fairmu/snakemake/plots/h_posZ_{n_hits}hit_{number_params}_{info}.pdf',n_hits=['0','1','2'],number_params=['0_5','5_10','10_15','15_20','20_25','25_32'],info=['defaultModificaRealGeom'])

rule rm_yaml:
	input:
		path_y='/home/espedica/modifica_fair_install/modificaFairRoot/share/MUonE/common/jobs/snakemake/try'
	shell:
		"""rm -rf {input.path_y}/*yaml"""

rule rm_files:
	input:
		path_f='/mnt/raid10/DATA/espedica/fairmu/reco/snakemake'
	shell:
		"""rm -rf {input.path_f}/*defaultModificaRealGeom.root"""

rule rm_plots:
	input:
		path_p="/home/espedica/macros_fairmu/snakemake/plots"
	shell:
		"""rm -rf {input.path_p}/*defaultModificaRealGeom.pdf {input.path_p}/*defaultModificaRealGeom.root {input.path_p}/results/*defaultModificaRealGeom.pdf {input.path_p}/results/*defaultModificaRealGeom.root"""


rule yaml_files:
	input:
		path='/home/espedica/modifica_fair_install/modificaFairRoot/share/MUonE/common/jobs/snakemake/try'
	params:
		nevents='100000',
		info_y='defaultModificaRealGeom'
	output:
		outfile='/home/espedica/modifica_fair_install/modificaFairRoot/share/MUonE/common/jobs/snakemake/try/only_reco_myTRMesmer_{number_params}_{n_hits}hit.yaml'
	shell:
		"""source {input.path}/yaml_creator.sh {input.path}/only_reco_myTRMesmer_{wildcards.number_params}_{wildcards.n_hits}hit.yaml {wildcards.number_params} {wildcards.n_hits} {params.info_y} {params.nevents}"""

rule reco:
	input:
		name='/home/espedica/modifica_fair_install/modificaFairRoot/share/MUonE/common/jobs/snakemake/try/only_reco_myTRMesmer_{number_params}_{n_hits}hit.yaml'
	output:
		out='/mnt/raid10/DATA/espedica/fairmu/reco/snakemake/range_{number_params}_{n_hits}hit_{info}.root'
	shell:
		"""root -q 'runProductionJob.C("{input.name}")'"""

rule validation_resolution:
	input:
		reconstructed_r='/mnt/raid10/DATA/espedica/fairmu/reco/snakemake/range_{number_params}_{n_hits}hit_{info}.root',
		#generated_r='/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_258faf6b_{number_params}mrad.root',
		#generated_r='/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_7e6f0f25_{number_params}mrad.root',
		#generated_r='/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_d652278e_{number_params}mrad_notilt.root',
		generated_r='/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_d652278e_{number_params}mrad_realGeom.root',
		path_r='/home/espedica/macros_fairmu/snakemake/plots'
	output:
		out_r='/home/espedica/macros_fairmu/snakemake/plots/umb_mu_bv_{n_hits}hit_index{number_params}_{info}.root'
	shell:
		"""root -q '/home/espedica/macros_fairmu/snakemake/resolutions.cpp("{input.reconstructed_r}","{input.generated_r}",{wildcards.n_hits},"{wildcards.number_params}","{wildcards.info}","{input.path_r}")'"""

rule validation_efficiency:
	input:
		reconstructed_e='/mnt/raid10/DATA/espedica/fairmu/reco/snakemake/range_{number_params}_{n_hits}hit_{info}.root',
		#generated_e='/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_258faf6b_{number_params}mrad.root',
		#generated_e='/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_7e6f0f25_{number_params}mrad.root',
		#generated_e='/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_d652278e_{number_params}mrad_notilt.root',
		generated_e='/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_d652278e_{number_params}mrad_realGeom.root',
		path_e='/home/espedica/macros_fairmu/snakemake/plots'
	output:
		out_e='/home/espedica/macros_fairmu/snakemake/plots/hits_shared_{n_hits}hit_{number_params}_{info}.pdf'		
	shell:
		"""root -q '/home/espedica/macros_fairmu/snakemake/efficiency_MC.cpp("{input.reconstructed_e}","{input.generated_e}",{wildcards.n_hits},"{wildcards.number_params}","{wildcards.info}","{input.path_e}")'"""

rule validation_fake:
	input:
		reconstructed_fk='/mnt/raid10/DATA/espedica/fairmu/reco/snakemake/range_{number_params}_{n_hits}hit_{info}.root',
		#generated_fk='/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_258faf6b_{number_params}mrad.root',
		#generated_fk='/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_7e6f0f25_{number_params}mrad.root',
		#generated_fk='/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_d652278e_{number_params}mrad_notilt.root',
		generated_fk='/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_d652278e_{number_params}mrad_realGeom.root',

		path_fk='/home/espedica/macros_fairmu/snakemake/plots'
	output:
		out_fk='/home/espedica/macros_fairmu/snakemake/plots/reco_el_{n_hits}hit_LO_{number_params}_{info}.root'
	shell:
		"""root -q '/home/espedica/macros_fairmu/snakemake/fake_rates_and_quality_prevrtx.cpp("{input.reconstructed_fk}","{input.generated_fk}",{wildcards.n_hits},"{wildcards.number_params}","{wildcards.info}","{input.path_fk}")'"""

#rule vrtx:
#	input:
#		reconstructed_e='/mnt/raid10/DATA/espedica/fairmu/reco/snakemake/range_{number_params}_{n_hits}hit_{info}.root',
#		generated_e='/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_2f4e96f4_{number_params}mrad.root'
#	output:
#		out_e='/home/espedica/macros_fairmu/snakemake/plots/h_posZ_{n_hits}hit_{number_params}_{info}.pdf'
#	shell:
#		"""root -q '/home/espedica/macros_fairmu/snakemake/vrtx.cpp("{input.reconstructed_e}","{input.generated_e}",{wildcards.n_hits},"{wildcards.number_params}","{wildcards.info}")'"""

rule check_res:
	input:
		path_p='/home/espedica/macros_fairmu/snakemake/plots'
	output:
		out_p='/home/espedica/macros_fairmu/snakemake/plots/results/sigma_log_prepost_defaultModificaRealGeom.pdf'
	shell:
		"""root -q '/home/espedica/macros_fairmu/snakemake/plots_fit_resolutions.cpp("{input.path_p}/","defaultModificaRealGeom")'"""

rule check_eff:
	input:
		path_h='/home/espedica/macros_fairmu/snakemake/plots'
	output:
		out_h='/home/espedica/macros_fairmu/snakemake/plots/results/eff_event_LO_012hit_defaultModificaRealGeom.pdf'
	shell:
		"""root -q '/home/espedica/macros_fairmu/snakemake/plots_012hits_eff.cpp("{input.path_h}/","defaultModificaRealGeom")'"""

rule check_fake:
	input:
		path_h='/home/espedica/macros_fairmu/snakemake/plots'
	output:
		out_h='/home/espedica/macros_fairmu/snakemake/plots/results/good_fake_rate_LO_012hit_defaultModificaRealGeom.pdf'
	shell:
		"""root -q '/home/espedica/macros_fairmu/snakemake/plots_fake.cpp("{input.path_h}/","defaultModificaRealGeom")'"""

rule check_res_pre_post0:
	input:
		path_h='/home/espedica/macros_fairmu/snakemake/plots'
	output:
		out_h='/home/espedica/macros_fairmu/snakemake/plots/results/res_mu_pre_0hit_defaultModificaRealGeom.pdf'
	shell:
		"""root -q '/home/espedica/macros_fairmu/snakemake/plots_pre_post_vrtx_0hit.cpp("{input.path_h}/","defaultModificaRealGeom")'"""

rule comparison0:
	input:
		path1='/home/espedica/macros_fairmu/snakemake/plots/results',
		path2='/home/espedica/macros_fairmu/snakemake/plots/results'
	output:
		out_c='/home/espedica/macros_fairmu/snakemake/plots/results/sigma_overlap_0hit_defaultModificaRealGeom.pdf'
	shell:
		"""root -q '/home/espedica/macros_fairmu/snakemake/compare_version_0hit.cpp("{input.path1}/","{input.path2}/","defaultModificaRealGeom","defaultModificaRealGeom")'"""

rule check_res_pre_post1:
	input:
		path_h='/home/espedica/macros_fairmu/snakemake/plots'
	output:
		out_h='/home/espedica/macros_fairmu/snakemake/plots/results/res_mu_pre_1hit_defaultModificaRealGeom.pdf'
	shell:
		"""root -q '/home/espedica/macros_fairmu/snakemake/plots_pre_post_vrtx_1hit.cpp("{input.path_h}/","defaultModificaRealGeom")'"""

rule comparison1:
	input:
		path1='/home/espedica/macros_fairmu/snakemake/plots/results',
		path2='/home/espedica/macros_fairmu/snakemake/plots/results'
	output:
		out_c='/home/espedica/macros_fairmu/snakemake/plots/results/sigma_overlap_1hit_defaultModificaRealGeom.pdf'
	shell:
		"""root -q '/home/espedica/macros_fairmu/snakemake/compare_version_1hit.cpp("{input.path1}/","{input.path2}/","defaultModificaRealGeom","defaultModificaRealGeom")'"""

rule check_res_pre_post2:
	input:
		path_h='/home/espedica/macros_fairmu/snakemake/plots'
	output:
		out_h='/home/espedica/macros_fairmu/snakemake/plots/results/res_mu_pre_2hit_defaultModificaRealGeom.pdf'
	shell:
		"""root -q '/home/espedica/macros_fairmu/snakemake/plots_pre_post_vrtx_2hit.cpp("{input.path_h}/","defaultModificaRealGeom")'"""

rule comparison2:
	input:
		path1='/home/espedica/macros_fairmu/snakemake/plots/results',
		path2='/home/espedica/macros_fairmu/snakemake/plots/results'
	output:
		out_c='/home/espedica/macros_fairmu/snakemake/plots/results/sigma_overlap_2hit_defaultModificaRealGeom.pdf'
	shell:
		"""root -q '/home/espedica/macros_fairmu/snakemake/compare_version_2hit.cpp("{input.path1}/","{input.path2}/","defaultModificaRealGeom","defaultModificaRealGeom")'"""


