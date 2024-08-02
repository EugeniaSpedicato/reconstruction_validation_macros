rule all:
	input:
		expand('/home/espedica/macros_fairmu/snakemake/plots/h_res_pmuin_{n_hits}hit_index{number_params}_{info}.pdf',n_hits=['2'],number_params=['0_5','5_10','10_15','15_20','20_25','25_32'],info=['first2fittedPosition']),
		expand('/home/espedica/macros_fairmu/snakemake/plots/hits_shared_{n_hits}hit_{number_params}_{info}.pdf',n_hits=['2'],number_params=['0_5','5_10','10_15','15_20','20_25','25_32'],info=['first2fittedPosition'])

rule yaml_files:
	input:
		path='/home/espedica/fair_install/instFairRoot/share/MUonE/common/jobs/snakemake/try'
	params:
		nevents='100000',
		info_y='first2fittedPosition',
	output:
		outfile='/home/espedica/fair_install/instFairRoot/share/MUonE/common/jobs/snakemake/try/only_reco_myTRMesmer_{number_params}_{n_hits}hit.yaml'
	shell:
		"""source {input.path}/prova.sh {input.path}/only_reco_myTRMesmer_{wildcards.number_params}_{wildcards.n_hits}hit.yaml {wildcards.number_params} {wildcards.n_hits} {params.info_y} {params.nevents}"""

rule reco:
	input:
		name='/home/espedica/fair_install/instFairRoot/share/MUonE/common/jobs/snakemake/try/only_reco_myTRMesmer_{number_params}_{n_hits}hit.yaml'
	output:
		out='/mnt/raid10/DATA/espedica/fairmu/reco/snakemake/range_{number_params}_{n_hits}hit_{info}.root'
	shell:
		"""root -q 'runProductionJob.C("{input.name}")'"""

rule validation_resolution:
	input:
		reconstructed_r='/mnt/raid10/DATA/espedica/fairmu/reco/snakemake/range_{number_params}_{n_hits}hit_{info}.root',
		generated_r='/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_258faf6b_{number_params}mrad.root'
	output:
		out_r='/home/espedica/macros_fairmu/snakemake/plots/h_res_pmuin_{n_hits}hit_index{number_params}_{info}.pdf'
	shell:
		"""root -q '/home/espedica/macros_fairmu/snakemake/resolutions.cpp("{input.reconstructed_r}","{input.generated_r}",{wildcards.n_hits},"{wildcards.number_params}","{wildcards.info}")'"""

rule validation_efficiency:
	input:
		reconstructed_e='/mnt/raid10/DATA/espedica/fairmu/reco/snakemake/range_{number_params}_{n_hits}hit_{info}.root',
		generated_e='/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_258faf6b_{number_params}mrad.root'
	output:
		out_e='/home/espedica/macros_fairmu/snakemake/plots/hits_shared_{n_hits}hit_{number_params}_{info}.pdf'		
	shell:
		"""root -q '/home/espedica/macros_fairmu/snakemake/efficiency_MC.cpp("{input.reconstructed_e}","{input.generated_e}",{wildcards.n_hits},"{wildcards.number_params}","{wildcards.info}")'"""

rule check_res:
	input:
		path_p='/home/espedica/macros_fairmu/snakemake/plots'
	output:
		out_p='/home/espedica/macros_fairmu/snakemake/plots/results/sigma_log_prepost_first2fittedPosition.pdf'
#		out_p='/home/espedica/macros_fairmu/snakemake/plots/results/sigma_prepost_{info}.pdf'
	shell:
		"""root -q '/home/espedica/macros_fairmu/snakemake/plots_fit_resolutions.cpp("{input.path_p}/","first2fittedPosition")'"""
#plots_pre_post_vrtx.cpp("{input.path_p}/")'"""

rule check_eff:
	input:
		path_h='/home/espedica/macros_fairmu/snakemake/plots'
	output:
		out_h='/home/espedica/macros_fairmu/snakemake/plots/results/eff_event_LO_012hit_first2fittedPosition.pdf'
	shell:
		"""root -q '/home/espedica/macros_fairmu/snakemake/plots_012hits_eff.cpp("{input.path_h}/","first2fittedPosition")'"""

rule check_res_pre_post:
	input:
		path_h='/home/espedica/macros_fairmu/snakemake/plots'
	output:
		out_h='/home/espedica/macros_fairmu/snakemake/plots/results/res_mu_pre_2hit_first2fittedPosition.pdf'
	shell:
		"""root -q '/home/espedica/macros_fairmu/snakemake/plots_pre_post_vrtx.cpp("{input.path_h}/","first2fittedPosition")'"""

rule comparison:
	input:
		path1='/home/espedica/macros_fairmu/snakemake/plots/results',
		path2='/home/espedica/macros_fairmu/snakemake/plots/results'
	output:
		out_c='/home/espedica/macros_fairmu/clean_codes/separate/validation/results/sigma_overlap_first2fittedPosition.pdf'
	shell:
		"""root -q '/home/espedica/macros_fairmu/snakemake/compare_version.cpp("{input.path1}/","{input.path2}/","first2","first2fittedPosition")'"""
