#conda activate snakemake
#snakemake rm_yaml --cores 1
#snakemake rm_files --cores 1
#snakemake rm_plots --cores 1
#snakemake --cores 6
#snakemake check_res --cores 1
#snakemake check_eff --cores 1
#snakemake check_fake --cores 1
#snakemake check_res_pre_post0 --cores 1
snakemake comparison0 --cores 1
#snakemake check_res_pre_post1 --cores 1
snakemake comparison1 --cores 1
#snakemake check_res_pre_post2 --cores 1
snakemake comparison2 --cores 1
