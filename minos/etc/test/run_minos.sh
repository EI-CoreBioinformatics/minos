#!/bin/bash -e

# MIKADO=/ei/software/cb/mikado/20190523_42fe752/x86_64/mikado-20190523_42fe752.img
#MIKADO=/ei/software/cb/mikado/20190605_2705cbe/x86_64/mikado-20190605_2705cbe.img
#MIKADO=/ei/software/cb/mikado/20190606_6c8d542/x86_64/mikado-20190606_6c8d542.img
#MIKADO=/ei/software/cb/mikado/2.0rc6_f7f086d_CBG/x86_64/mikado-2.0rc6_f7f086d_CBG.img
MIKADO=/ei/software/cb/mikado/2.0rc6_d094f99_CBG/x86_64/mikado-2.0rc6_d094f99_CBG.img

RUN_DIR=$1
TPM_FOR_PICKING="--use-tpm-for-picking"

etc=$(python -c "import pkg_resources; print(pkg_resources.resource_filename(\"minos\", \"etc\"))")
SCORING_TEMPLATE=${etc}/scoring_template.yaml
HPC_CONFIG=${etc}/hpc_config.json
INPUT_MODELS=list_gtf.txt


minos configure --mikado-container $MIKADO -o ${RUN_DIR} --external external.yaml --external-metrics external_metrics.txt $TPM_FOR_PICKING --genus-identifier ABC --annotation-version EIv2 --busco-lineage /ei/workarea/group-pb/BUSCO_DATABASES/odb10/embryophyta_odb10 $INPUT_MODELS $SCORING_TEMPLATE Melia_azedarach.genome.uc.fasta 
#minos configure --mikado-container $MIKADO -o ${RUN_DIR} --external external.yaml --external-metrics external_metrics.txt $TPM_FOR_PICKING --genus-identifier ABC --annotation-version EIv2 --busco-level off --busco-genome-run $(pwd)/minos_run4/busco/runs/genome/genome/run_embryophyta_odb10 list.txt /ei/software/testing/minos/dev/src/minos/minos/etc/scoring_template.yaml Melia_azedarach.genome.uc.fasta 
#minos configure --mikado-container $MIKADO -o ${RUN_DIR} --external external.yaml --external-metrics external_metrics.txt $TPM_FOR_PICKING --genus-identifier ABC --annotation-version EIv2 --busco-level off --busco-genome-run $(pwd)/minos_run4/busco/runs/genome/genome/run_embryophyta_odb10 list.txt /ei/software/testing/minos/dev/src/minos/minos/etc/scoring_template.yaml Melia_azedarach.genome.uc.fasta 
#/minos_run3/busco_old/runs/genome/genome/run_embryophyta_odb10 
# echo NO RECONFIGURATION!
#minos_run_sub --mikado-container $MIKADO --partition ei-cb -o ${RUN_DIR} --hpc_config $HPC_CONFIG
minos_run_sub --mikado-container $MIKADO --partition ei-cb -o ${RUN_DIR} --hpc_config $(pwd)/hpc_config.json 
