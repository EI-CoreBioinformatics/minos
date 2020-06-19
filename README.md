# minos - a gene model consolidation pipeline for genome annotation projects

minos is a Python3/Snakemake - based pipeline that generates and utilises metrics derived from protein, transcript and expression data sets to consolidate gene models obtained from gene annotation workflows. 

For the majority of the computational work, minos utilises Mikado (https://github.com/EI-CoreBioinformatics/Mikado). The pipeline runs Mikado `prepare` on provided gene sets, and generates external metrics such as blastp/blastx alignments, busco assessments and kallisto expression quantification. These metrics are then passed on to Mikado `serialise` and `pick`. In a final set of steps, models are filtered according to user-provided criteria and annotated release gene/transcript sets are generated.

## Installation 

minos requires Python 3.6 at the very least (better Python 3.7+ as it is essential that dictionary insertion order is preserved.) 

### Python dependencies
These dependencies should be installed automatically if not present.

* pyyaml >= 5 (in order to not have the yaml configuration files ordered alphabetically)
* Snakemake >= 5.14.0 (to make use of later Snakemake performance and stability features)
* drmaa (for hpc environments)

### Installation from GitHub
    git clone https://github.com/EI-CoreBioinformatics/minos.git
    cd minos
    python setup.py bdist_wheel;
    pip install dist/*whl

### Installation from conda
    TBD

### Workflow dependencies
minos makes extensive use of 3rd party software, most of which can be installed from conda.
The following tools are required.

* kallisto 0.44.0 (generation of expression based metrics from RNAseq data)
* blast+ (generation of protein similarity metrics from blastp/x analysis of protein evidence data)
* prinseq (fasta processing)
* genometools (fasta validation)
* gffread (extraction of transcript parameters and CDS, cDNA, protein sequences)
* bedtools (transposable element analysis)

The following three dependencies require special treatment (s. below):

* Mikado *
* BUSCO >= 4 (busco protein based metrics, full busco statistics for annotation quality control) **
* CPC2 (http://cpc2.cbi.pku.edu.cn) ***

#### Registering workflow dependencies with minos
Before running minos, the user has to register all dependencies in the `program_calls` section of the run configuration file (etc/minos_config.yaml). This allows users to manage their own installations. For convenience, we provide singularity container definitions for most dependencies as well as a conda recipe. Unfortunately, due to our hpc environment, we cannot provide solutions utilising the native conda/singularity support in current Snakemake version. (However, we would be glad for any community contributions in that regard).

#### Special case: Mikado
- preferred: singularity container due to development cycle and ease of installation
- Due to its complexity and rapid development cycle Mikado should have its own container
- this should all not matter anymore when there is a stable Mikado

#### Special case: BUSCO
- BUSCO requires an older version of blast due to issues with multithreading tblastn
- this is not compatible with having a newer blast+ for the minos protein blast analyses
- our solution: individual BUSCO singularity container
- potential solution: conda environment

#### Special case: CPC2
- CPC2 is not on conda
- CPC2 is not on GitHub
- CPC2 requires Python2
- instructions for a Python3 port of CPC2 are available below (also part of our main dependency container definition)
- alternatively, the user has to ensure a local CPC2 installation (but still should port CPC2 to Python3)
- There is a CPC2 installation/patch script `minos/scripts/install_cpc2.py`, which can be invoked with `install_cpc2 PATH_TO_INSTALLATION` after installation of minos. 

      install_cpc2 PATH_TO_INSTALLATION
      export CPC_HOME=PATH_TO_INSTALLATION/CPC2-beta
      export PATH=$PATH/CPC_HOME/bin

##### Manually porting CPC2 to Python3
    # patch compress.py
    sed -i "s/\bfile(/open(/g" compress.py
    # patch CPC2.py
    sed -i "s/commands/subprocess/g" CPC2.py
    sed -i "s/\bfile(/open(/g" CPC2.py
    sed -i "s/triplet_got.next()/next(triplet_got)/g" CPC2.py

    head -n 355 CPC2.py >> CPC2.py.1
    echo $'\t'$'\t'"os.system('mv ' + outfile + '.txt ' + outfile)" >> CPC2.py.1
    tail -n +356 CPC2.py >> CPC2.py.1
    mv CPC2.py.1 CPC2.py
    # patch seqio.py
    sed -i "s/print \([^ $]\+\)/print(\1)/g" seqio.py
    chmod 775 *.py

## Getting started

minos requires a set of input files and data. To keep minos-based analyses structured it is recommended to collect all input in a project directory, e.g. by soft-linking or copying.

A minos run consists of two steps:
1. Generation of the run configuration (`minos configure`)
2. Running the minos workflow (`minos run`) driven by the run configuration generated in Step 1. 


### minos configure

`minos configure` takes as input a set of configuration files:

1. **list_file** A tab-separated file describing the set of transcript models for the run described in Section transcript_model_info. The files with the transcript models need to be soft-linked into the project directory.
2. **scoring_template** A yaml file containing the scoring settings. This can be copied from the minos repo (minos/etc/scoring_template.yaml) and modified if required. There is no need to add in the external metrics and scoring sections as this will be done automatically by minos configure!
3. **genome_reference** A fasta file containing the genome reference. This should be softlinked instead of copied.
4. **external_metrics_configuration** A Mikado configuration file (e.g. https://github.com/EI-CoreBioinformatics/mikado/blob/master/sample_data/plant_external.yaml)
5. **external_metrics** A tab-separated file describing the metrics data to be used in the minos run. The column descriptions can be found in Section metrics_info.
6. **configuration** A yaml file with parameters to control minos run. This can be obtained from the minos repo (https://github.com/EI-CoreBioinformatics/minos/blob/master/etc/minos_config.yaml).

After these input files have been generated/obtained, `minos configure` can be run. Items 1,2,3 from the above list are positional arguments, items 4,5,6 are optional. [TBC!]

#### Command line arguments

    --outdir, -o

Sets the output directory (default: `minos_run`)

    --blastmode {blastp,blastx}

Controls whether transcript models (blastx) or their translated peptide sequences (blastp, default) are compared against the provided protein evidence.

    --annotation-version ANNOTATION_VERSION, --genus-identifier GENUS_IDENTIFIER

Final output files are named with the prefix `<GENUS_IDENTIFIER>_<ANNOTATION_VERSION>`, e.g. `QUISA32244_EIv1`.

    --use-tpm-for-picking

Controls whether RNA-seq data is used for metrics generation in addition to classification. (default: off) Caution: `--use-tpm-for-picking` activates using expression metrics for the picking stage. *This is supposed to be optional but might be required in the current version, as one of the downstream scripts might not have it coded as optional (which I just realized…or maybe I am mixing this up with the post-picking stage).*

    --force-reconfiguration, -f

Controls whether an existing configuration file in the chosen project directory (set with `-o` or `--outdir`) will be overwritten. By default, this behaviour is turned off.

    --mikado-container

Due to the ongoing Mikado development and frequent updates, the Mikado version (currently required to be served from a Singularity container) used by minos has to be submitted via command line option rather than via the configuration file. This allows for more flexibility in swapping out Mikado versions. As soon as Mikado 2.0 has reached a stable state, we will take all efforts so that this option will no longer be mandatory.

##### BUSCO Options

minos supports core gene assessments with BUSCO4 for metrics generation and quality assessment.

    --busco-level BUSCO_LEVEL

A comma-separated string to select the BUSCO mode. Valid levels are:

    {"proteins", "proteome", "transcripts", "transcriptome", "genome", "none", "off", "p", "t", "g", "a", "all", "prot", "tran", "geno"}.

The default mode is “proteins” (or all in a development version). “none” and “off” disable BUSCO analyses.

    --busco-lineage BUSCO_LINEAGE

This needs to be either a path to the odb10 database used for the BUSCO analysis or a valid database identifier. This option is required if `--busco-level` is not in `{none,off}`. Note, the latter option requires an internet connection and removal of the `--offline` parameter in the run configuration file. 

    --busco-genome-run BUSCO_GENOME_RUN

As an alternative to computing the genome BUSCO analysis, a path can be specified pointing to a directory containing the short_summary.txt and full_table.tsv from a precomputed BUSCO genome run on the reference. Using this option will prompt minos’s BUSCO genome rule to copy the input files to the output folder instead of running BUSCO in genome mode.

BUSCO copy number assessment can be configured by passing/modifying the `--limit` parameters in the `params: busco: <busco_run>` section of the configuration file. Additionally, the maximum copy number that is reported in an individual category in the final BUSCO output table can be set via `misc: busco_max_copy_number`.

**minos configure example run:**

    minos configure --mikado-container </path/to/mikado/container> -o <output-directory> --external external.yaml --external-metrics external_metrics.txt --use-tpm-for-picking --genus-identifier <GENUS_ID> --annotation-version <ANN_VERSION> list.txt </path/to/scoring-template.yaml> </path/to/genome_reference>

EI_internal:
Latest tested mikado-container: `/ei/software/cb/mikado/2.0rc6_d094f99_CBG/x86_64/mikado-2.0rc6_d094f99_CBG.img`

### minos run

`minos run` starts the minos pipeline, controlled by a run configuration file. It is highly recommended to use `minos configure` to generate this configuration file as this will also generate the required Mikado configuration. 

EI_internal:
Note that for convenience, `minos run` has an NBI cluster wrapper called `minos_run_sub`.


**minos run example run:**

    minos_run_sub --mikado-container </path/to/mikado/container> --partition <partition> -o <output-directory> --hpc_config /ei/software/testing/minos/dev/src/minos/minos/etc/hpc_config.json

## Configuration files

### run configuration

The minos workflow is controlled by a run configuration file in yaml format. This file is generated by `minos configure` from a configuration template (https://github.com/EI-CoreBioinformatics/minos/blob/master/etc/minos_run_config.yaml), the command line options, and [...] and saved to the run directory.

The run configuration template contains the following information:

1. `params`

This section allows the user to specify additional parameters for the tools used in the workflow (e.g. score, evalue cutoffs etc). 

2. `program_calls`

This section contains the instructions for running the minos dependencies. If a dependency is installed in the environment, this could just be the tool name. Otherwise, this might include a `singularity exec </path/to/container>` call, a source or module load command, etc.

3. `collapse_metrics_threshold`

This section contains the classification rules for the collapse_metrics rule.

4. `misc`

This section contains parameters for rules not involving 3rd party dependencies and other unsorted things.


### metrics info
A tab-separated file with the following columns:

1. `metric_name_prefix` 
Short name for the metric. 

Currently, the following convention has to be followed:
kallisto metrics need to have a suffix according to their strandedness:

`unstranded: _xx, reverse-forward: _rf, forward-reverse: _fr`

The rf/fr suffixes follow the kallisto command line parameter naming scheme.

2. `metric_class`
Can be one of the following: `expression`, `aln_tran`, `aln_prot`, `seq_prot` (for protein blast db generation), `junction`, `repeat`, `*blastdb_prot` (* not implemented)

3. `multiplier`
The multiplier can be given either as a single value or a comma-separated key:value list (e.g. `aF1:X,eF1:Y,jF1:Z,nF1:W`)

4. `not_fragmentary_min_value`


5. `file_path`
Full path to the metrics file.
- Paired-end reads for kallisto need to be provided as `<path/to/R1>,<path/to/R2>`
- Multiple samples for expression metrics can be given as one line each with the same `metric_name_prefix` (in such a case, the `multiplier` will be read from the line with the first sample)


## Citing minos

TBD

Since minos uses Mikado, please consider to cite:

> Venturini L., Caim S., Kaithakottil G., Mapleson D.L., Swarbreck D. Leveraging multiple transcriptome assembly methods for improved gene structure annotation. GigaScience, Volume 7, Issue 8, 1 August 2018, giy093, [doi:10.1093/gigascience/giy093](https://doi.org/10.1093/gigascience/giy093)


If you also use Portcullis junctions as input metric, please consider to cite:


> Mapleson D.L., Venturini L., Kaithakottil G., Swarbreck D. Efficient and accurate detection of splice junctions from RNAseq with Portcullis. GigaScience, Volume 7, Issue 12, 12 December 2018, giy131, [doi:10.1093/gigascience/giy131](https://doi.org/10.1093/gigascience/giy131)


