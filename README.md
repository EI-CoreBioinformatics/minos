# minos - a gene model consolidation pipeline for genome annotation projects

minos is a Python3/Snakemake - based pipeline that generates and utilises metrics derived from protein, transcript and expression data sets to consolidate gene models obtained from gene annotation workflows. 

For the majority of the computational work, minos utilises [Mikado](https://github.com/EI-CoreBioinformatics/Mikado). The pipeline runs Mikado `prepare` on provided gene sets, and generates external metrics such as blastp/blastx alignments, busco assessments and kallisto expression quantification. These metrics are then passed on to Mikado `serialise` and `pick`. In a final set of steps, models are filtered according to user-provided criteria and annotated release gene/transcript sets are generated.

## Workflow
![Alt text](/minos/doc/Minos.png)
Figure 1. The overview of MINOS pipeline

## Installation 

minos requires Python 3.6 at the very least (better Python 3.7+ as it is essential that dictionary insertion order is preserved.) 

### Python dependencies
These dependencies should be installed automatically if not present.

* pyyaml >= 5 (in order to not have the yaml configuration files ordered alphabetically)
* Snakemake >= 5.14.0 (to make use of later Snakemake performance and stability features)
* drmaa (for hpc environments)

### Installation from GitHub
```console
git clone https://github.com/EI-CoreBioinformatics/minos.git
cd minos
python setup.py bdist_wheel;
pip install dist/*whl
```

### Installation from conda
    TBD

### Workflow dependencies
minos makes extensive use of 3rd party software, most of which can be installed from conda.
The following tools are required and these tools can be installed using singularity container definitions provided ([Singularity.tools.def](minos/etc/Singularity.tools.def).

* kallisto 0.44.0 (generation of expression based metrics from RNAseq data)
* blast+ (generation of protein similarity metrics from blastp/x analysis of protein evidence data)
* seqkit (fasta processing)
* genometools (GFF3 validation)
* gffread (extraction of transcript parameters and CDS, cDNA, protein sequences)
* bedtools (transposable element analysis)
* CPC2 (http://cpc2.cbi.pku.edu.cn) (assess the protein-coding potential of transcripts)

The following three dependencies require special treatment (s. below):

* Mikado *
* BUSCO >= 4 (busco protein based metrics, full busco statistics for annotation quality control) **
* CPC2 (http://cpc2.cbi.pku.edu.cn) ***

#### Registering workflow dependencies with minos
Before running minos, the user has to register all dependencies in the `program_calls` section of the run configuration file (etc/minos_config.yaml). This allows users to manage their own installations. For convenience, we provide singularity container definitions ([Singularity.tools.def](minos/etc/Singularity.tools.def)) for most dependencies, except BUSCO, as well as a conda recipe. Unfortunately, due to our hpc environment, we cannot provide solutions utilising the native conda/singularity support in current Snakemake version. (However, we would be glad for any community contributions in that regard).

#### Special case: Mikado
- preferred: singularity container due to development cycle and ease of installation
- Due to its complexity and rapid development cycle Mikado should have its own container
- this should all not matter anymore when there is a stable Mikado

#### Special case: BUSCO
- BUSCO requires an older version of blast due to issues with multithreading tblastn
- this is not compatible with having a newer blast+ for the minos protein blast analyses
- our solution: individual BUSCO singularity container [Singularity.busco.def](minos/etc/Singularity.busco.def)
- potential solution: conda environment

#### Special case: CPC2
- CPC2 is forked to our repository and ported to Python3 and can be installed from here - https://github.com/EI-CoreBioinformatics/CPC2
- Installation instructions are also added to the singularity container definition file ([Singularity.tools.def](minos/etc/Singularity.tools.def)

## Getting started

minos requires a set of input files and data. To keep minos-based analyses structured it is recommended to collect all input in a project directory, e.g. by soft-linking or copying.

A minos run consists of two steps:
1. Generation of the run configuration (`minos configure`)
2. Running the minos workflow (`minos run`) driven by the run configuration generated in Step 1. 


### minos configure

`minos configure` takes as input a set of configuration files:

1. **list_file** A tab-separated file describing the set of transcript models. An example file is provided here [list.txt](data/list.txt). A detailed description for each column header can be found here - [Mikado list file format](https://mikado.readthedocs.io/en/stable/Usage/Configure/#input-annotation-files)
2. **scoring_template** A yaml file containing the scoring settings. This can be copied from the minos repo [minos/etc/scoring_template.yaml](minos/etc/scoring_template.yaml) and modified if required. There is no need to add in the external metrics and scoring sections as this will be done automatically by minos configure!
3. **genome_reference** A fasta file containing the genome reference. This should be softlinked instead of copied.
4. **external_metrics_configuration** A Mikado configuration file (e.g. [sample_data/plant_external.yaml](https://github.com/EI-CoreBioinformatics/mikado/blob/master/sample_data/plant_external.yaml))
5. **external_metrics** A tab-separated file describing the metrics data to be used in the minos run. The column descriptions can be found in [Section metrics_info](#metrics-info).
6. **configuration** A yaml file with parameters to control minos run. This can be obtained from the minos repo ([minos/etc/minos_config.yaml](minos/etc/minos_config.yaml)).

After these input files have been generated/obtained, `minos configure` can be run. Items 1,2,3 from the above list are positional arguments, items 4,5,6 are optional. [TBC!]

#### Command line arguments

    --outdir, -o

Sets the output directory (default: `minos_run`)

    --blastmode {blastp,blastx}

Controls whether transcript models (blastx) or their translated peptide sequences (blastp, default) are compared against the provided protein evidence.

    --annotation-version ANNOTATION_VERSION, --genus-identifier GENUS_IDENTIFIER

Final output files are named with the prefix `<GENUS_IDENTIFIER>_<ANNOTATION_VERSION>`, e.g. `QUISA32244_EIv1`.

    --use-tpm-for-picking

Controls whether RNA-seq data is used for metrics generation in addition to classification. (default: off) Caution: `--use-tpm-for-picking` activates using expression metrics for the picking stage.

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

```console
minos configure --mikado-container </path/to/mikado/container> \
    -o <output-directory> --external external.yaml \
    --external-metrics external_metrics.txt --use-tpm-for-picking \
    --genus-identifier <GENUS_ID> --annotation-version <ANN_VERSION> \
    list.txt </path/to/scoring-template.yaml> </path/to/genome_reference>
```

### minos run

`minos run` starts the minos pipeline, controlled by a run configuration file. It is highly recommended to use `minos configure` to generate this configuration file as this will also generate the required Mikado configuration. 

EI_internal:
Note that for convenience, `minos run` has an NBI cluster wrapper called `minos_run_sub`.


**minos run example run:**

From local machine:

```console
minos run --mikado-container </path/to/mikado/container> \
    --no_drmaa --scheduler NONE -o <output-directory>
```

On HPC (default `SLURM`, an example HPC config JSON we use is here [hpc_config.json](minos/etc/hpc_config.json)):

```console
minos_run_sub --mikado-container </path/to/mikado/container> \
    --partition <partition> --hpc_config /path/to/hpc_config.json \
    -o <output-directory>
```

## Configuration files

### run configuration

The minos workflow is controlled by a run configuration file in yaml format. This file is generated by `minos configure` from a configuration template ([minos_config.yaml](minos/etc/minos_config.yaml)), and is saved to the output directory.

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


