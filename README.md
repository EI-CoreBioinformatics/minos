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

##### Porting CPC2 to Python3
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


## Citing minos

TBD

Since minos uses Mikado, please consider to cite:

> Venturini L., Caim S., Kaithakottil G., Mapleson D.L., Swarbreck D. Leveraging multiple transcriptome assembly methods for improved gene structure annotation. GigaScience, Volume 7, Issue 8, 1 August 2018, giy093, [doi:10.1093/gigascience/giy093](https://doi.org/10.1093/gigascience/giy093)


If you also use Portcullis junctions as input metric, please consider to cite:


> Mapleson D.L., Venturini L., Kaithakottil G., Swarbreck D. Efficient and accurate detection of splice junctions from RNAseq with Portcullis. GigaScience, Volume 7, Issue 12, 12 December 2018, giy131, [doi:10.1093/gigascience/giy131](https://doi.org/10.1093/gigascience/giy131)


