# Comprehensive, residue-level mapping of polyclonal neutralizing antibody responses in BG505 SOSIP trimer vaccinated rabbits
Adam S. Dingens, Payal Pratap, Keara Malone, Thomas Ketas, Sarah Hilton, P.J Klasse, John Moore, Andrew Ward, Jesse D. Bloom (and likely others)

We performed mutational antigenic profiling of BG505 SOSIP trimer vaccinated rabbit serum, provided by John Moore and PJ Klasse. Our first mutational antigenic profiling analysis of escape from PGT151 using the BF520 env libraries was published [here](http://dx.doi.org/10.1016/j.chom.2017.05.003) in June 2017, and this original analysis is located [in this ipython notebook](https://github.com/adingens/BF520_MutationalAntigenicProfiling_PGT151).

Here, we used BG505.T332N mutant Env libraries, first described and characterized in [Haddox, Dingens et al 2018](https://elifesciences.org/articles/34420). Note that this is matched to the BG505 immunogens, enabling profiling of autologous responses. 

## Quick summary
Look at the [notebook results](results/analysis_notebook.md) for an overview of the results.

## Running the analysis
The main analysis is performed by the Jupyter notebook [analysis_notebook.ipynb](analysis_notebook.ipynb).
The Markdown results of the notebook are shown in [results/analysis_notebook.md](results/analysis_notebook.md).

To run [analysis_notebook.ipynb](analysis_notebook.ipynb) and generate the [Markdown results](results/analysis_notebook.md), run the bash script [run_notebook.bash](run_notebook.bash) with:

    ./run_notebook.bash
    
On the Hutch cluster, you probably want to submit this job using [slurm](https://slurm.schedmd.com/), which you can simply do with:

    sbatch -p largenode -c 16 --mem=100000 run_notebook.bash

## Configuring the analysis
The configuration for the analysis is in a separate file, [config.yaml](config.yaml). 
This file defines key variables for the analysis, and should be self-explanatory. 
The [config.yaml](config.yaml) file points to several files in the [./data/](data) subdirectory that specify essential data for the analysis:

  - [data/serum_info.yaml](data/serum_info.yaml):
    YAML file that gives information on all of the serum samples used for selections.
    For each serum there is an entry with the label used in the experiments, then:
      - *name*: a more informative name used when displaying results
      - *description*: description of the serum
      - *group*: group of samples to which serum belong
      - *species*: species from which serum is derived (if relevant)
      - *vaccination*: information of vaccination status (if relevant)

  - [data/sample_list.csv](data/sample_list.csv):
    CSV file giving each sample that was deep sequenced.
    Columns are:
      - *sample*: sample label used in experiments
      - *serum*: serum used for selection
      - *library*: viral library, using simple 1, 2, 3 naming rather than the more confusing library codes used to label experiments
      - *date*: day when sequencing was done
      - *serum_dilution*: dilution of serum used; this includes the 1:4 dilution used during the RDE treatment of the serum. For antibodies, it is the concentration in ug/ml.
      - *percent_infectivity*: percent of viral library retaining infectivity
      - *R1*: glob pattern to R1 FASTQ files on Hutch server; the R2 file names are guessed from the R1 names. If [config.yaml](config.yaml) sets *seq_data_source* to *R1* then there must be a valid R1 file glob for all samples; otherwise this column is ignored.
      - *SRA_accession*: the accession number on the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) (SRA) for the sequencing for this sample. If [config.yaml](config.yaml) sets *seq_data_source* to *SRA_accession* then there must be a valid accession for all samples; otherwise this column is ignored.
  
  
## Results
Results are placed in the [./results/](results) subdirectory.
Many of the results files are **not** tracked in this GitHub repo since they are very large.
However, the following results are tracked:

  - [results/analysis_notebook.md](results/analysis_notebook.md): Markdown rendering of results of the Jupyter notebook [analysis_notebook.ipynb](analysis_notebook.ipynb)
