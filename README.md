# High-resolution mapping of the neutralizing and binding specificities of polyclonal rabbit serum elicited by HIV trimer immunizations
Adam S. Dingens, Payal Pratap, Keara Malone, Thomas Ketas, Sarah K. Hilton, Christopher Cottrell, P.J. Klasse, John Moore, Andrew Ward, Jesse D. Bloom 

We performed mutational antigenic profiling of BG505 SOSIP trimer vaccinated rabbit serum, provided by John Moore and PJ Klasse. In collaboration with Payal Pratap and Adrew Ward's group, we also performed emPEM with these matched sera samples. This analysis focuses only on the mutational antigenic profiling data. 

Here, we used BG505.T332N mutant Env libraries, first described and characterized in [Haddox, Dingens et al 2018](https://elifesciences.org/articles/34420). Note that this is matched to the BG505 immunogens, enabling profiling of autologous responses. 

## Quick summary
Look at the [notebook results](results/analysis_notebook.md) for an overview of the results.

## Running the analysis
The main analysis is performed by the Jupyter notebook [analysis_notebook.ipynb](analysis_notebook.ipynb).
The Markdown results of the notebook are shown in [results/analysis_notebook.md](results/analysis_notebook.md).

To run [analysis_notebook.ipynb](analysis_notebook.ipynb) and generate the [Markdown results](results/analysis_notebook.md), run the bash script [run_notebook.bash](run_notebook.bash) with:

    ./run_notebook.bash
    
On the Hutch cluster, you probably want to submit this job using [slurm](https://slurm.schedmd.com/), which you can simply do with:

    sbatch --partition=campus-new -c 6 run_notebook.bash

## Configuring the analysis
The configuration for the analysis is in a separate file, [config.yaml](config.yaml). 
This file defines key variables for the analysis, and should be self-explanatory. 
The [config.yaml](config.yaml) file points to several files in the [./data/](data) subdirectory that specify essential data for the analysis:

  - [data/serum_info.yaml](data/serum_info.yaml):
    YAML file that gives information on all of the serum samples used for selections.
    For each serum there is an entry with the label used in the experiments, then:
      - *name*: a more informative name used when displaying results
      - *description*: description of the serum
      - *group*: group of samples to which serum belong, based on preliminary point mutant mapping only
      - *vaccination*:  timing of sera "pre" (wk0) or "post" vaccination
      - *immunogen*: the BG505 SOSIP trimer version used as immunogen, and how many times it was administered
      - *dilution*: the serum dilution used

  - [data/sample_list.csv](data/sample_list.csv):
    CSV file giving each sample that was deep sequenced.
    Columns are:
      - *sample*: sample label used in experiments
      - *rabbit_id*: 4 digit rabbit ID
      - *vaccine_status*: timing of sera "pre" (wk0) or "post" vaccination
      - *library*: viral library
      - *dilution*: dilution of serum used
      - *fraction_infectivity*: fraction of viral library retaining infectivity
      - *percent_infectivity*: percent of viral library retaining infectivity
      - *R1*: glob pattern to R1 FASTQ files on Hutch server; the R2 file names are guessed from the R1 names. If [config.yaml](config.yaml) sets *seq_data_source* to *R1* then there must be a valid R1 file glob for all samples; otherwise this column is ignored. UPDATE THIS!
      - *SRA_accession*: the accession number on the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) (SRA) for the sequencing for this sample. If [config.yaml](config.yaml) sets *seq_data_source* to *SRA_accession* then there must be a valid accession for all samples; otherwise this column is ignored.
  
  
  
## Results
Results are placed in the [./results/](results) subdirectory.
Many of the results files are **not** tracked in this GitHub repo since they are very large.
However, the following results are tracked:

  - [results/analysis_notebook.md](results/analysis_notebook.md): Markdown rendering of results of the Jupyter notebook [analysis_notebook.ipynb](analysis_notebook.ipynb)
