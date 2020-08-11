# High-resolution mapping of the neutralizing and binding specificities of polyclonal rabbit serum elicited by HIV trimer immunizations
Adam S. Dingens, Payal Pratap, Keara Malone, Thomas Ketas, Sarah K. Hilton, Christopher Cottrell, P.J Klasse, John Moore, Andrew Ward, Jesse D. Bloom 

Below are links to interactive visualizations of the HIV mutational antigenic profiling of sera from rabbits vaccinated with BG505 trimers, enabled by [dms-view](https://dms-view.github.io/docs/). This data is posted in [this pre-print - broken link](link).

[Figure 2 - Plotting both the glycan hole and C3/V5 epitopes]

[Figure 4 - Plotting the C3/V5 epitope]

[Figure 5 - Plotting the glycan hole epitope]

For dms-view visualizations, logo plots are colored as in Figure 2 of the paper. The glycan hole epitope is green, the C3/V5 epitope is blue, and all other sites are grey. These residue-level epitope annotations are arbitrarily defined based on all data, we encourage readers to further explore additional sites. Mutations tested in preliminary point mutant mapping (Figure 1A) are shown in black. Selection is mapped onto the [5fyl BG505 trimer or monomer structure](https://www.rcsb.org/structure/5FYL).

[This page](https://jbloomlab.github.io/dms_tools2/diffsel.html) documents the differential selection statsitics we use to analyze these data.

The site-metrics (dot plot) include:

- **positive diffsel**: The sum of all positive differential selection values at a site (site level). This gives a sense to the total amount of escape/selective pressure at each site.
- **negative diffsel**: The sum of all negative differential selection values at a site (site level). This gives a sense for mutations that are depleted, rather than enriched, during serum selection relative to a non-selected control library. It is intriguing that many of these potential serum sensitizing mutations cluster and are consistent across sera.
- **max diffsel**: The value of the largest effect mutation (largest mutation differential selection) at each site (site level).

The mutation-metrics (logoplot) include

- **diffsel**
- **pos diffsel**: only the mutations with positive differential selection values. 

Additionally,

- The frequency at which each amino acid is found in nature (**Natural Frequencies**), accessed from [LANL's filtered web alignment](https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html]) (2017 version).
- The BG505 amino-acid preferences (**DMS preferences**), determined using the same BG505.T332N mutant virus libraries in [Haddox, Dingens et al 2018](https://elifesciences.org/articles/34420). Here, the height of each amino acid is proportional to how well that virus replicates in cell culture. This statistic can crudely be used to examine what mutations are viable and in our mutant virus libraries before serum selection.