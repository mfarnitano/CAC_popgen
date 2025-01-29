# CAC_popgen

Scripts associated with the manuscript "Fluctuating reproductive isolation and stable ancestry structure in a fine-scaled mosaic of hybridizing Mimulus monkeyflowers" by Matthew Farnitano, Keith Karoly, and Andrea Sweigart.

An earlier version of the manuscript is available on bioRxiv at https://doi.org/10.1101/2024.09.18.613726

The current version of the manuscript is currently in review/revision.

1. Creating a reference panel and SNP panels: folder /reference_panels/

- Reference panel samples aligned using fastq_align_wrapper.sh, which calls fastq_align_each.sh
- Reference panel samples genotyped using genotype_wrapper.sh, which calls genotype_each.sh
- Sample genotypes combined into VCF using genotype_combine.sh
- Genotypes filtered to create SNP panels using genotype_filter.sh
- Allele counts within panel subgroups counted using genocounts_groups.py
- Recombination rate added to ancestry-informative sites panel using get_recomb_dist.py for ancestryinfer
- Ancestry-informative sites thinned using thin_positions.py for borice input

2. Aligning low-coverage samples to reference genome: folder /alignment/

- Samples aligned using fastq_align_wrapper.sh, which calls fastq_align_each.sh
- Samples genotyped using genotype_wrapper.sh, which calls genotype_each.sh
- Sample genotypes combined into VCF using genotype_combine.sh
- Genotypes filtered using genotype_filter.sh
- Genome coverage calculated using do_qualimap.sh and summarized with summarize_coverage.sh

3. Using angsd for analyses of population structure: folder /angsd_analyses/

- Genotype likelihoods calculated using do_angsd_GLs.sh
- PCA analysis using do_PCA.sh
- NGSAdmix structure analysis using do_structure.sh

4. Hybrid local ancestry with ancestryinfer and AncestryHMM: folder /ancestry/

- ancestryinfer pipeline called using ancestryinfer_[group].sh using config file ancestryinfer_[group].cfg
- HMM completed to bypass pipeline errors using just_hmm.sh and just_hmm_really.sh
- HMM outputs combined and processed using post_hmm_processing.hmm and combine_ancestry_subset.sh
- Genome-wide hybrid index and heterozygosity per sample calculated with summarize_ancestry_bysample.py
- Site-level ancestry frequencies calculated with summarize_ancestry_bysite.py
- Maternal-offspring pair ancestry calls compared by site using compare_samples_ancestry_submit.sh, which calls compare_samples_ancestry.sh

5. Organellar haplotype networks: folder /organelles/

- Samples aligned to organellar reference sequences using fastq_align_wrapper_organelles.sh, which calls fastq_align_each.sh
- Samples genotyped using genotype_wrapper.sh, which calls genotype_each.sh
- Sample genotypes combined into VCF using genotype_combine_organelles.sh
- Genotypes filtered and combined with previous dataset using genotype_filter_knownsites.sh
- Genotypes processed into gapless alignment using phylofasta.sh

6. Selfing rate estimation with borice: folder /borice/

- Genotype likelihoods converted into borice input format using do_prepborice.sh, which calls likes2borice.py
- Borice input files at [Group]_borice_families.txt, [Group]_borice_families_key.txt, [Group]_borice_subpops.txt
- Borice control files at [Group]_Control.v3.txt
- Borice model run using run_borice.sh

7. Data processing, statistical analysis, and figure generation in R

- Initial data processing and combining data across multiple analyses using Data_prep_and_processing.R
- Statistical models run using models_revision.R
- Phenological RI calculated using RI_calcs.R
- Figures 1,2,4,5, created using Revision_Figure[X].R (supplemental figure code also included in these files)



Note: this repository is meant as a record of the analyses performed, and not as a user-ready package for running these analyses on other datasets. During reorganization for publication, some input and output folder paths and filenames may have been modified from what is listed in the scripts.

Feel free to borrow individual parts of these analyses for your own use, but we warned that modifications will likely be necessary to fit your needs. In particular, changes to filenames and folder organization, cluster-specific parameters and module names, and dataset-specific inputs will be necessary. These scripts are provided with no warranty.
