## Antibodies to the SARS-CoV-2 receptor-binding domain that maximize breadth and resistance to viral escape 

For experimental background, see our preprints **[here](https://www.biorxiv.org/content/10.1101/2021.04.06.438709v1)** and **[here](https://www.biorxiv.org/content/10.1101/2021.04.07.438818v1)**.

### What data are shown here?
We are showing mutations to the SARS-CoV-2 RBD that escape antibody binding as measured using mutational antigenic profiling. Raw data is available [here](https://github.com/jbloomlab/SARS-CoV-2-RBD_MAP_Vir_mAbs/blob/main/results/supp_data/vir_antibodies_raw_data.csv).
The drop-down menus can be used to select the escape-mutation maps for each antibody.

When you click on sites, they will be highlighted on the protein structure of the ACE2-bound RBD ([PDB 6M0J](https://www.rcsb.org/structure/6M0J) or antibody-bound RBD structures.

At the site level you can visualize one of two quantities:

 - *total escape* is the sum of the escape from all mutations at a site.

 - *max escape* is the magnitude of the largest-effect escape mutation at each site.

At the mutation level, the height of each letter is proportional to the extent to which that amino-acid mutation escapes antibody binding.
You can color the logo plot letters in four ways:

 - *escape color ACE2 bind* means color letters according to how that mutation affects ACE2 binding as measured in our prior deep mutational scanning ([Starr et al. (2020)](https://doi.org/10.1016/j.cell.2020.08.012), with yellow meaning highly deleterious, and brown meaning neutral or beneficial for ACE2 binding. See Figure 1 in the preprint linked above for color scale.
 
 - *escape color RBD expr* means color letters according to how that mutation affects RBD expression as measured in [Starr et al. (2020)](https://doi.org/10.1016/j.cell.2020.08.012).

 - *escape color gray* means color all letters gray.

 - *escape color func group* means color letters by functional group.
