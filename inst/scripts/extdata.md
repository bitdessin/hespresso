# extdata

Files in the `inst/extdata` directory include model definitions for the HOBIT method,
as well as example datasets for a quick overview of package functionality.

## STAN Files

The Stan files `HOBIT.NB.stan` and `HOBIT.ZINB.stan`,
located in the `inst/extdata` directory,
define the core model structures used by the HOBIT algorithm.
These models are designed to detect homeologs with significant changes
in homeolog expression ratios (HERs) across experimental conditions.

## Expression Seed Matrices

Files prefixed with `seed_matrix.` in the `inst/extdata` directory
contain gene expression matrices used 
to simulate artificial read counts for allopolyploid species.

`seed_matrix.C_flexuosa.tsv.gz` is an expression matrix
derived from RNA-seq data of *Cardamine flexuosa* (Akiyama et al., 2021).
Expression levels were quantified using the HomeoRoq pipeline and normalized via the TMM method.
Genes with zero counts across all replicates, variance >10^9,
or mean >10^5 were filtered out before saving.

`seed_matrix.T_aestivum.tsv.gz` is an expression matrix from RNA-seq data
of allotetraploid wheat (*T. aestivum*) (Yang et al., 2021).
Raw reads were processed using the EAGLE-RC pipeline, normalized with the TMM method,
and filtered to remove genes with zero counts or variance >10^9.

## Sample Datasets

Three types of sample datasets are included in the `inst/extdata` directory
to demonstrate package usage.

`C_flexuosa.tsv.gz` and `C_flexuosa.homeolog.tsv.gz` are expression data and mapping table
from *C. flexuosa* RNA-seq (Akiyama et al., 2021), respectively.
While the original dataset includes various habitats and collection dates with biological replicates,
this sample contains 100 randomly selected homeologs collected on May 16, 2013.
The mapping table is derived from the original study.

`C_insueta.tsv.gz` and `C_insueta.homeolog.tsv.gz` are expression data and mapping table
from *C. insueta* RNA-seq (Sun et al., 2020), respectively.
The full study analyzed leaflets floated on water for 96 hours across nine time points.
This sample dataset includes 100 randomly selected homeologs across the same time points.

`T_aestivum.tsv.gz` and `T_aestivum.homeolog.tsv.gz` are expression data and mapping table
from *T. aestivum* RNA-seq (Yang et al., 2021), respectively.
The original study analyzed two tissues across developmental stages and strains.
The sample dataset here includes only the wild strain at stage W3.5
for comparison of leaf and shoot tissues, with 100 randomly selected homeologs.


# References

- Akiyama R, Sun J, Hatakeyama M, et al. (2021). Fine-scale empirical data on niche divergence and homeolog expression patterns in an allopolyploid and its diploid progenitor species. *New Phytologist*, 229(6):3587–3601. doi:10.1111/nph.17101
- Sun J, Shimizu-Inatsugi R, Hofhuis H, et al. (2020). A Recently Formed Triploid *Cardamine insueta* Inherits Leaf Vivipary and Submergence Tolerance Traits of Parents. *Frontiers in Genetics*, 11:567262. doi:10.3389/fgene.2020.567262
- Yang Y, Zhang X, Wu L, et al. (2021). Transcriptome profiling of developing leaf and shoot apices to reveal the molecular mechanism and co-expression genes responsible for the wheat heading date. *BMC Genomics*, 22(1):468. doi:10.1186/s12864-021-07797-7
