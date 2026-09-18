# Too Many Isoforms
This code is associated with the bioinformatic methods and analysis presented in the manuscript:

**Gene Body Architecture Predicts Isoform Diversity Across Eukaryotic Genomes**

Quinton Brohart<sup>1,2</sup>, Arvind Mer<sup>1,2,3,4,#</sup>

1. Department of Biochemistry, Microbiology & Immunology, University of Ottawa, Ottawa, Canada
2. Ottawa Institute of Systems Biology, Ottawa, Ontario, Canada
3. School of Electrical Engineering & Computer Science, University of Ottawa, Ottawa, Canada
4. Ottawa Hospital Research Institute, Ottawa, Ontario, Canada

\# Corresponding author

**ABSTRACT**

> *Alternative splicing generates transcript and isoform diversity, but the gene-level features that distinguish isoform-rich genes from simpler genes remain incompletely defined. We examined whether structural, regulatory, sequence-derived, and functional features jointly explain isoform complexity in the human genome and whether these relationships are conserved across mouse, Drosophila, and zebrafish genomes as well. Using curated transcript annotations and gene-level biomolecular features, we found that isoform diversity is pervasive across all four species but varies by genome and biotype. Across species, isoform-rich genes were consistently associated with longer introns, shorter exons, and increased exon-intron organization. UTR associations were more variable, with 3′ UTR length showing stronger positive associations in humans and Drosophila. In humans, isoform count was negatively associated with GC content and CpG island density and positively associated with protein-protein interaction network size. Machine learning based modeling predicted human isoform counts with strong performance and retained predictive signal in mouse. PCA and UMAP further showed that isoform-associated gene features form a continuous, structured manifold rather than discrete gene classes. These findings support a model in which isoform diversity is shaped by conserved gene-structural features and modulated by lineage- and biotype-specific regulatory architecture. Integrating predictive modeling, feature attribution, and dimensionality reduction provides complementary insight into the genomic organization underlying transcript diversity.*

## Steps to Reproduce the Results:
1. Download the manuscript-associated data files from figshare at https://doi.org/10.6084/m9.figshare.32748741.
2. Place the downloaded files in the `./Data/` folder of the repository, corresponding to the species of interest.
3. Run the scripts in the `./Scripts/` folder. 
4. Run the script pertaining to the figure or analysis of interest.

## Example
To reproduce Figure 3 in the manuscript, run the NCAPD2_Transcript_Breakdown.R script. Run the script from the repository root directory (the ./tooManyIso folder) with the appropriate annotated data files.

```
Rscript ./Scripts/NCAPD2_Transcript_Breakdown.R
```

Results:
<img width="1760" height="1760" alt="NCAPD2_Robust_Manuscript" src="https://github.com/user-attachments/assets/a3034b47-c3bb-4776-aaaa-d2c3f37cd5f4" />
