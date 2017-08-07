[Acanthaster](http://marinegenomics.oist.jp/cots/viewer/download?project_id=46) is seafish with eyes located on its terminal tube feet. In efforts to study the development of vision and light sensing in a opsin dependent manner, we sequenced several tissue specific samples (w/ 3 replicates each) for differential expression (DE) analysis using Illimina sequencers.

The raw data for the reads has been uploaded to NCBI under the accession number TBD. Scripts used for this analysis are located in the scripts folder.

- reads were quasi-mapped using [salmon](http://salmon.readthedocs.io/en/latest/) against the gbr Acanthaster transcriptome.
- Differentially expressed transcripts were identified using [DESeq2](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).
- opsin sequences were idenified by reciprical best hit blast using sequences collected from various sources.
- opsin sequences were then aligned using [MAFFT](http://mafft.cbrc.jp/alignment/software/) v7.215 and phylogenomic trees were then generated using [FastTree](http://www.microbesonline.org/fasttree/) (v1.4.3) to better classify the retrived opsins.

