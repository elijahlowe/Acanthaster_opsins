[Acanthaster](http://marinegenomics.oist.jp/cots/viewer/download?project_id=46) is starfish with eyes located on its terminal tube feet. In efforts to study the development of vision and light sensing in a opsin dependent manner, we sequenced tissues for A. planci eyes, radial nerves, tube feet and a mix of other organ tissue (w/ 3 replicates each) for differential expression (DE) analysis using Illimina sequencers.

The preprint for this analysis and results can be found on [biorxiv](http://www.biorxiv.org/content/early/2017/09/13/173187)

The raw data for the reads has been uploaded to NCBI under the accession number TBD. Scripts used for this analysis are located in the scripts folder.

## Phylogenomic analysis
- opsin sequences were idenified by reciprical best hit blast using sequences collected from various sources. The
- opsin sequences were then aligned using [MAFFT](http://mafft.cbrc.jp/alignment/software/) v7.215 and phylogenomic trees were then generated using [FastTree](http://www.microbesonline.org/fasttree/) (v1.4.3) to better classify the retrived opsins.
- Additional trees were then generated using [iqtree](http://www.iqtree.org/) to confirm the placement of the opsin sequences. Tree support was 10,000 [Ultrafast Bootstrap Approximation](https://academic.oup.com/mbe/article-lookup/doi/10.1093/molbev/mst024) with a LG+F+L6 amino acid substution model identifed by [ModelFinder](https://www.nature.com/nmeth/journal/v14/n6/full/nmeth.4285.html).

## Read quanification
- reads were quasi-mapped using [salmon](http://salmon.readthedocs.io/en/latest/) against the [gbr Acanthaster transcriptome](http://marinegenomics.oist.jp/cots/viewer/download?project_id=46).
- Differentially expressed transcripts were identified using [DESeq2](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).
- Additional visualization was done using the [karyoploteR](https://bernatgel.github.io/karyoploter_tutorial/) R package.
