# neoantigen-toolkit

A place for somatic variant calling and neoantigen related scripts


## Includes

### `mhcnuggets-neoantigen_prediction.sh`
For calling HLA typing from RNA-seq data, extracting HLA peptides from a somatic variant VCF and predicting HLA binding using [mhcnuggets](https://github.com/KarchinLab/mhcnuggets)

#### Notes:
-  requires a mhcnuggets env `make_mhcnuggets_env.sh` script to make a conda env. Otherwise, use pip to install and make sure your `pyensembl` is up to date and matches your reference release (eg `pyensembl install --release 108`) 
- [arcasHLA requires](https://github.com/RabadanLab/arcasHLA) pandas, kallisto, bedtools, pigz, samtools, python >= 3.6, biopython, numpy, scipy. The script requires a conda environment for running arcasHLA.
