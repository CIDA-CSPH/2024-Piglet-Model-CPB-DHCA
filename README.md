# 2024-Piglet-Model-CPB-DHCA

The data in this Zenodo entry corresponds to the data and results associated with the manuscript titled [Cardiopulmonary bypass with deep hypothermic circulatory arrest results in organ-specific transcriptomic responses in pediatric swine](https://www.sciencedirect.com/science/article/abs/pii/S193152442500012X) and the [Zenodo record](https://zenodo.org/records/14713666). 

To download file from Zenodo which will generate the `results` and `data` directories run
```
sh get_data.sh
```
or directly download the data from Zenodo, unzip the file and place the resulting directories here.

Once the Zenodo data is placed here, this repository will have the following structure

- **data** 
  - **mapped_reads** : generated by nf-core pipeline that is used in the analyses
  - **networks** : data used to generate network clusters
  - **orthologs**: data used to convert swine genes to human genes
  - **Sample info + gender with heart.csv** : sample sheet describing labels of the samples 
- **results**
  - **Make_DEG_Lists** : results generated for the entire gene list (includes results for only looking at all up- or all down-regulated genes, which aren't included in the main manuscript)
  - **Make_Cluster** : results of the network clustering analysis (includes results for only looking at all up- or all down-regulated genes, which aren't included in the main manuscript)
  - **nfcore_count_tables** : count tables generated by the rnaseq pipeline run using nf-core
- **src**
  - **map_reads** : script used to run nf-core pipeline
  - **analyses** : script run to generate the results
  - **renv** : contains information on R packages used
  - **requirements.txt** : contains information on python packages used

Note: the names of the groups in these files are slightly different than in the manuscript  

Controlnodrug = Control  
CPBnodrug = CPB/DHCA  
CPBLind = CPB/DHCA + Lin (Linrodostat analyses ended up not included in the published version of the manuscript)  

### GEO/SRA

The raw fastq files are deposited in GEO/SRA under [GSE294547](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294547).

### Funding

This study was funded by the National Institutes of Health/NHLBI R01HL156936
