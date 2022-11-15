## Example Usage for The scDrug

### Preprocessing

The example data is composed of a random 10% subdata from [GSE156625: Onco-fetal reprogramming of endothelial cells drives immunosuppressive macrophages in Hepatocellular Carcinoma (scRNA-seq)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156625).

Before going through the steps in the scDrug, download and uncompress the example data `data.zip ` from [figshare](https://figshare.com/articles/dataset/example_data_zip/20003180) and put the files under the `example` folder. 



### Single-Cell Data Analysis

- First, we execute **Single-Cell Data Analysis** on the 10x-Genomics-formatted mtx directory `data/10x_mtx`, with batch correction of `PatientID` in the metadata `data/metadata.csv`, and clustering at resolution 0.6. Additionally, we assign arguments `--annotation` and `--gsea` to perform cell type annotation and Gene Set Enrichment Analysis (GSEA).
- Required memory: 2.5GB

Note: This step may take a few minutes.

```
mkdir write/clustering

python3 ../script/single_cell_analysis.py --input data/10x_mtx --output write/clustering \
--metadata data/metadata.csv --batch PatientID --resolution 0.6 \
--annotation --gsea --GEP False
```

- Inspecting the preceding output stored in `write/clustering/scanpyobj.h5ad`, we regard the clusters with tumor cell percentages over twice the normal cell percentages, which consist of clusters 1, 5 and 9, as the tumor clusters. Then, we apply **Single-Cell Data Analysis** once again to carry out sub-clustering on the tumor clusters at automatically determined resolution.

Note: To accelerate the process of automatically determined resolution, increase the number of CPU with arguments `--cpus CPUS`.

```
mkdir write/subclustering

python3 ../script/single_cell_analysis.py --input write/clustering/scanpyobj.h5ad --output write/subclustering \
--format h5ad  --clusters '1,5,9' --auto-resolution --cpus 4
```

### Drug Response Prediction

- Based on the sub-clustering result `write/subclustering/scanpyobj.h5ad`, we run **Drug Response Prediction** to predict clusterwise IC50 and cell death percentages to drugs in the GDSC database.
- Required memory: 2GB

```
mkdir write/drug_response_prediction

python3 ../script/drug_response_prediction.py --input write/subclustering/scanpyobj.h5ad --output write/drug_response_prediction
```


### Treatment Selection

In **Treatment Selection**, we first **impute cell fractions** of bulk GEPs from the LINCS L1000 database with single-cell GEP `write/subclustering/GEP.txt`. Then, we **selecte treatment combinations** from the LINCS L1000 database with the CIBERSORTx result, and visualize the result treatment effect.


#### Impute Cell Fractions

- Since it takes several hours to **impute cell fractions**, the result of CIBERSORTx/fractions, `data/CIBERSORTx_Adjusted.txt` and the L1000 instance info file `data/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt`, is provided for the next step.
- Required memory: 18GB

Note: With `USERNAME` and `TOKEN` acquired from [CIBERSORTx](https://cibersortx.stanford.edu), we could also run the following commands to **impute cell fractions** on previously generated `write/subclustering/GEP.txt` with celltype HEPG2 assigned. Notice that this could take several hours.

```
mkdir write/CIBERSORTx_fractions

python3 ../script/CIBERSORTx_fractions.py --input write/subclustering/GEP.txt --output write/CIBERSORTx_fractions \
--username USERNAME --token TOKEN --celltype HEPG2
```

#### Select Treatment Combinations

- With the CIBERSORTx/fractions result `data/CIBERSORTx_Adjusted.txt` and the L1000 instance info file `data/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt` as input, we **select treatment combinations** from the LINCS L1000 database with celltype HEPG2 assigned.
- Required memory: <1GB


```
mkdir write/treatment_selection

python3 ../script/treatment_selection.py --input data/CIBERSORTx_Results.txt --output write/treatment_selection \
--celltype HEPG2 --metadata data/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt
```
