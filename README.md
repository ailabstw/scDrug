# scDrug: From scRNA-seq to Drug Repositioning

The scDrug constructed a workflow for comprehensive analysis on single-cell RNA sequencing (scRNA-seq) data. It provided a powerful tool with various functions, from fundamental data analysis to drug response prediction, and treatment suggestions.

The scDrug went through three parts on raw scRNA-seq data investigation: **Single-Cell Data Analysis**, **Drug Response Prediction**, and **Treatment Selection**.

- **Single-Cell Data Analysis** performed data preprocessing, clustering, cell type annotation, Gene Set Enrichment Analysis (GSEA), and survival analysis. 

- **Drug Response Prediction** estimated the half maximal inhibitory concentration (IC50) of cell clusters, and reported the cell death percentages as drug kill efficacy.

- **Treatment Selection** listed treatment combinations of given cell clusters.


## Download and Installation

1.  Clone the repository to local directory, e.g., `./scDrug`.

    ```
    git clone https://github.com/ailabstw/scDrug.git ./scDrug
    ```


2.  Build the Docker image tagged `sc-drug`.

    ```
    docker build -t sc-drug ./scDrug
    ```


3.  Run the Docker container named `scDrug` with `/docker/path` mounted to `/server/path` to access files within the Docker container.
    
    ```
    docker run -it --name scDrug -v /server/path:/docker/path --privileged sc-drug
    ```

    
4.  In the Docker container `scDrug`, pull the Docker image `cibersortx/fractions` used in treatment selection.

    ```
    /etc/init.d/docker start
    docker pull cibersortx/fractions
    ```
    
    Note 1: Get `CONTAINER_ID` with command `docker ps -a` and start the container with `docker start -i CONTAINER_ID`.
    Note 2: If docker-in-docker cannot be operated on the user's computer, the user can pull and run the CIBERSORTx container outside the scDrug container as long as a shared folder is mounted on both containers for file sharing.

## Usage

Note: Refer to [example](example) for a detail illustration of the usage for the scDrug.

### Single-Cell Data Analysis

**Single-Cell Data Analysis** took the scRNA-seq data in a 10x-Genomics-formatted mtx directory or a CSV file as input, performed fundamental data analysis, and output a Scanpy Anndata object `scanpyobj.h5ad`, a UMAP `umap_cluster.png` and differentially expressed genes (DEGs) `cluster_DEGs.csv` of the clustering result, and a gene expression profile (GEP) file `GEP.txt`.

Optionally, **Single-Cell Data Analysis** carried out batch correction, cell type annotation and Gene Set Enrichment Analysis (GSEA), and provided additional UMAPs showing batch effects and cell types (`umap_batch.png` and `umap_cell_type.png`), and the GSEA result `GSEA_results.csv`. For cell type annotation, we used [scMatch: a single-cell gene expression profile annotation tool using reference datasets](https://github.com/asrhou/scMatch).

Furthermore, **Single-Cell Data Analysis** could take previously produced Anndata as input and applied sub-clustering on specified clusters.


- Run `python3 single_cell_analysis.py -h` to show the help messages as follow for **Single-Cell Data Analysis**.

```
usage: single_cell_analysis.py [-h] -i INPUT [-f FORMAT] [-o OUTPUT] [-r RESOLUTION] [--impute] [--auto-resolution] [-m METADATA]
                               [-b BATCH] [-c CLUSTERS] [--cname CNAME] [--GEP] [--annotation] [--gsea] [--cpus CPUS] [--survival]
                               [--tcga TCGA] [--id ID] [--prefix PREFIX] [--not_treated]

scRNA-seq data analysis

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to input 10x directory or CSV file
  -f FORMAT, --format FORMAT
                        input format, 10x (default) | csv | h5ad (Anndata object for subclustering with --clusters CLUSTERS)
  -o OUTPUT, --output OUTPUT
                        path to output directory, default='./'
  -r RESOLUTION, --resolution RESOLUTION
                        resolution for clustering, default=0.6
  --impute              do imputation. default: no
  --auto-resolution     automatically determine resolution for clustering
  -m METADATA, --metadata METADATA
                        path to metadata CSV file for batch correction (index as input in first column)
  -b BATCH, --batch BATCH
                        column in metadata (or adata.obs) for batch correction, e.g. 'PatientID'
  -c CLUSTERS, --clusters CLUSTERS
                        perform single cell analysis only on specified clusters, e.g. '1,3,8,9'
  --cname CNAME         which variable should be used when selecting clusters; required when clusters are provided. Default:
                        'louvain'
  --GEP                 generate Gene Expression Profile file.
  --annotation          perform cell type annotation
  --gsea                perform gene set enrichment analysis (GSEA)
  --cpus CPUS           number of CPU used for auto-resolution and annotation, default=1
  --survival            perform survival analysis
  --tcga TCGA           path to TCGA data
  --id ID               Specify TCGA project id in the format "TCGA-xxxx", e.g., "TCGA-LIHC"
  --prefix PREFIX       Any prefix before matrix.mtx, genes.tsv and barcodes.tsv.
  --not_treated         only consider untreated samples from TCGA for survival analysis.
```

- Apply **Single-Cell Data Analysis** with batch correction, clustering resolution 1.0, cell type annotation and GSEA.

```
python3 single_cell_analysis.py --input INPUT --metadata METADATA --batch BATCH --resolution 1.0 --annotation --gsea
```

- **Single-Cell Data Analysis** for sub-clustering on specified clusters at automatically determined resolution run under 4 cpus.

```
python3 single_cell_analysis.py -f h5ad --input scanpyobj.h5ad --clusters CLUSTERS --auto-resolution --cpus 4
```


### Drug Response Prediction

**Drug Response Prediction** examined  `scanpyobj.h5ad` generated in **Single-Cell Data Analysis**, reported clusterwise IC50 and cell death percentages to drugs in the GDSC database via [CaDRReS-Sc](https://github.com/CSB5/CaDRReS-SC) (a recommender system framework for *in silico* drug response prediction), or drug sensitivity AUC in the PRISM database from [DepMap Portal PRISM-19Q4] (https://doi.org/10.1038/s43018-019-0018-6). The output the prediction results are `IC50_prediction.csv` and `drug_kill_prediction.csv` while using parameter `--model GDSC`, and `AUC_prediction.csv` whlie using parameter `--model PRISM`.

- Run `python3 drug_response_prediction.py -h` to show the help messages as follow for **Drug Response Prediction**.

```
usage: drug_response_prediction.py [-h] -i INPUT [-o OUTPUT] [-c CLUSTERS] [-m MODEL] [--n_drugs N_DRUGS]

Drug response prediction

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to input Anndata object (h5ad file)
  -o OUTPUT, --output OUTPUT
                        path to output directory, default='./'
  -c CLUSTERS, --clusters CLUSTERS
                        perform sensitivity prediction on specified clusters, e.g. '1,3,8,9', default='All'
  -m MODEL, --model MODEL
                        the sensitivity screening is from GDSC ic50/PRISM auc, e.g. GDSC, PRISM
  --n_drugs N_DRUGS     the number of drugs to visualize for each cluster
```

- Predict drug response on specified clusters (here for default all clusters) with **Drug Response Prediction**.

```
python3 drug_response_prediction.py --input scanpyobj.h5ad
```


### Treatment Selection

In **Treatment Selection**, we first **imputed cell fractions** of bulk GEPs from the LINCS L1000 database with single-cell GEP `GEP.txt` created in **Single-Cell Data Analysis** via Docker version of [CIBERSORTx Cell Fractions](https://cibersortx.stanford.edu), which enumerated the proportions of distinct cell subpopulations in tissue expression profiles. Then, we **selected treatment combinations** from the LINCS L1000 database with the CIBERSORTx result, and generated plots and a dataframe to show the drug effect.

#### Impute Cell Fractions

**Impute Cell Fractions** took the reference sample file `GEP.txt` as input to run CIBERSORTx Cell Fractions with bulk GEP of user specified or automatically determined cell type, and output CIBERSORTx result files to the output directory, including `CIBERSORTx_Adjusted.txt`. The cell type for bulk GEP involved A375 (malignant melanoma),  A549 (non-small cell lung carcinoma),  HCC515 (non-small cell lung adenocarcinoma),  HEPG2 (hepatocellular carcinoma), HT29 (colorectal adenocarcinoma),  MCF7 (breast adenocarcinoma),  PC3 (prostate adenocarcinoma),  YAPC (Pancreatic carcinoma).

- Run `python3 CIBERSORTx_fractions.py -h` to show the help messages as follow for **Impute Cell Fractions**.

```
usage: CIBERSORTx_fractions.py [-h] -i INPUT [-o OUTPUT] [-l LINCS] [-c CLUSTERS] -u USERNAME -t TOKEN
                               [--celltype CELLTYPE] [--develop]

impute the fractions of previous identified cell subsets under each bulk sample in the LINCS L1000 database.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to input single-cell GEP file
  -o OUTPUT, --output OUTPUT
                        path to output directory, default='./'
  -l LINCS, --lincs LINCS
                        path to LINCS data directory
  -c CLUSTERS, --clusters CLUSTERS
                        perform combined treatment prediction on specified clusters, e.g. '1,3,8,9'
  -u USERNAME, --username USERNAME
                        email address registered on CIBERSORTx website
  -t TOKEN, --token TOKEN
                        token obtained from CIBERSORTx website
  --celltype CELLTYPE   choose a cell line from the options. If no name is provided, we will automatically
                        determine the cell type. Options: A375 (malignant melanoma), A549 (non-small cell lung
                        carcinoma), HCC515 (non-small cell lung adenocarcinoma), HEPG2 (hepatocellular carcinoma),
                        HT29 (colorectal adenocarcinoma), MCF7 (breast adenocarcinoma), PC3 (prostate
                        adenocarcinoma), YAPC (Pancreatic carcinoma)
  --develop             Only for development version.

```

-  **Impute Cell Fractions** via CIBERSORTx Cell Fractions with single-cell GEP `GEP.txt` and LINCS L1000 bulk GEP of automatically determined cell type.

```
python3 CIBERSORTx_fractions.py --input GEP.txt --username USERNAME --token TOKEN
```

Note: To obtain `USERNAME` and `TOKEN`, register and request for access to CIBERSORTx Docker on [CIBERSORTx](https://cibersortx.stanford.edu) website.

#### Select Treatment Combinations

**Select Treatment Combinations** takes the CIBERSORTx result `CIBERSORTx_Results.txt` and the L1000 instance info file as input, selects treatment combinations for the given cell type from the LINCS L1000 database, and output the report of the identified treatment combinations (`treatment_combinations.pdf`).

- Run `python3 treatment_selection.py -h` to show the help messages as follow for **Select Treatment Combinations**.

```
usage: treatment_selection.py [-h] -i INPUT [-o OUTDIR] [-t THRESHOLD] [-c CON_THRESHOLD] --celltype CELLTYPE
                              [--metadata METADATA]

Select treatment combination from the LINCS L1000 database.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        CIBERSORTx output file.
  -o OUTDIR, --outdir OUTDIR
                        path to output directory, default='./'
  -t THRESHOLD, --threshold THRESHOLD
                        Sensitivity threshold. Range: [-1,0), default:-0.9
  -c CON_THRESHOLD, --con_threshold CON_THRESHOLD
                        Consistency threshold. Range: [-1,0), default:-0.75
  --celltype CELLTYPE   Same as the cell type for decomposition. Options: A375 | A549 | HEPG2 | HT29 | MCF7 | PC3 | YAPC
  --metadata METADATA   the L1000 instance info file, e.g., 'GSE70138_Broad_LINCS_inst_info_2017-03-06.txt'
```

- **Select Treatment Combinations** with the L1000 metadata.

```
python3 treatment_selection.py --input CIBERSORTx_Adjusted.txt --celltype CELLTYPE --metadata METADATA
```

