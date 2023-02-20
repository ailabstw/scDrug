import argparse, os, sys
from os import listdir
from os.path import isfile, join

import pickle
import math
import collections
import time

import numpy as np
import pandas as pd
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

from sklearn.metrics import silhouette_score
import multiprocess as mp
from functools import partial

import kaplanmeier as km
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

mpl.rcParams['figure.autolayout'] = True
mpl.rcParams['font.size'] = 10
mpl.rcParams['figure.dpi'] = 300
sc.set_figure_params(dpi=300)

def runGSEAPY(adata, group_by='louvain', cutoff=0.05, logfc_threshold=2, outdir='.'):
    import gseapy as gp

    print('Running GSEAPY...')
    start = time.time()
    
    with open('/scDrug/data/GO_Biological_Process_2021.pkl', 'rb') as handle:
        gene_sets = pickle.load(handle)

    df_list = []
    cluster_list = []
    celltypes = sorted(adata.obs[group_by].unique())

    for celltype in celltypes:
        degs = sc.get.rank_genes_groups_df(adata, group=celltype, key='rank_genes_groups', log2fc_min=logfc_threshold, 
                                    pval_cutoff=cutoff)['names'].squeeze()
        if isinstance(degs, str):
            degs = [degs.strip()]
        else:
            degs = degs.str.strip().tolist()
        
        if not degs:
            continue

        enr = gp.enrichr(gene_list=degs,
                gene_sets=gene_sets,
                no_plot=True
                )
        if (enr is not None) and hasattr(enr, 'res2d') and (enr.res2d.shape[0] > 0):
            df_list.append(enr.res2d)
            cluster_list.append(celltype)

    columns = ['Cluster', 'Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Genes']

    df = pd.DataFrame(columns = columns)
    for cluster_ind, df_ in zip(cluster_list, df_list):
        df_ = df_[df_['Adjusted P-value'] <= cutoff]
        df_ = df_.assign(Cluster = cluster_ind)
        if (df_.shape[0] > 0):
            df = pd.concat([df, df_[columns]], sort=False)
            df_tmp = df_.loc[:, ['Term', 'Adjusted P-value']][:min(10, df_.shape[0])]
            df_tmp['Term'] = [x.split('(',1)[0] for x in df_tmp['Term']]
            df_tmp['-log_adj_p'] = - np.log10(df_tmp['Adjusted P-value'])
            df_tmp = df_tmp.sort_values(by='-log_adj_p', ascending=True)
            ax = df_tmp.plot.barh(y='-log_adj_p', x='Term', legend=False, grid=False, figsize=(12,4))
            ax.set_title('Cluster {}'.format(cluster_ind))
            ax.set_ylabel('')
            ax.set_xlabel('-log(Adjusted P-value)')
            pp.savefig(ax.figure, bbox_inches='tight')
            plt.close()
        else:
            print('No pathway with an adjusted P-value less than the cutoff (={}) for cluster {}'.format(cutoff, cluster_ind))
    df.to_csv(os.path.join(args.output, 'GSEA_results.csv'))

    end = time.time()
    print('time: {}', end-start)


def getBulkProfile(bulkpath, gencode_table):
    bk_gep = pd.read_csv(bulkpath, sep=',', compression='gzip')
    # covert gene ID to gene name
    bk_gep_name = bk_gep.merge(gencode_table, how='left', on='id')
    bk_gep_name.drop_duplicates(subset='name', inplace=True)
    bk_gep_name = bk_gep_name.set_index('name')
    bk_gep_name = bk_gep_name.drop(columns=['id'])
    return bk_gep_name

def getSpecCellDict(bk_gep_name, dict_deg):
    score_table = pd.DataFrame(index = bk_gep_name.columns)
    for gene in bk_gep_name.index:
        median = bk_gep_name.loc[gene, :].median()
        score_table[gene] = list(bk_gep_name.loc[gene, :] >= median)
    score_table = score_table.T

    spec_score_table = pd.DataFrame(index = dict_deg.keys())
    for sample in score_table.columns:
        values = []
        for _, genes in dict_deg.items():
            values.append(score_table.loc[genes, sample].sum())
        spec_score_table[sample] = values
    spec_score_table = spec_score_table.T
    
    dict_low  = {}
    dict_high = {}
    for col in spec_score_table.columns:
        dict_low[col] = spec_score_table[col].quantile(q=0.25)
        dict_high[col] = spec_score_table[col].quantile(q=0.75)
        
    dict_celltype = collections.defaultdict(lambda : collections.defaultdict(dict))
    for c in spec_score_table.columns.tolist():
        dict_celltype[c]['high'] = [x for x in spec_score_table.index if spec_score_table.loc[x,c] >= dict_high[c]]
        dict_celltype[c]['low'] = [x for x in spec_score_table.index if spec_score_table.loc[x,c] <= dict_low[c]]
    return dict_celltype

def drawSurvivalPlot(dict_celltype, clinical_df, project_id):
    n_types = len(dict_celltype.keys())
    n_cols = 3
    n_rows = math.ceil(n_types/n_cols)
    dict_group = {'high':1, 'low':0}
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(2.5*n_cols, 2.5*n_rows))
    i_col = 0
    i_row = 0
    alpha_val = 0.05
    cph_selected_cols = ['project', 'cell type', 'coef', 'exp(coef)', 'se(coef)', 'coef lower 95%', 'coef upper 95%','exp(coef) lower 95%', 'exp(coef) upper 95%', 'z', 'p', '-log2(p)']
    df_hazard = pd.DataFrame(columns = cph_selected_cols)

    for c in sorted(dict_celltype.keys(), key=int):
        ax = axs[i_row][i_col]
        dict_km = {'time':[], 'died':[], 'group':[], 'abundance':[]}
        for g in ['high', 'low']:
            death_time = [int(x) for x in clinical_df.loc[dict_celltype[c][g], '_event_time_']]
            dict_km['time'].extend(death_time)
            dict_km['died'].extend([not x.startswith('\'') for x in clinical_df.loc[dict_celltype[c][g], 'days_to_death']])
            dict_km['group'].extend([g]*len(death_time))
            dict_km['abundance'].extend([dict_group[g]]*len(death_time))

        df = pd.DataFrame(dict_km)
        T = df['time']
        E = df['died']
        dem = (df['group'] == 'high')

        # Compute survival
        kmf = KaplanMeierFitter(alpha=alpha_val)
        kmf.fit(T[dem], event_observed=E[dem], label='High')
        kmf.plot_survival_function(ax=ax, ci_show=False)
        kmf.fit(T[~dem], event_observed=E[~dem], label='Low')
        kmf.plot_survival_function(ax=ax, ci_show=False)
        out_lr = logrank_test(T[dem], T[~dem], E[dem], E[~dem], alpha=1-alpha_val)
        p_value = out_lr.p_value

        # Hazard Rate
        df.drop(['group'], axis=1, inplace=True)
        cph = CoxPHFitter()
        cph.fit(df, duration_col='time', event_col='died')
        hr = cph.hazard_ratios_.values[0]
        dict_tmp = cph.summary.to_dict(orient='records')[0]
        dict_tmp['project'] = project_id
        dict_tmp['cell type'] = c
        df_hazard = df_hazard.append(dict_tmp, ignore_index=True)

        ax.set_title('C{} (Log-rank P={:.2f}, HR={:.2f})'.format(c, p_value, hr), fontsize=8, y=1.08)

        # update i_col, i_row
        i_col += 1
        if i_col >= n_cols: 
            i_col = 0
            i_row += 1

    # remove empty axes
    while (i_col < n_cols and i_row < n_rows):
        fig.delaxes(axs[i_row][i_col])
        i_col += 1

    fig.suptitle('Survival Analysis: '+project_id.rsplit('-',1)[1], fontsize=12, y=1)
    fig.tight_layout()
    plt.tight_layout()
    pp.savefig(fig, bbox_inches='tight')
    plt.close('all')
    
    return df_hazard

def readFile():
    # Output: adata with raw count and metadata
    print('Reading files...')
    start = time.time()
    # check metadata
    if not args.metadata is None:
        if not os.path.exists(args.metadata):
            sys.exit("The metadata file does not exist.")
        else:
            metadata_df = pd.read_csv(args.metadata, index_col=0)
    # read main file
    if args.format == 'h5ad':
        adata = sc.read(args.input)
        if not adata.raw is None:
            adata = adata.raw.to_adata()
        if adata.X.max() < 50:
            print("Warning: The input h5ad doesn't contain raw count information.")
            adata.X = np.expm1(adata.X)
        # if not adata.obs is None and not args.metadata is None:
        #     newinfo = metadata_df.columns.difference(adata.obs.columns)
        #     if len(newinfo) > 0:
        #         adata.obs = adata.obs.merge(metadata_df.loc[:, newinfo], left_index=True, right_index=True)
    elif args.format == 'csv':
        adata = sc.read_csv(args.input)
        if not args.metadata is None:
            adata.obs = metadata_df
    else:
        adata = sc.read_10x_mtx(args.input, var_names='gene_symbols', cache=True, prefix=args.prefix)
        if not args.metadata is None:
            adata.obs = metadata_df
    adata.var_names_make_unique()
    # filter clusters 
    if args.clusters:
        clusters = [x.strip() for x in args.clusters.split(',')]
        if args.cname in adata.obs:
            adata = adata[adata.obs[args.cname].isin(clusters)]
        else:
            sys.exit(f"{args.cname} cannot be found in data.")
    end = time.time()
    print('time: {}', end-start)
    return adata

def writeGEP(adata_GEP):
    print('Exporting GEP...')
    sc.pp.normalize_total(adata_GEP, target_sum=1e6)
    mat = adata_GEP.X.transpose()
    if type(mat) is not np.ndarray:
        mat = mat.toarray()
    GEP_df = pd.DataFrame(mat, index=adata_GEP.var.index)
    GEP_df.columns = adata.obs['louvain'].tolist()
    # GEP_df = GEP_df.loc[adata.var.index[adata.var.highly_variable==True]]
    GEP_df.dropna(axis=1, inplace=True)
    GEP_df.to_csv(os.path.join(args.output, 'GEP.txt'), sep='\t')

def preprocessing(adata):
    print('Preprocessing...')
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    if not (adata.obs.pct_counts_mt == 0).all():
        adata = adata[adata.obs.pct_counts_mt < 30, :]
    sc.pp.normalize_total(adata, target_sum=1e4)
    adata.raw = adata.copy()
    sc.pp.log1p(adata)
    if args.impute:
        print('Doing imputation...')
        start = time.time()
        bk_adata_raw = adata.raw.to_adata()
        sc.external.pp.magic(adata)
        adata.raw = bk_adata_raw
        end = time.time()
        print('time: {}', end-start)
    sc.pp.highly_variable_genes(adata)
    if 'features' in adata.raw.var.columns:
        adata.raw.var.index = adata.raw.var['features']
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver='arpack')
    return adata

def batchCorrect(adata):
    if args.batch in adata.obs.columns:
        print('Batch correction...')
        start = time.time()
        adata.obs[args.batch] = adata.obs[args.batch].astype(str)
        sc.external.pp.harmony_integrate(adata, args.batch, adjusted_basis='X_pca')
        end = time.time()
        print('time: {}', end-start)
    else:
        print('Skip batch correction...')
    return adata

def autoResolution(adata):
    print("Automatically determine clustering resolution...")
    start = time.time()
    def subsample_clustering(adata, sample_n, subsample_n, resolution, subsample):
        subadata = adata[subsample]
        sc.tl.louvain(subadata, resolution=resolution)
        cluster = subadata.obs['louvain'].tolist()
        
        subsampling_n = np.zeros((sample_n, sample_n), dtype=bool)
        coclustering_n = np.zeros((sample_n, sample_n), dtype=bool)
        
        for i in range(subsample_n):
            for j in range(subsample_n):
                x = subsample[i]
                y = subsample[j]
                subsampling_n[x][y] = True
                if cluster[i] == cluster[j]:
                    coclustering_n[x][y] = True
        return (subsampling_n, coclustering_n)
    rep_n = 5
    subset = 0.8
    sample_n = len(adata.obs)
    subsample_n = int(sample_n * subset)
    resolutions = np.linspace(0.4, 1.4, 6)
    silhouette_avg = {}
    np.random.seed(1)
    best_resolution = 0
    highest_sil = 0
    for r in resolutions:
        r = np.round(r, 1)
        print("Clustering test: resolution = ", r)
        sub_start = time.time()
        subsamples = [np.random.choice(sample_n, subsample_n, replace=False) for t in range(rep_n)]
        p = mp.Pool(args.cpus)
        func = partial(subsample_clustering, adata, sample_n, subsample_n, r)
        resultList = p.map(func, subsamples)
        p.close()
        p.join()
        
        subsampling_n = sum([result[0] for result in resultList])
        coclustering_n = sum([result[1] for result in resultList])
        
        subsampling_n[np.where(subsampling_n == 0)] = 1e6
        distance = 1.0 - coclustering_n / subsampling_n
        np.fill_diagonal(distance, 0.0)
        
        sc.tl.louvain(adata, resolution=r, key_added = 'louvain_r' + str(r))
        silhouette_avg[str(r)] = silhouette_score(distance, adata.obs['louvain_r' + str(r)], metric="precomputed")
        if silhouette_avg[str(r)] > highest_sil:
            highest_sil = silhouette_avg[str(r)]
            best_resolution = r
        print("robustness score = ", silhouette_avg[str(r)])
        sub_end = time.time()
        print('time: {}', sub_end - sub_start)
        print()
    adata.obs['louvain'] = adata.obs['louvain_r' + str(best_resolution)]
    print("resolution with highest score: ", best_resolution)
    res = best_resolution
    # write silhouette record to uns and remove the clustering results except for the one with the best resolution
    adata.uns['sihouette score'] = silhouette_avg
    # draw lineplot
    df_sil = pd.DataFrame(silhouette_avg.values(), columns=['silhouette score'], index=[float(x) for x in silhouette_avg.keys()])
    df_sil.plot.line(style='.-', color='green', title='Auto Resolution', xticks=resolutions, xlabel='resolution', ylabel='silhouette score', legend=False)
    pp.savefig()
    plt.close()
    end = time.time()
    print('time: {}', end-start)
    return adata, res

def clustering(adata):
    sc.pp.neighbors(adata, n_pcs=20)
    sc.tl.umap(adata)

    if args.auto_resolution:
        adata, res = autoResolution(adata)
    else:
        print("Clustering with resolution = ", args.resolution)
        sc.tl.louvain(adata, resolution=args.resolution)
        res = args.resolution

    print('Exporting UMAP...')
    sc.settings.autosave = True
    sc.settings.figdir = args.output
    fig, ax = plt.subplots(dpi=300, figsize=(12,12))
    sc.pl.umap(adata, color='louvain', ax=ax, legend_loc='on data',
                frameon=False, legend_fontsize='small', legend_fontoutline=2, legend_fontweight='normal',
                use_raw=False, show=False, title='louvain, resolution='+str(res))
    pp.savefig(fig, bbox_inches='tight')
    plt.close()

    if not args.batch is None:
        fig, ax = plt.subplots(dpi=300, figsize=(12,12))
        sc.pl.umap(adata, color=[args.batch], ax=ax, legend_loc='right margin',
                    frameon=False, legend_fontsize='medium', legend_fontoutline=2, legend_fontweight='normal',
                    use_raw=False, show=False, title=args.batch)
        ax.legend(
            frameon=False,
            loc='center left',
            bbox_to_anchor=(1, 0.5),
            ncol=1,
            fontsize=16
        )
        plt.title(args.batch, fontsize=24)
        pp.savefig(fig, bbox_inches='tight')
        plt.close()
    return adata

def annotation(adata, groups):
    if not adata.raw:
        print('Skip annotation since the process needs the expression data for all genes.')
    else:
        print('Cell type annotation...')
        
        # Export csv used by scMatch
        mat = np.zeros((len(adata.raw.var.index), len(groups)), dtype=float)
        for group in groups:
            mat[: , int(group)] = adata.raw.X[adata.obs['louvain']==group].mean(axis=0)
        dat = pd.DataFrame(mat, index = adata.raw.var.index, columns = groups)
        dat.to_csv(os.path.join(args.output, 'cluster_mean_exp.csv'))
        
        os.system('python /opt/scMatch/scMatch.py --refDS /opt/scMatch/refDB/FANTOM5 \
                --dFormat csv --refType {} --testType {} --testDS {} --coreNum {}'.format(
                args.species, args.species, os.path.join( args.output, 'cluster_mean_exp.csv'), args.cpus))
        
        # Cell annotation result
        scMatch_cluster_df = pd.read_csv(os.path.join(args.output, 'cluster_mean_exp') + '/annotation_result_keep_all_genes/{}_Spearman_top_ann.csv'.format(args.species))
        scMatch_cluster_names = [group + " " + scMatch_cluster_df.loc[scMatch_cluster_df['cell']==int(group)]\
                                ['cell type'].tolist()[0] for group in groups]
        adata.obs['cell_type'] = adata.obs['louvain'].cat.rename_categories(scMatch_cluster_names)
        scMatch_candidate_df = pd.read_excel(os.path.join(args.output, 'cluster_mean_exp') + '/annotation_result_keep_all_genes/{}_Spearman.xlsx'.format(args.species), skiprows=4, header=None, index_col=0)
        for i in range(len(scMatch_candidate_df.columns)):
            if i%2 == 0:
                scMatch_candidate_df.iloc[:, i] = [x.split(',',1)[0].split(':',1)[0] for x in scMatch_candidate_df.iloc[:, i]]
        dict_candidates = {}
        for i in range(int(len(scMatch_candidate_df.columns)/2)):
            candidates = list(dict.fromkeys(scMatch_candidate_df.iloc[:5, i*2]))
            idx = 5
            while len(candidates) < 5:
                cell = scMatch_candidate_df.iloc[idx, i*2]
                if not cell in candidates:
                    candidates.append(cell)
                idx += 1
            dict_candidates[str(i)] = candidates
        df_candidate = pd.DataFrame(dict_candidates).T.reset_index().rename(columns={'index':'cluster'})
        del scMatch_candidate_df

        fig, ax = plt.subplots(dpi=300, figsize=(12,12))
        sc.pl.umap(adata, color='louvain', ax=ax, legend_loc='on data', frameon=False, legend_fontsize='small', legend_fontoutline=2, legend_fontweight='normal', show=False)
        labels = sorted(adata.obs['cell_type'].unique(), key=lambda y: int(y.split(' ',1)[0]))
        for label in labels:
            ind = int(label.split(' ', 1)[0])
            ax.scatter([], [], c=adata.uns['louvain_colors'][ind], label=label)
            ax.legend(
                frameon=False,
                loc='center left',
                bbox_to_anchor=(1, 0.5),
                ncol=1,
                fontsize=16
            )
        fig.tight_layout()
        pp.savefig(fig, bbox_inches='tight')
        plt.close()

        fig, ax = plt.subplots(dpi=300)
        # hide axes
        fig.patch.set_visible(False)
        ax.axis('off')
        ax.axis('tight')
        table = ax.table(cellText=df_candidate.values, colLabels=df_candidate.columns, loc='center')
        # table.auto_set_font_size(False)
        # table.set_fontsize(8)
        table.auto_set_column_width(col=list(range(len(df_candidate.columns))))
        table.scale(1, 1.8)

        for cell in table._cells:
            if cell[0] == 0:
                table._cells[cell].set_color('lightblue')
                table._cells[cell].set_height(.05)
        
        ax.set_title('Top 5 annotation', y=1.08)
        fig.tight_layout()
        pp.savefig(fig, bbox_inches='tight')
        plt.close('all')

    return adata

def findDEG(adata, groups):
    # Finding differentially expressed genes
    print('Finding Differentially Expressed Genes...')
    method = "t-test"
    if adata.raw:
        if adata.raw.X.max() >= 50:
            # haven't been log-transform
            adata_raw = adata.raw.to_adata()
            adata_raw_bk = adata_raw.copy()
            adata_raw.X = np.log1p(adata_raw.X)
            adata.raw = adata_raw
            sc.tl.rank_genes_groups(adata, 'louvain', method=method, pts=True, use_raw=True)
            adata.raw = adata_raw_bk
    else:
        sc.tl.rank_genes_groups(adata, 'louvain', method=method, pts=True)


    # cluster DEGs
    result = adata.uns['rank_genes_groups']
    dat = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'logfoldchanges','scores','pvals']})
    dat.to_csv(os.path.join(args.output, 'cluster_DEGs.csv'))

    return adata

def survivalAnalysis(adata, clinicalpath, gencode):
    print('Survival Analysis...')
    if args.not_treated:
        treatment_status = 'no'
    else:
        treatment_status = ''
    project_ids = [f.rsplit('.csv',1)[0] for f in listdir(args.tcga) if isfile(join(args.tcga, f)) and f.startswith('TCGA')]
    if args.id:
        if not args.id in project_ids:
            print(f"Cannot find the TCGA file for the specified project_id: {args.id}\nCandidates:{project_ids}")
            return adata
    gencode_table = pd.read_csv(gencode, names=['id','name'])
    clinical_df_all = pd.read_csv(clinicalpath, sep='\t', index_col=0)
    cph_selected_cols = ['project', 'cell type', 'coef', 'exp(coef)', 'se(coef)', 'coef lower 95%', 'coef upper 95%','exp(coef) lower 95%', 'exp(coef) upper 95%', 'z', 'p', '-log2(p)']
    df_hazard = pd.DataFrame(columns = cph_selected_cols)
    if not args.id: # user didn't specify tumor type
        run_project_ids = project_ids
    else:
        run_project_ids = [args.id]

    for project_id in run_project_ids:
        print(project_id)
        
        # Read bulk GEP
        bk_gep_name = getBulkProfile(f'{args.tcga}/{project_id}.csv.gz', gencode_table)
        
        # Extract DEGs (masked: mitochondrial genes and genes that not exist in bulk GEPs)
        dict_deg = {}
        all_degs = []
        for celltype in adata.obs['louvain'].unique().tolist():
            degs = [x for x in adata.uns['rank_genes_groups']['names'][celltype] if not (x.startswith('MT-') or not x in bk_gep_name.index)][:20]
            dict_deg[celltype] = degs
            all_degs.extend(degs)
        all_degs = list(set(all_degs))

        # Keep expression of mutual genes in bulk GEP
        bk_gep_name = bk_gep_name.loc[all_degs,:]
        
        # Add clinical info
        if treatment_status == 'no':
            clinical_df = clinical_df_all[(clinical_df_all['project_id'] == project_id) & (clinical_df_all['isTreated'] == treatment_status)]
        else:
            clinical_df = clinical_df_all[clinical_df_all['project_id'] == project_id]
        
        # Remove treated/untreated samples
        bk_gep_name = bk_gep_name.loc[:, [x for x in clinical_df['case_submitter_id'].tolist() if x in bk_gep_name.columns]]
        if bk_gep_name.shape[1] < 10:
            print(f'Skipped: the number of {project_id} patients (treated: {treatment_status}) is less than 10.')
            continue
        
        clinical_df['_event_time_'] = clinical_df['days_to_death']
        for ind in clinical_df.index:
            if not clinical_df.loc[ind, '_event_time_'].isnumeric():
                clinical_df.loc[ind, '_event_time_'] = clinical_df.loc[ind, 'days_to_last_follow_up']
        
        # Build the score table
        # Find patients that show high/low expression for each cell type
        dict_celltype = getSpecCellDict(bk_gep_name, dict_deg)
            # Compute and draw plot
        df_new_hazard = pd.DataFrame()
        try:
            df_new_hazard = drawSurvivalPlot(dict_celltype, clinical_df, project_id)
            df_new_hazard.iloc[:,2:] = df_new_hazard.iloc[:,2:].round(2)
            #df_new_hazard = df_new_hazard.iloc[:, [0, 1, 2, 3, 7, 8, 9, 10]]
            fig, ax = plt.subplots()
            # hide axes
            fig.patch.set_visible(False)
            ax.axis('off')
            ax.axis('tight')
            table = ax.table(cellText=df_new_hazard.iloc[:,1:].values, colLabels=df_new_hazard.columns[1:], loc='center')
            table.auto_set_font_size(False)
            table.auto_set_column_width(col=list(range(len(df_new_hazard.columns[1:]))))
            table.set_fontsize(8)
            for cell in table._cells:
                if cell[0] == 0:
                    table._cells[cell].set_fontsize(6)
                    table._cells[cell].set_color('lightblue')
                    table._cells[cell].set_height(.05)
            
            ax.set_title(project_id, y=1.08)
            fig.tight_layout()
            pp.savefig(fig, bbox_inches='tight')
            plt.close('all')
        except Exception as e:
            print(e)
            print(f'Survival analysis failed for {project_id}. The reason might be that the abundances, when conditioned on the presence or absence of death events, have very low variance and thus fail to converge.')
            continue
        else:
            df_hazard = pd.concat([df_hazard, df_new_hazard], ignore_index=True, join="inner")

    adata.uns['survival_analysis'] = df_hazard
    df_hazard.to_csv(f'{args.output}/HR_{treatment_status}.csv')
    return adata
## Parse command-line arguments
# process arguments
parser = argparse.ArgumentParser(description="scRNA-seq data analysis")

parser.add_argument("-i", "--input", required=True, help="path to input 10x directory or CSV file")
parser.add_argument("-f", "--format", default='10x', help="input format, 10x (default) | csv | h5ad (Anndata object for subclustering with --clusters CLUSTERS)")
parser.add_argument("-o", "--output", default='./', help="path to output directory, default='./'")
parser.add_argument("-r", "--resolution", type=float, default=0.6, help="resolution for clustering, default=0.6")
parser.add_argument("--impute", action="store_true", help="do imputation. default: no")
parser.add_argument("--auto-resolution", action="store_true", help="automatically determine resolution for clustering")
parser.add_argument("-m", "--metadata", default=None, help="path to metadata CSV file for batch correction (index as input in first column)")
parser.add_argument("-b", "--batch", default=None, help="column in metadata (or adata.obs) for batch correction, e.g. 'PatientID'")
parser.add_argument("-c", "--clusters", default=None, help="perform single cell analysis only on specified clusters, e.g. '1,3,8,9'")
parser.add_argument("--cname", default='louvain', help="which variable should be used when selecting clusters; required when clusters are provided. Default: 'louvain'")
parser.add_argument("--GEP", action="store_true", help=' generate Gene Expression Profile file.')
parser.add_argument("--annotation", action="store_true", help="perform cell type annotation")
parser.add_argument("--species", default="human", help="sample species. Options: human (default) | mouse")
parser.add_argument("--gsea", action="store_true", help="perform gene set enrichment analysis (GSEA)")
parser.add_argument("--cpus", default=1, type=int, help="number of CPU used for auto-resolution and annotation, default=1")
parser.add_argument("--survival", action="store_true", help="perform survival analysis")
parser.add_argument("--tcga", default='/scDrug/data/TCGA/', help="path to TCGA data")
parser.add_argument("--id", default=None, help='Specify TCGA project id in the format "TCGA-xxxx", e.g., "TCGA-LIHC"')
parser.add_argument("--prefix", default=None, help='Any prefix before matrix.mtx, genes.tsv and barcodes.tsv.')
parser.add_argument("--not_treated", action="store_true", help='only consider untreated samples from TCGA for survival analysis.')

global args
args = parser.parse_args()

# check format, input and clusters
if not os.path.exists(args.input):
    sys.exit("The input path does not exist.")
if args.format == 'csv':
    if args.input[-4:] != '.csv':
        sys.exit("The input file is not a CSV file.")
elif args.format == '10x':
    if args.prefix:
        prefixes = f'{args.input}/{args.prefix}'
    else:
        prefixes = f'{args.input}/'
    if not os.path.exists(f'{prefixes}matrix.mtx') and not os.path.exists(f'{prefixes}matrix.mtx.gz'):
        sys.exit("Cannot find 'matrix.mtx' file in the input directory.")
    if not os.path.exists(f'{prefixes}genes.tsv') and not os.path.exists(f'{prefixes}genes.tsv.gz'):
        if not os.path.exists(f'{prefixes}features.tsv') and not os.path.exists(f'{prefixes}features.tsv.gz'):
            sys.exit("Cannot find 'genes.tsv' or 'features.tsv' file in the input directory.")
    if not os.path.exists(f'{prefixes}barcodes.tsv') and not os.path.exists(f'{prefixes}barcodes.tsv.gz'):
        sys.exit("Cannot find 'barcodes.tsv' file in the input directory.")
elif args.format == 'h5ad':
    if args.input[-5:] != '.h5ad':
        sys.exit("The input file is not a h5ad file.")
else:
     sys.exit("The format can only be '10x' or 'csv'.")

# check output
if not os.path.isdir(args.output):
    sys.exit("The output directory does not exist.")

# check option
if args.survival:
    if not os.path.isdir(args.tcga):
            sys.exit("The path to TCGA files does not exist.")
    if not os.path.isfile(f'{args.tcga}/clinical.tsv'):
        sys.exit("The TCGA clinical file does not exist.")
    else:
        clinicalpath = f'{args.tcga}/clinical.tsv'
    if not os.path.isfile(f'{args.tcga}/gencode.v22.annotation.id.name.gtf'):
        sys.exit("The gencode file for id conversion does not exist.")
    else:
        gencode = f'{args.tcga}/gencode.v22.annotation.id.name.gtf'

# read input file
global pdfname, pp, results_file

pdfname = f'{args.output}/results.pdf'
pp = PdfPages(pdfname)
results_file = os.path.join(args.output, 'scanpyobj.h5ad')

# Main process
adata = readFile()
adata = preprocessing(adata)
if args.batch:
    adata = batchCorrect(adata)
adata = clustering(adata)
groups = sorted(adata.obs['louvain'].unique(), key=int)
if args.GEP: 
    writeGEP(adata.raw.to_adata())
if args.annotation:
    if args.species not in ['human', 'mouse']:
        print('Invalid species: {}. Use human instead.'.format(args.species))
        args.species = 'human'
    adata =  annotation(adata, groups)
adata = findDEG(adata, groups)
if args.gsea:
    runGSEAPY(adata)
if args.survival:
    adata = survivalAnalysis(adata, clinicalpath, gencode)
pp.close()
if args.clusters:
    results_file = '{}.sub.h5ad'.format(results_file.rsplit('.',1)[0])
if 'features' in adata.var.columns:
    adata.var.index = adata.var['features']
adata.write(results_file)
