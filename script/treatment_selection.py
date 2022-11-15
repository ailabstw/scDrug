import argparse
import csv
import itertools
import math
import os
import pickle
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt

plt.rcParams['figure.dpi'] = 300
pd.options.mode.chained_assignment = None

global threshold
global con_threshold
global pp

class Drug:
    
    def __init__(self, name, dose, time, inst, effect):
        self.name = name
        self.dose = dose
        self.time = time
        self.id = '{}_{:.3f}_{}'.format(self.name, self.dose, self.time)
        self.inst = inst
        self.effect = np.array(effect)
        self.sum_effect = self.effect.sum()

    def __eq__(self, other):
        return self.sum_effect == other.sum_effect and self.dose == other.dose

    def __lt__(self, other):
        return (self.sum_effect > other.sum_effect)

    def __gt__(self, other):
        return (self.sum_effect < other.sum_effect)

    def __str__(self):
        return '{} {}uM {}h'.format(self.name, self.dose, self.time)


def append_df_effect(df, drug_id, result_effect):
    n_killed_all = len([x for x in result_effect if x <= threshold])
    df.loc[drug_id] = result_effect + [n_killed_all]
    return df


def read_metadata(filename, celltype=None):
    if '.csv' in filename:
        df_metadata = pd.read_csv(filename, sep=',', index_col=0)
    else:
        if '.gz' in filename:
            df_metadata = pd.read_csv(filename, sep='\t', index_col=0, compression='gzip')
        else:
            df_metadata = pd.read_csv(filename, sep='\t', index_col=0)

    if celltype:
        df_metadata = df_metadata[df_metadata['cell_id'] == celltype]
        if(df_metadata.shape[0] < 1):
            sys.exit("Error: Perturbation for the selected cell type is less than 1.")
    print('{} perturbations. '.format(df_metadata.shape[0]))
    grouped_plate = df_metadata.groupby(['det_plate'])[['pert_iname','pert_dose','pert_dose_unit','pert_time','pert_time_unit']]
    list_plate_df = [grouped_plate.get_group(x) for x in grouped_plate.groups]
    return list_plate_df


def add_consistency_info(df, DICT_DRUG):
    for cluster_id in df.columns[:-1]:
        colname = 'con_{}'.format(cluster_id)
        df[colname] = 0
        for drug_id in df.index:
            if(df.loc[drug_id, cluster_id] <= threshold):
                drug = DICT_DRUG[drug_id]
                keys = [x for x in DICT_DRUG.keys() if DICT_DRUG[x].name == drug.name and DICT_DRUG[x].dose >= drug.dose and DICT_DRUG[x].time >= drug.time]
                if(all(df.loc[keys,cluster_id] <= con_threshold)):
                    df.loc[drug_id, colname] = 1
    return df


def cal_effect_consistency(df):
    mid_col = int(len(df.columns)/2)
    df_eff = df[df.columns[:mid_col+1]]
    df_con = df[df.columns[mid_col+1:]]
    for i, j in zip(df_eff.columns.to_list()[:-1], df_con.columns.to_list()):
        tmp = df_eff[i] * df_con[j]
        df_eff[i] = tmp.values
    df_eff['kill_all_count'] = df_eff.apply(lambda row: len([x for x in row if x <= threshold]), axis=1)
    return df_eff


def update_df_effect(df, removed_clusters = []):
    #print('\nremoved_clusters: ', removed_clusters)
    df = df.drop(columns = removed_clusters)
    df['kill_all_count'] = df.apply(lambda row: len([x for x in row if x <= threshold]), axis=1)
    return df

def choose_strongest(drug_ids):
    if len(drug_ids) <= 1 :
        return drug_ids
    else:
        drug_list = [DICT_DRUG[x] for x in drug_ids]
        least_survival = max(drug_list).sum_effect
        return [x.id for x in drug_list if x.sum_effect == least_survival]

def select_candidate_drugs(df, _max):
    if(_max <= 0): return []
    #print('\n1st step: kill same amount clusters ( %d clusters) : ' %  _max)
    same_ability_drugs = df[df['kill_all_count'].values == _max].index.to_list()
    # choose the one with strongest killing ability
    same_ability_drugs = choose_strongest(same_ability_drugs)
    
    #same_ability_drugs = sorted(same_ability_drugs, reverse=True, key=lambda d: (d.split('_',2)[0], float(d.split('_',2)[1])))
    # Only the one with min concentration will be kept in candidates among the drugs with same name
    dict_filtered_drugs = {}
    for d in same_ability_drugs:
        name = d.split('_',2)[0]
        if not name in dict_filtered_drugs.keys():
            dict_filtered_drugs[name] = d
    return list(dict_filtered_drugs.values())


def find_drug(df, solution=[], LIST_SOLUTION=[]):
    # All cell types were killed
    if len(df.columns) == 1 or df.min().min() > threshold :
        LIST_SOLUTION.append(sorted(solution))
    else :
        candidates = select_candidate_drugs(df, df['kill_all_count'].values.max())
        # sort candidates by its concentration
        #candidates = sorted(candidates, reverse=False, key=lambda d: float(d.split('_',2)[1]))[:len([x for x in candidates if float(x.split('_',2)[1]) == float(candidates[0].split('_',2)[1])])]
        for drug in candidates:
            solution.append(drug)
            killed_clusters = [x for x in df.columns if df.loc[drug,x] <= threshold]
            df_tmp = update_df_effect(df, killed_clusters)
            find_drug(df_tmp, solution, LIST_SOLUTION)      
            solution.pop()



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Select treatment combination from the LINCS L1000 database.")

    parser.add_argument("-i", "--input", required=True, help="CIBERSORTx output file.")
    parser.add_argument("-o", "--outdir", default='./', help="path to output directory, default='./'")
    parser.add_argument("-t", "--threshold", default=-0.9, type=float, help="Sensitivity threshold. Range: [-1,0), default=-0.9")
    parser.add_argument("-c", "--con_threshold", default=-0.75, type=float, help="Consistency threshold. Range: [-1,0), default=-0.75")
    parser.add_argument("--celltype", required=True, type=str, help="Same as the cell type for decomposition. Options: A375 | A549 | HEPG2 | HT29 | MCF7 | PC3 | YAPC")
    parser.add_argument("--metadata", default='./GSE70138_Broad_LINCS_inst_info_2017-03-06.txt', help="the L1000 instance info file, e.g., 'GSE70138_Broad_LINCS_inst_info_2017-03-06.txt'")

    args = parser.parse_args()
    cell_types = ['A375','A549','HCC515','HEPG2','HT29','MCF7','PC3','YAPC']
    
    # check arguments
    if not os.path.isfile(args.input):
        sys.exit("The input file does not exist.")
    if not os.path.isdir(args.outdir):
        sys.exit("The output directory does not exist.")
    if args.threshold >= 0 and args.threshold < -1:
        sys.exit("Unacceptable sensitivity threshold. The threshold range is [-1,0).")
    if args.con_threshold >= 0 and args.con_threshold < -1:
        sys.exit("Unacceptable consistency threshold. The threshold range is [-1,0).")
    if not args.celltype in cell_types:
        sys.exit("The provided cell type does not exist in the LINCS L1000 database.")
    if not os.path.isfile(args.metadata):
        if not os.path.isfile(args.metadata+'.gz'):
            sys.exit("The metadata file for LINCS L1000 does not exist.")


    print('--------Preprocessing--------')

    threshold = args.threshold
    con_threshold = args.con_threshold
    print('sensitivity threshold: {}, consistency threshold: {} '.format(threshold, con_threshold))
    
    # read metadata
    print('Reading the metadata file and the decomposition results...')
    list_df_metadata = read_metadata(args.metadata, celltype = args.celltype)

    # read decomposition results
    if( '.csv' in args.input):
        df_comp = pd.read_csv(args.input, index_col=0)
    else:
        df_comp = pd.read_csv(args.input, index_col=0, sep='\t')
    df_comp = df_comp.reindex(sorted(df_comp.columns), axis=1)
    df_percent = pd.DataFrame(data=0, index = df_comp.index, columns = sorted(df_comp.columns[:-3].to_list()))
    n_clusters = len(df_comp.columns) - 3 
    
    col_names = df_comp.columns[:-3].to_list()+ ['kill_all_count']

    df_effect = pd.DataFrame(columns=col_names)
    DICT_DRUG_PRE = {}

    for df_metadata in list_df_metadata:
        df_result = pd.merge(df_metadata, df_percent, how='left', left_index=True, right_index=True, validate='one_to_one')
        
        ctrl_inst_ids = df_result[df_result['pert_dose']<0].index
        pert_inst_ids = df_result[df_result['pert_dose']>0].index
        df_result = df_result[df_result['pert_dose']>0]

        # average ctrl results
        df_ctrl_composition = df_comp.loc[ctrl_inst_ids,df_comp.columns[:n_clusters]]
        ctrl_composition = df_ctrl_composition.mean()
        ctrl_composition = ctrl_composition/np.sum(ctrl_composition.values)
        
        # fill the result dataframe
        for perturbation in df_result.index:
            drug_composition = df_comp.loc[perturbation,df_comp.columns[:n_clusters]]
            change_list = (drug_composition-ctrl_composition)/(drug_composition+ctrl_composition)

            # fill nan (occurs when 0/0) with 2 
            change_list = [ 2 if math.isnan(x) else x for x in change_list]
            df_result.loc[perturbation,df_result.columns[-n_clusters:]] = change_list

            drug = Drug(df_metadata.loc[perturbation,'pert_iname'],\
                        df_metadata.loc[perturbation,'pert_dose'],\
                        df_metadata.loc[perturbation,'pert_time'],\
                        perturbation,\
                        change_list)
            drug_id = '{}_{:.3f}_{}'.format(\
                df_metadata.loc[perturbation,'pert_iname'],\
                df_metadata.loc[perturbation,'pert_dose'],\
                df_metadata.loc[perturbation,'pert_time'])

            if drug_id not in DICT_DRUG_PRE:
                DICT_DRUG_PRE[drug_id] = [drug]
            else:
                DICT_DRUG_PRE[drug_id].append(drug)

    #print('DICT_DRUG_PRE: ', len(DICT_DRUG_PRE))

    print('Calculating drug effects...')
    # average replicates
    DICT_DRUG = {}
    for drug in DICT_DRUG_PRE:
        # only keep those have replicates
        if (len(DICT_DRUG_PRE[drug]) <=2):
            continue
        avg_effect = np.array([0.0]*n_clusters)
        inst_list = []
        counts = np.array([0.0]*n_clusters)
        
        for x in DICT_DRUG_PRE[drug]:
            for i in range(0,len(x.effect)):
                if(x.effect[i] < 2): # The subpopulation existed in the control samples
                    avg_effect[i] += x.effect[i]
                    counts[i] += 1
        for i in range(0, len(avg_effect)):
            if counts[i] >= 3:
                avg_effect[i] /= counts[i]
            else:
                avg_effect[i] = 0
        
        x = DICT_DRUG_PRE[drug][0]
        avg_drug = Drug(x.name, x.dose, x.time, inst_list, avg_effect)
        DICT_DRUG[avg_drug.id] = avg_drug
        df_effect = append_df_effect(df_effect, avg_drug.id, avg_effect.tolist())
    df_effect = add_consistency_info(df_effect, DICT_DRUG)
    


    # store output
    df_effect_bk = df_effect.iloc[:, :int(len(df_effect.columns)/2)]
    df_effect_bk = df_effect_bk.reindex(sorted(df_effect_bk.columns, key=int), axis=1)
    outputname = args.input.rsplit('/',1)[1].rsplit('.')[0]
    # fh = open(file = '{}/{}_t{}_ct{}_df_effect.pickle'.format(args.outdir, outputname, threshold, con_threshold), mode='wb')
    # pickle.dump(df_effect, fh) 
    # fh.close()
    # fh = open(file = '{}/{}_t{}_ct{}_DICT_DRUG_PRE.pickle'.format(args.outdir, outputname, threshold, con_threshold), mode='wb')
    # pickle.dump(DICT_DRUG_PRE, fh)
    # fh = open(file = '{}/{}_t{}_ct{}_DICT_DRUG.pickle'.format(args.outdir, outputname, threshold, con_threshold), mode='wb')
    # pickle.dump(DICT_DRUG, fh)
    # fh.close()


    print('\n--------Subpopulation Analysis--------')
    # subpopulation analysis
    df_effect = cal_effect_consistency(df_effect)
    grouped = df_effect.groupby('kill_all_count')
    
    for name, group in grouped:
        print('{} perturbations (e.g., {}) can kill {} cluster(s).'.format(group.shape[0], group.index[0], name))
    
    # subpopulation analysis
    resistant_clusters = []
    for cluster in sorted(df_effect.columns[:-1], key=int):
        pert = df_effect[df_effect[cluster] <= threshold].index.to_list()
        if (len(pert) > 0):
            print('Cluster {} can be killed by {} perturbations.'.format(cluster, len(pert)))
            print('    peturbation with best efficacy: {} -> {:.2f}'.format(df_effect[cluster].idxmin(), df_effect[cluster].min()))
        else:
            resistant_clusters.append(cluster)
    if(len(resistant_clusters)>0):
        print('\nClusters cannot be killed: ', sorted(resistant_clusters, key=int))
        print('\n')


    print('--------Treatment Selection--------')
    LIST_SOLUTION = []

    # find a cocktail therapy
    find_drug(df_effect, LIST_SOLUTION=LIST_SOLUTION)
    
  
    LIST_SOLUTION = sorted(LIST_SOLUTION)
    list_result = list(LIST_SOLUTION for LIST_SOLUTION,_ in itertools.groupby(LIST_SOLUTION))

    print('--------Generating Figures--------')
    def draw_lincs_heatmap(df, drugs, name=""):
        if len(drugs[0].split('_')) > 1 :
            selected_index = drugs
        else:
            selected_index = [x for x in df.index if x.split('_',1)[0] in drugs]
        if len(selected_index) < 1:
            print('Selected drugs were not used in the treatment selection:')
            print(*drugs)
        else:
            subdf = df.loc[selected_index,:]
            subdf = subdf.assign(sum=subdf.sum(axis=1)).sort_values(by='sum', ascending=True).drop(columns='sum')
            for i in subdf.index:
                for j in subdf.columns:
                    if(subdf.loc[i,j] > 0): subdf.loc[i,j] = 0
            sns.set(rc={'figure.figsize':(12, int(len(selected_index))/2)})
            subdf.index = [x.rsplit('_',1)[0] for x in subdf.index]
            total_row = pd.Series(subdf.min(),name='combined treatment')
            subdf = subdf.append(total_row)
            ax = sns.heatmap(subdf, center=0, cmap='RdBu', vmin=-1, vmax=0, \
                    linewidths=0.5, linecolor='lightgrey', cbar=True, cbar_kws={'label': 'cell survival rate'})
            for _, spine in ax.spines.items():
                spine.set_visible(True)
                spine.set_color('lightgrey')
            lab = ax.get_yticklabels()[-1]
            lab.set_weight('bold')
            lab.set_color('red')
            ax.set(xlabel='Cluster')

            c_bar = ax.collections[0].colorbar
            c_bar.set_ticks([threshold, -0.5, 0])
            c_bar.set_ticklabels([str(threshold), str(-0.5), '\u2265 0'])
            pp.savefig(bbox_inches='tight')
            plt.close()
            subdf.to_csv('{}'.format(os.path.join(args.outdir, 'treatment_combinations.csv')))

    def draw_lincs_cons_plot(index_list, df_effect, df_all):
        index_list = [x for x in index_list if x in df_effect.index]
        index_list = index_list[::-1]
        df = df_effect.loc[index_list,:]
        df_all = df_all.append(df)
        df = df.T
        df.plot(kind='bar', ylim=(-1.0,1.0), rot=0, colormap='summer_r', width=0.7, figsize=(15,5), yticks=[-1.0,-0.5,0,0.5,1.0], xlabel='cluster', ylabel='survival rate')
        plt.legend(labels=[x.split('_',2)[1]+' \u03BCM' for x in df.columns], loc='center left', bbox_to_anchor=(1.0, 0.5))
        plt.title(index_list[0].split('_',1)[0], fontweight='bold')
        plt.axhline(y=0, color='black', linestyle='-', lw=0.8)
        plt.axhline(y=threshold, color='red', linestyle='dotted', lw=0.8)
        n_c = len(df_effect.columns)
        plt.text(n_c-0.3, threshold+0.01, 'threshold={}'.format(threshold), fontsize=8, color='red')
        plt.axhline(y=con_threshold, color='blue', linestyle='dotted', lw=0.8)
        plt.text(n_c-0.3, con_threshold+0.01, 'con. threshold={}'.format(con_threshold), fontsize=8, color='blue')
        pp.savefig(bbox_inches='tight')
        plt.close()
        return df_all

    #print(list_result)
    keywords_long = [item for sublist in list_result for item in sublist]

    keywords = list(set([x.split('_',1)[0] for x in keywords_long if isinstance(x, str) ]))
    #print('{} unique drugs.'.format(len(keywords)))

    DICT_INDEX = {}

    pdfname = os.path.join(args.outdir, 'treatment_combinations.pdf')
    pp = PdfPages(pdfname)
    df_all = pd.DataFrame(columns = df_effect_bk.columns)
    for drug_id in sorted(keywords):
        DICT_INDEX[drug_id] = [x for x in DICT_DRUG_PRE.keys() if drug_id == x.split('_',1)[0]]
        df_all = draw_lincs_cons_plot(DICT_INDEX[drug_id], df_effect_bk, df_all)
    
    # consistency table
    fig, ax = plt.subplots()
    # hide axes
    fig.patch.set_visible(False)
    ax.axis('off')
    ax.axis('tight')
    df_cons = df_all.loc[keywords_long,:].round(2).reset_index()
    table = ax.table(cellText=df_cons.values, colLabels=df_cons.columns, loc='center')
    table.auto_set_font_size(False)
    table.auto_set_column_width(col=list(range(len(df_cons.columns))))
    table.set_fontsize(8)
    for cell in table._cells:
        if cell[0] == 0:
            table._cells[cell].set_fontsize(6)
            table._cells[cell].set_color('lightblue')
            table._cells[cell].set_height(.05)
    
    ax.set_title(args.celltype)
    pp.savefig(fig, bbox_inches='tight')
    plt.close('all')

    for i, comb in enumerate(list_result):
        draw_lincs_heatmap(df_all, [x.strip() for x in comb], name=str(i))

    pp.close()
    print('Figures are stored in \'{}\'.'.format(pdfname))





