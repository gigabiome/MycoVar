
########################################################################################################
## Collection of functions used by three notebooks - phylo, lineage and resisto - to perform analysis ##
########################################################################################################

######### SETUP #########

def renameFiles(files,filetype):
    import os,sys
    from pathlib import Path
    
    for old_path in files:
        dir_name = os.path.dirname(old_path)
        base_name = os.path.basename(old_path)
        if os.path.splitext(base_name)[1] == '.csv':
            new_base = base_name.split('_')[0] + '.'+filetype
            new_path = os.path.join(dir_name, new_base)
            os.rename(old_path, new_path)
            
def checkCSV(csv_files):
    import csv
    false_csvs = []
    for file in csv_files:
        try:
            with open(file, newline='', encoding="utf-8") as f:
                reader = csv.reader(f)
                for _ in range(5):  
                    next(reader, None)
        except Exception as e:
            print(f"❌ {file} is not a valid CSV: {e}")
            false_csvs.append(file)
    return false_csvs
    
# check all input files and rename
def setupInput(input_dir):
    import os,sys
    from pathlib import Path
    from itertools import chain
    
    false_csvs = []
    sub_dirs = ['snp','cov','stat']
    for sub_dir in sub_dirs:
        directory = Path(input_dir+sub_dir)
        files = [f for f in directory.iterdir() if f.is_file()]
        renameFiles(files,sub_dir)
        false_csvs.append(checkCSV(files))
    
    false_csvs = list(chain.from_iterable(false_csvs))
    if len(false_csvs) > 0:
        return(print('The following files are not valid CSV files: '+str('\t'.join(false_csvs))))
    else:
        return(print('✅ Input files successfully renamed and all files valid CSV files\n'))
    
def setupOutput(output_dir,analysis_type):
    import os
    import shutil
    path = output_dir + analysis_type
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)
    
    if analysis_type == 'resisto':
        os.makedirs(path+'/sample_resisto/')
        return(print('✅ Resisto output directories created\n'))
    elif analysis_type == 'phylo':
        os.makedirs(path+'/iTol/')
        return(print('✅ Phylo output directories created\n'))
    
    
######### PHYLOGENETICS #########

def setupPhylo(BIN):
    import pandas as pd
    import plotly.express as px
    
    # drug resistance positions
    DB_DR = []
    df = pd.read_csv(BIN+'drdb.txt',sep='\t')
    df = df[df['Reference_Position'] != '-']
    for snp in df['Reference_Position']:
        snp = snp.split('/')
        for i in snp:
            DB_DR.append(float(i))

    # repeat positions and ne genes
    df1 = pd.read_csv(BIN+'repeats.txt',sep='\t')
    df2 = pd.read_csv(BIN+'ne_genes.txt',sep='\t')
    df = pd.concat([df1,df2])
    DB_REP_NE = df.apply(lambda row:list(range(row['Start'],row['End']+1)),axis=1).explode().astype(int).tolist()

    # lineages
    df_lin = pd.read_csv(BIN+'lineages.bed',sep='\t',index_col=1)
    DB_LIN = df_lin.drop(columns=['Chromosome','End'])
    colors = px.colors.qualitative.Dark24
    colors.append(px.colors.qualitative.Light24[1])
    colors = dict(zip(DB_LIN['Lineage_name'].unique(),colors))
    DB_LIN['Lineage_color'] = DB_LIN['Lineage_name'].map(colors)
    print('✅ Databases successfully loaded\n')
    return[DB_DR,DB_REP_NE,DB_LIN]

def mainPhylo(snp_files,setup_params,phylo_list):
    import pandas as pd
    import numpy as np
    
    QUAL,FRB,FREQ,COV,INPUT,OUTPUT,BIN = setup_params
    DB_DR,DB_REP_NE,DB_LIN = phylo_list
    
    ref_dict = dict()
    snp_dfs = dict()
    snp_dfs_filtered = dict()

    for file in snp_files:
        sample = file.replace('.snp','')
        df = pd.read_csv(INPUT+'snp/'+file)
        snp_dfs[sample] = df

        df['Sample'] = sample

        # check if SNV and not other variant
        df = df[df['Type'] == 'SNV']

        # normal QC
        df = df[(df['Frequency'] >= FREQ) & (df['Coverage'] >= COV) & (df['Average quality'] >= QUAL)]

        # forward reverse more than cut-off (keep to 0 unless for specific reasons)
        df = df[(df['Forward/reverse balance'] > FRB)]

        # exclude resistant genes
        df = df[~df['Reference Position'].isin(DB_DR)]

        # exclude repeats
        df = df[~df['Reference Position'].isin(DB_REP_NE)]

        # exclude exons
        df = df[df['Overlapping annotations'].notna()]
        df = df[~df['Amino acid change'].notna()]

        # exlude proximity snps
        arr_sorted = np.sort(df['Reference Position'])
        diff_next = np.diff(arr_sorted)
        too_close = diff_next <= 50
        mask = np.ones_like(arr_sorted, dtype=bool)
        mask[:-1] &= ~too_close
        mask[1:]  &= ~too_close
        df = df[df['Reference Position'].isin(arr_sorted[mask])]

        # add all reference alleles to ref_dict
        for index,row in df.iterrows():
            ref_dict[row['Reference Position']] = row['Reference']
        snp_dfs_filtered[sample] = df

    # order dictionary
    ref_dict = dict(sorted(ref_dict.items()))
    
    # this next section creates alignment files and runs the necessary tools to perform multiple
    # sequence alignment and build the phylogenetic tree
    
    # create fasta for alignment
    fasta_path = buildFasta(OUTPUT,snp_dfs_filtered,ref_dict)
    # use mafft for alignment
    msa_path = buildMSA(OUTPUT,fasta_path)
    print('✅ Multiple sequence alignment completed!\n')
    
    # build tree
    buildTree(OUTPUT,msa_path)
    print('✅ Phylogenetic tree built!\n')
    

def buildFasta(OUTPUT,dfs,ref):
    
    fasta_path = OUTPUT+'phylo/snps.fasta'
    ffasta = open(fasta_path,'w')
    # create reference for tree root
    ffasta.write('>MTB\n'+''.join(ref.values())+'\n')

    # read in from filtered dataframes
    for sample in dfs:
        df = dfs[sample]
        sample_dict = dict(zip(df['Reference Position'],df['Allele']))
        missing_keys = set(ref.keys()) - set(sample_dict.keys())
        for key in missing_keys:
            sample_dict[key] = ref[key]
        sample_dict = dict(sorted(sample_dict.items()))
        seq = ''.join(sample_dict.values())
        ffasta.write('>'+sample+'\n'+str(seq)+'\n')
    ffasta.close()
    return(fasta_path)

def buildMSA(OUTPUT,path):
    import subprocess
    import sys
    import time
    msa_path = OUTPUT+'phylo/aligned.fasta'
    parts = []
    parts.append('mafft')
    parts.append('--auto')
    parts.append('--quiet')
    parts.append(path)
    parts.append('>')
    parts.append(msa_path)
    mcmd = str(' '.join(parts))

    try:
        process = subprocess.Popen(mcmd,shell=True)
        process.wait()
    except:
        IOError('Could not run mafft: check fasta files')
        print("❌ Could not run mafft: check fasta files")
        sys.exit(1)
    time.sleep(1)
    return(msa_path)

def buildTree(OUTPUT,path):
    import subprocess
    import sys
    import time
    tree_path = OUTPUT+'phylo/aligned.fasta'

    parts = []
    parts.append('iqtree2')
    parts.append('-nt')
    parts.append('AUTO')
    parts.append('-m')
    parts.append('MFP')
    parts.append('-quiet')
    parts.append('-s')
    parts.append(tree_path)
    #parts.append('-b')
    #parts.append('100')
    icmd = str(' '.join(parts))
    
    try:
        process = subprocess.Popen(icmd,shell=True)
        process.wait()
    except:
        IOError('Could not run iqtree: check fasta files')
        print("❌ Could not run iqtree: check fasta files")
        sys.exit(1)

######### RESISTOTYPING #########

def setupResisto(bin_path):
    import pandas as pd
    from collections import defaultdict
    # parse databases
    DR_WHO = pd.read_excel(bin_path+'WHO-UCN-TB-2023.5-eng.xlsx',skiprows=2)
    DR_WHO = DR_WHO[['drug','gene','mutation','variant','tier','genomic position',
                    'Present_SOLO_SR','Present_SOLO_R','Present_SOLO_S','FINAL CONFIDENCE GRADING']]
    DR_CTB = pd.read_excel(bin_path+'CTB_Mutation_Catalogue.xlsx')
    DR_CTB.columns = ['CTB_variant','CTB_drug','CTB_FINAL_CONFIDENCE_GRADING']
    
    DR_WHO_promotors = parseWHOForPromotors(DR_WHO)
    
    # create dictionaries to use for resistotyping
    dCTB = defaultdict(dict)
    for index,row in DR_CTB.iterrows():
        dCTB[row['CTB_variant']][row['CTB_drug']] = row['CTB_FINAL_CONFIDENCE_GRADING']
    dWHO = defaultdict(dict)
    for index,row in DR_WHO.iterrows():
        dWHO[row['variant']][row['drug']] = row['FINAL CONFIDENCE GRADING']

    dGene2Drug = defaultdict(list)
    dDrug2Gene = defaultdict(list)

    for index,row in DR_WHO.iterrows():
        dGene2Drug[row['gene']].append(row['drug'])
        dDrug2Gene[row['drug']].append(row['gene'])
    for index,row in DR_CTB.iterrows():
        dGene2Drug[row['CTB_variant'].split('_')[0]].append(row['CTB_drug'])
        dDrug2Gene[row['CTB_drug']].append(row['CTB_variant'].split('_')[0])

    # set lists
    for gene in dGene2Drug:
        dGene2Drug[gene] = list(set(dGene2Drug[gene]))
    for drug in dDrug2Gene:
        dDrug2Gene[drug] = list(set(dDrug2Gene[drug]))
        
    DR_GENES = list(DR_WHO['gene'].unique())
    DR_GENES = list(set(DR_GENES + list(DR_CTB['CTB_variant'].str.split('_').str[0].unique())))
    print('✅ Resistance databases successfully loaded\n')
    return[DR_WHO,DR_CTB,DR_GENES,dCTB,dWHO,dGene2Drug,dDrug2Gene,DR_WHO_promotors]

def parseWHOForPromotors(df):
    import pandas as pd
    import warnings
    warnings.simplefilter(action="ignore", category=pd.errors.SettingWithCopyWarning)
    df = df[df['genomic position'].apply(lambda x: isinstance(x, int))]
    for index,row in df.iterrows():
        if '>' in row['mutation']:
            parts = row['mutation'].split('>')
            ref = parts[0][-1]
            allele = parts[1]
            df.loc[index,'Reference'] = ref
            df.loc[index,'Allele'] = allele
        elif 'ins' in row['mutation'] and 'del' in row['mutation']: 
            # need to deal with insertion-deletions here - check WHO catalog
            pass
        elif 'ins' in row['mutation']:
            parts = row['mutation'].split('ins')
            allele = parts[1]
            df.loc[index,'Reference'] = '-'
            df.loc[index,'Allele'] = allele
        elif 'del' in row['mutation']:
            parts = row['mutation'].split('del')
            ref = parts[1]
            df.loc[index,'Reference'] = ref
            df.loc[index,'Allele'] = '-'
        else:
            # need to deal with dup here - check WHO catalog
            pass
    return df

def mainResistotyping(snp_files,setup_params,resisto_list):
    import pandas as pd
    from collections import defaultdict
    
    WHORules = {'katG':   ['Isoniazid'],
                'pncA':   ['Pyrazinamide'],
                'ddn':    ['Pyrazinamide','Delamanid'],
                'fbiA':   ['Pyrazinamide','Delamanid'],
                'fbiB':   ['Pyrazinamide','Delamanid'],
                'fbiC':   ['Pyrazinamide','Delamanid'],
                'fgd1':   ['Pyrazinamide','Delamanid'],
                'Rv2983': ['Pyrazinamide','Delamanid'],
                'gid':    ['Streptomycin'],
                'ethA':   ['Prothionamide','Ethambutol'],
                'tlyA':   ['Capreomycin'],
                'pepQ' :  ['Bedaquiline','Clofazimine']}
    
    QUAL,FRB,FREQ,COV,INPUT,OUTPUT,BIN = setup_params
    DR_WHO,DR_CTB,DR_GENES,dCTB,dWHO,dGene2Drug,dDrug2Gene,DR_WHO_promotors = resisto_list
    
    dRes = defaultdict(lambda: defaultdict(list))
    columns = []
    for drug in dDrug2Gene:
        columns.append( ('Resistotype',drug,'') )

    queries = ['Coding region change','Amino acid change','Resistance level']
    for drug in dDrug2Gene:
        for gene in dDrug2Gene[drug]:
            for query in queries:
                columns.append( (drug,gene,query) )

    columns = pd.MultiIndex.from_tuples(columns)
    mainOut = pd.DataFrame(columns=columns)
    
    for file in snp_files:
        indexes = []
        sample = file.replace('.snp','')
        dfResSample = pd.DataFrame(columns=['Drug','Gene','Coding region change','Amino acid change','Level'])
        counter = 0

        # read in files and do QC 
        snp = pd.read_csv(INPUT+'snp/'+sample+'.snp') #[SNP_COLUMNS]
        snp = snp[snp['Average quality'] >= QUAL]
        snp = snp[snp['Forward/reverse balance'] > FRB]
        snp = snp[snp['Frequency'] >= FREQ]
        snp = snp[snp['Coverage'] >= COV]

        cov = pd.read_csv(INPUT+'cov/'+sample+'.cov') #[COV_COLUMNS]
        stat = pd.read_csv(INPUT+'stat/'+sample+'.stat') #[STAT_COLUMNS]

        # split column to get gene
        df = snp.copy()
        df = df[df['Overlapping annotations'].notna()]
        df = df[df['Overlapping annotations'].str.contains('Gene')]
        df['Gene'] = df['Overlapping annotations'].str.split('Gene: ').str[1].str.split().str[0]

        # get only DR genes
        df = df[df['Gene'].isin(DR_GENES)]
        # merge snp and cov
        df = pd.merge(df, cov, left_on='Gene', right_on='Name')
        df = df.drop(columns=['Gene'])
        # merge with stat
        for col in stat.columns:
            df[col] = stat.at[0,col]

        # add sample name as a column and move it to first column
        df['Sample'] = sample
        df = df[[df.columns[-1]] + df.columns[:-1].tolist()]
        
         # create resistance dictionary
        
        # add in the grading here for the mutations based on genomic position (promotor regions)
        df2 = df[df['Reference Position'].isin(list(DR_WHO_promotors['genomic position']))]
        for index,row in df2.iterrows():
            position = row['Reference Position']
            subdf = DR_WHO_promotors[DR_WHO_promotors['genomic position'] == position]
            for index2,row2 in subdf.iterrows():
                if row['Allele'] == row2['Allele'] and row['Reference'] == row2['Reference']:
                    dRes[sample][row2['drug']].append(row2['FINAL CONFIDENCE GRADING'])
                    dfResSample.loc[len(dfResSample)] = [row2['drug'],row2['gene'],str(row2['variant']),
                                                         str(row2['variant']),row2['FINAL CONFIDENCE GRADING']]
            
        # check for ungraded mutations -  add ungraded for all genes associated with resistance
        for index,row in df.iterrows():
            if row['Name'] in dGene2Drug:
                for drug in dGene2Drug[row['Name']]:
                    dfResSample.loc[len(dfResSample)] = [drug,str(row['Name']),str(row['Coding region change']),
                                                         str(row['Amino acid change']),'Ungraded']
                    if 'Ungraded' not in dRes[sample][drug]:
                        dRes[sample][drug].append('Ungraded')
                        indexes.append(index)

        # for CTB
        df['tmp1'] = df['Name']+'_'+df['Coding region change'].str.split('c.').str[1] 
        df['tmp2'] = df['Name']+'_'+df['Amino acid change'].str.split('p.').str[1]
        for index,row in df.iterrows():
            if row['tmp1'] in dCTB:
                for drug in dCTB[row['tmp1']]:
                    dRes[sample][drug].append(dCTB[row['tmp1']][drug])
                    indexes.append(index)
                    dfResSample.loc[len(dfResSample)] = [drug,row['Name'],str(row['tmp1']),
                                                         str(row['tmp2']),dCTB[row['tmp1']][drug]]
            if row['tmp2'] in dCTB:
                for drug in dCTB[row['tmp2']]:
                    dRes[sample][drug].append(dCTB[row['tmp2']][drug])
                    indexes.append(index)
                    dfResSample.loc[len(dfResSample)] = [drug,row['Name'],str(row['tmp1']),
                                                         str(row['tmp2']),dCTB[row['tmp2']][drug]]
        df = df.drop(columns=['tmp1','tmp2'])

        # for WHO
        df['tmp1'] = df['Name']+'_'+df['Coding region change'].str.split(':').str[1] 
        df['tmp2'] = df['Name']+'_'+df['Amino acid change'].str.split(':').str[1]
        for index,row in df.iterrows():
            if row['tmp1'] in dWHO:
                for drug in dWHO[row['tmp1']]:
                    dRes[sample][drug].append(dWHO[row['tmp1']][drug])
                    indexes.append(index)
                    dfResSample.loc[len(dfResSample)] = [drug,row['Name'],str(row['tmp1']),
                                                         str(row['tmp2']),dWHO[row['tmp1']][drug]]
            if row['tmp2'] in dWHO:
                for drug in dWHO[row['tmp2']]:
                    dRes[sample][drug].append(dWHO[row['tmp2']][drug])
                    indexes.append(index)
                    dfResSample.loc[len(dfResSample)] = [drug,row['Name'],str(row['tmp1']),
                                                         str(row['tmp2']),dWHO[row['tmp2']][drug]]
        df = df.drop(columns=['tmp1','tmp2'])

        # for WHO expert
        df2 = df[(df['Type'].isin(['Insertion','Deletion'])) | 
                 (df['Amino acid change'].str.contains('*',na=False,regex=False))]
        for index,row in df2.iterrows():
            if row['Name'] in WHORules:
                for drug in WHORules[row['Name']]:
                    dRes[sample][drug].append('Expert_R')
                    dfResSample.loc[len(dfResSample)] = [drug,row['Name'],str(row['Coding region change']),
                                                         str(row['Amino acid change']),'Expert_R']
        if 'Rv0678' in list(df2['Name']):
            if 'mmpL5' in list(df2['Name']):
                pass
            else:
                dRes[sample]['Bedaquiline'].append('Expert_R')
                dRes[sample]['Clofazimine'].append('Expert_R')
                df3 = df2[df2['Name'] == 'Rv0678']
                for index,row in df3.iterrows():
                    dfResSample.loc[len(dfResSample)] = ['Bedaquiline',row['Name'],
                                                         str(row['Coding region change']),
                                                         str(row['Amino acid change']),
                                                         'Expert_R']
                    dfResSample.loc[len(dfResSample)] = ['Clofazimine',row['Name'],
                                                         str(row['Coding region change']),
                                                         str(row['Amino acid change']),
                                                         'Expert_R']

        df = df.loc[list(set(indexes))]
        dfResSample.to_excel(OUTPUT+'resisto/sample_resisto/'+sample+'.xlsx')
        mainOut = resistotyping(dfResSample,mainOut,sample)
        
    dfRes = buildResistogram(dRes)
    for i in dfRes.index:
        for c in dfRes.columns:
            mainOut.at[i,('Resistotype',c,'')] = dfRes.at[i,c]
    mainOut = mainOut.dropna(axis=1, how='all')
    mainOut.to_excel(OUTPUT+'resisto/mainResistotyping.xlsx')
    print('✅ Resistotyping complete and output files written\n')
    return(mainOut)

def buildResistogram(d):
    import pandas as pd
    dfRes = pd.DataFrame()
    for sample in d:
        for drug in d[sample]:
            if '1) Assoc w R' in d[sample][drug]:
                dfRes.at[sample,drug] = 'AWR'
            elif '2) Assoc w R - Interim' in d[sample][drug]:
                dfRes.at[sample,drug] = 'AWR-I'
            elif 'Expert_R' in d[sample][drug]:
                dfRes.at[sample,drug] = 'AWRI-ER'
            elif '3) Uncertain significance' in d[sample][drug]:
                dfRes.at[sample,drug] = 'US'
            elif '5) Not assoc w R' in d[sample][drug]:
                dfRes.at[sample,drug] = 'NAWR'
            elif '4) Not assoc w R - Interim' in d[sample][drug]:
                dfRes.at[sample,drug] = 'NAWR-I'
            elif 'R' in d[sample][drug]:
                dfRes.at[sample,drug] = 'CTB-R'
            elif 'Ungraded' in d[sample][drug]:
                dfRes.at[sample,drug] = 'UNGR'
    dfRes = dfRes.fillna('WT')
    return dfRes

def resistotyping(df,outdf,sample):
    import pandas as pd
    priority = ["1) Assoc w R",
                "2) Assoc w R - Interim",
                "Expert_R",
                "3) Uncertain significance",
                "5) Not assoc w R",
                "4) Not assoc w R - Interim",
                "R",
                "Ungraded"]
    
    for drug in df['Drug'].unique():
        for gene in df['Gene'].unique():
            df2 = df[(df['Drug'] == drug) & (df['Gene'] == gene)]
            for level in priority:
                if (df2['Level'] == level).any():
                    df2 = df2[df2['Level'] == level]
                    break
                
            #print(df2)
            nucl_change = str(' / '.join(df2['Coding region change']))
            aa_change = str(' / '.join(df2['Amino acid change']))
    
            if len(df2.index) > 0:
                level = list(set(df2['Level']))[0]
                # Ungraded will have several mutations (can remove if don't want all the mutations in there)
                if level == 'Ungraded': 
                    outdf.at[sample,(drug,gene,'Coding region change')] = nucl_change
                    outdf.at[sample,(drug,gene,'Amino acid change')] = aa_change
                    outdf.at[sample,(drug,gene,'Resistance level')] = level
                else:
                    outdf.at[sample,(drug,gene,'Coding region change')] = nucl_change
                    outdf.at[sample,(drug,gene,'Amino acid change')] = aa_change
                    outdf.at[sample,(drug,gene,'Resistance level')] = level
            
    return(outdf)















































