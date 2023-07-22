import json, os, glob, sys
import pandas as pd
import numpy as np
import argparse

def search_files(directory, pattern):
    file_paths = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if not file.startswith('.'):
                file_path = os.path.join(root, file)
                if glob.fnmatch.fnmatch(file_path, pattern):
                    file_paths.append(file_path)
    return file_paths

def compare_variants(refdf,querydf,onfields):

    if 'variant' not in refdf.columns or 'variant' not in querydf.columns:
        sys.exit("variant comparison dataframes must have a 'variant' field")

#    if all(f in refdf.columns for f in onfields) or any(onfields not in querydf.columns):
#        sys.exit("variant comparison dataframes do not contain joining fields")

    merged = pd.merge(refdf,
                      querydf,
                      on=onfields,how='outer',suffixes=['.ref','.query'])

    cond = [
        merged['variant.ref'].notna() & merged['variant.query'].notna(),
        merged['variant.ref'].isna() & merged['variant.query'].notna(),
        merged['variant.ref'].notna() & merged['variant.query'].isna()
    ]

    vals = [
        'TP',
        'FP',
        'FN'
    ]

    merged['type'] = pd.Categorical(np.select(cond, vals, default='Unknown'),categories=['TP','FP','FN'])

    out = dict(zip(['num_reference','num_query','TP','FN','FP','unique_reference','unique_query'],
                   [refdf.shape[0],
                   querydf.shape[0],
                   merged['type'].value_counts().to_dict()['TP'],
                   merged['type'].value_counts().to_dict()['FN'],
                   merged['type'].value_counts().to_dict()['FP'],
                   ';'.join(merged[merged['type']=='FN']['variant.ref']),
                   ';'.join(merged[merged['type']=='FP']['variant.query'])]))
    
    if out['num_reference'] > 0:
        out['recall'] = round(out['TP']/(out['TP']+out['FN'])*100,1)
    else:
        out['recall'] = pd.NA

    if out['TP']+out['FP'] > 0:
        out['precision'] = round(out['TP']/(out['TP']+out['FP'])*100,1)
    else:
        out['precision'] = pd.NA

    if out['precision'] is not pd.NA and out['recall'] is not pd.NA:
        out['f1'] = round(2 * (out['precision'] * out['recall']) / (out['precision'] + out['recall']),1)
    else:
        out['f1'] = pd.NA

    return(pd.DataFrame.from_dict([out]))

def compare_cases(row):

    print(row['path.ref'])
    print(row['path.query'])

    refj = pd.read_json(row['path.ref'])
    queryj = pd.read_json(row['path.query'])

    # small variants
    ref_variants = pd.DataFrame() 
    query_variants = pd.DataFrame() 
    if 'columns' in refj['VARIANTS']['PASS'].keys():
        ref_variants = pd.DataFrame(columns=refj['VARIANTS']['PASS']['columns'],data=refj['VARIANTS']['PASS']['data'])
        ref_variants['variant'] = ref_variants.apply(lambda x: x['gene']+':'+x['psyntax']+'['+x['vaf']+']',axis=1)

    if 'columns' in queryj['VARIANTS']['PASS'].keys():
        query_variants = pd.DataFrame(columns=queryj['VARIANTS']['PASS']['columns'],data=queryj['VARIANTS']['PASS']['data'])
        query_variants['variant'] = query_variants.apply(lambda x: x['gene']+':'+x['psyntax']+'['+x['vaf']+']',axis=1)

    out_variants = compare_variants(ref_variants,query_variants,['chrom','pos','ref','alt'])
    out_variants.columns = 'variants_' + out_variants.columns
#    out_variants['case'] = row['case']
#    out_variants.set_index('case')

    # CNVs
    ref_cnvs = pd.DataFrame() 
    query_cnvs = pd.DataFrame()
    if 'columns' in refj['CNV']['PASS'].keys():
        ref_cnvs = pd.DataFrame(columns=refj['CNV']['PASS']['columns'],data=refj['CNV']['PASS']['data'])
        if ref_cnvs.shape[0] > 0:
            ref_cnvs['variant'] = ref_cnvs.apply(lambda x: x['gene']+':'+x['type']+'(n='+str(x['copynumber'])+')',axis=1)
        else:
            ref_cnvs['variant'] = 'NA'

    if 'columns' in queryj['CNV']['PASS'].keys():
        query_cnvs = pd.DataFrame(columns=queryj['CNV']['PASS']['columns'],data=queryj['CNV']['PASS']['data'])
        if query_cnvs.shape[0] > 0:
            query_cnvs['variant'] = query_cnvs.apply(lambda x: x['gene']+':'+x['type']+'(n='+str(x['copynumber'])+')',axis=1)
        else:
            query_cnvs['variant'] = 'NA'
            
    out_cnvs = compare_variants(ref_cnvs,query_cnvs,['type','gene','copynumber'])
    out_cnvs.columns = 'cnvs_' + out_cnvs.columns
#    out_cnvs['case'] = row['case']
#    out_cnvs.set_index('case')

    # SVs/Fusions
    ref_svs = pd.DataFrame() 
    query_svs = pd.DataFrame()
    if 'columns' in refj['FUSIONS']['PASS'].keys():
        ref_svs = pd.DataFrame(columns=refj['FUSIONS']['PASS']['columns'],data=refj['FUSIONS']['PASS']['data'])
        if ref_svs.shape[0] > 0:
            ref_svs['variant'] = ref_svs.apply(lambda x: x['gene']+':'+x['psyntax']+'(n='+x['abundance']+')',axis=1)
        else:
            ref_svs['variant'] = 'NA'

    if 'columns' in queryj['FUSIONS']['PASS'].keys():
        query_svs = pd.DataFrame(columns=queryj['FUSIONS']['PASS']['columns'],data=queryj['FUSIONS']['PASS']['data'])
        if query_svs.shape[0] > 0:          
            query_svs['variant'] = query_svs.apply(lambda x: x['gene']+':'+x['psyntax']+'(n='+x['abundance']+')',axis=1)
        else:
            query_svs['variant'] = 'NA'

    out_svs = compare_variants(ref_svs,query_svs,['type','strand1','gene1','region1','strand2','gene2','region2'])
    out_svs.columns = 'svs_' + out_svs.columns
#    out_svs['case'] = row['case']
#    out_svs.set_index('case')

    out = pd.DataFrame.from_dict([{'case':row['case'],'path.ref': os.path.abspath(row['path.ref']),'path.query':os.path.abspath(row['path.query'])}])

    out = pd.concat([out,out_variants,out_cnvs,out_svs],axis=1)

    return(out.reset_index())

#
# Script
#

parser = argparse.ArgumentParser(description='Compare GatewaySeq cases')
parser.add_argument('referenceDir',help='Reference case dir')
parser.add_argument('queryDir',help='Query case dir')
parser.add_argument('-o','--out',default=None,help='Query case dir')

args = parser.parse_args()

outfile = sys.stdout
if args.out is not None:
    outfile = args.out

refDat = pd.DataFrame()
refDat['path'] = search_files(args.referenceDir, "*.report.json")
refDat['case'] = refDat.apply(lambda x: os.path.splitext(os.path.basename(x['path']))[0].replace(".report",""),axis=1)

queryDat = pd.DataFrame()
queryDat['path'] = search_files(args.queryDir, "*.report.json")
queryDat['case'] = queryDat.apply(lambda x: os.path.splitext(os.path.basename(x['path']))[0].replace(".report",""),axis=1)

dat = pd.merge(refDat,queryDat,on='case',suffixes=['.ref','.query'],how='inner')

if dat.shape[0] < 1:
    sys.exit("No reports for matching cases in reference and query directories")

compared = pd.DataFrame()
for index, row in dat.iterrows():
    compared = pd.concat([compared,compare_cases(row)],axis=0)

#if dat.shape[0] <= 1:
#    compared = compare_cases(dat.iloc[0])
#else:
#    compared = dat.apply(lambda row: compare_cases(row),axis=1,result_type='expand')

compared.iloc[:,1:].to_csv(outfile,sep="\t",index=False,na_rep='NA')

