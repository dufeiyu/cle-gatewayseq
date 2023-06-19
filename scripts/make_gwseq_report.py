#!/usr/bin/env python

import sys, os, re, tempfile, csv, pysam, json, binascii, argparse, subprocess
import sqlite3
import pandas as pd
import pyranges as pr
import numpy as np
from time import gmtime, strftime
from cyvcf2 import VCF
from pathlib import Path

def get_latest_tag(git_directory):
    command = ['git', '--git-dir', os.path.join(git_directory,".git"), 'describe', '--tags', '--abbrev=0']
    try:
        output = subprocess.check_output(command, stderr=subprocess.STDOUT).decode().strip()
        return output
    except subprocess.CalledProcessError as e:
        # Handle the case when no tags are found or an error occurs
        print(f"Error: {e.output.decode()}")
        return None

def revcomp(dna):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}  # DNA complement pairs
    reverse_complement = "".join(complement.get(base, base) for base in reversed(dna))
    return reverse_complement

def pos2codon(exonstart,exonend,cpos,pos,strand):
    if pos<=exonend and pos>=exonstart:
        if strand == '+':
            return(int((pos-exonstart + cpos) / 3 + .99))

        elif strand=='-':
             return(int((exonend-pos + cpos) / 3 + .99))

def make_ranges(lst):
    s = e = None
    r = []
    for i in sorted(lst):
        if s is None:
            s = e = i
        elif i == e or i == e + 1:
            e = i
        else:
            r.append((s, e))
            s = e = i
    if s is not None:
        r.append((s, e))

    out = ','.join(['%d' % i if i == j else '%d-%d' % (i, j) for (i, j) in r])
    
    return(out)

def decode_hex(string):
    hex_string = string.group(0).replace('%', '')
    return binascii.unhexlify(hex_string).decode('utf-8')

def convert_aa(codon):
    three = ["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "Ter"]
    one  = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*"]
    
    for i in range(0,len(three)):
        p = re.compile(three[i])
        codon = p.sub(one[i],codon)

    return codon


def check_qc_reference_ranges(value, minimum, maximum, unit):
    if minimum != '' and maximum != '':
        range_str = '>' + minimum + unit + ' ' + '<' + maximum + unit
        if value < float(minimum) or value > float(maximum):
            range_str = range_str + ' (!)'
    elif minimum != '':
        range_str = '>' + minimum + unit
        if value < float(minimum):
            range_str = range_str + ' (!)'
    elif maximum != '':
        range_str = '<' + maximum + unit
        if value > float(maximum):
            range_str = range_str + ' (!)'
    else:
        sys.exit("Reference range must have minimum and or maximum")
    return range_str

# add filter if the call is the wrong direction of what is clinically actionable
def testcnv(row,loss,gain,ratios):
    if row['Segment_Call']!='0' and ((row['Segment_Name'] in loss and row['Segment_Call']!='-') or (row['Segment_Name'] in gain and row['Segment_Call']!='+')):
        if row['Filter']=='PASS':
            return 'WrongSign'
        else:
            return ';'.join([row['Filter'],'WrongSign'])
        
    elif (row['Segment_Name'] in loss and row['Segment_Call']=='-' and row['Segment_Mean']>ratios[0]) or (row['Segment_Name'] in gain and row['Segment_Call']=='+' and row['Segment_Mean']<ratios[1]):
        if row['Filter']=='PASS':
            return 'cnvRatio'
        else:
            return ';'.join([row['Filter'],'cnvRatio'])

    else:
        return row['Filter']


def dbAnnotateVariants(db,queryPass,queryFiltered,params):                
    cursor = db.execute(queryPass, params)
    countPass = cursor.fetchone()[0]

    cursor = db.execute(queryFiltered, params)
    countFiltered = cursor.fetchone()[0]

    return 'PASS='+str(countPass)+'/'+'Filtered='+str(countFiltered)


def dbUpdateVariants(db,table,info,vardf):
     # first remove existing records for this case
    query = "DELETE FROM " + table + " WHERE name = ? AND mrn = ? AND accession = ?"
    cursor = db.cursor()
    cursor.execute(query, (info['name'],info['mrn'],info['accession']))
    db.commit()
    # now add records
    df = vardf.copy()
    for k in info.keys():
        df[k] = info[k]

    df.to_sql(table, db, if_exists='append', index=False)

#
# Script
#

parser = argparse.ArgumentParser(description='Make GatewaySeq report')
parser.add_argument('-n','--name',required=True,help='Sample name')
parser.add_argument('-d','--dir',required=True,help='Output directory')

parser.add_argument('-q','--qcrangejsonfile',required=False,help='QCReferenceRanges.json')
parser.add_argument('-m','--mrn',default='NONE',help='Sample MRN number')
parser.add_argument('-a','--accession',default='NONE',help='Sample accession number')
parser.add_argument('-s','--specimen',default='NONE',help='Sample specimen type')
parser.add_argument('-b','--DOB',default='NONE',help='Date of birth')
parser.add_argument('-e','--exception',default='NONE',help='Exception')
parser.add_argument('-i','--runinfostr',default='NONE',help='Illumina Run Information String')
parser.add_argument('--variantdb',required=False,default=None,help='Sqlite database file for variant lookup and recording')
parser.add_argument('--noupdate',required=False,default=False,action='store_true',help='Dont update variant database')
parser.add_argument('--excludevaf',default=None,type=float,help='Minimum VAF for reporting, below this will not be shown')
parser.add_argument('--minvaf',default=None,type=float,help='Minimum validated VAF')
parser.add_argument('--minreads',default=None,type=float,help='Minimum alt reads to report variant')
parser.add_argument('--maxaf',default=None,type=float,help='Maximum population allele frequency (in %) for potential somatic variants')
parser.add_argument('--maxhotspotaf',default=None,type=float,help='Maximum population allele frequency (in %) for hotspot variants')
parser.add_argument('--maxsvlen',default=None,type=float,help='Maximum SV length by breakend analysis')
parser.add_argument('--minsvreads',default=None,type=float,help='Min SR reads for SVs')
parser.add_argument('--minsvpairedreads',default=None,type=float,help='Min PR reads for SVs')
parser.add_argument('--minsvabundance',default=None,type=float,help='Min % abundance for SVs')

args = parser.parse_args()

# get repo location. 
repoLocation = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

caseinfo = {}
caseinfo['version'] = get_latest_tag(repoLocation)
caseinfo['name'] = args.name
caseinfo['mrn'] = args.mrn
caseinfo['DOB'] = args.DOB
caseinfo['accession'] = args.accession
caseinfo['specimen'] = args.specimen
caseinfo['casedir'] = args.dir
caseinfo['exception'] = args.exception
caseinfo['run_info_str'] = args.runinfostr
caseinfo['date'] = strftime("%Y-%m-%d_%H:%M:%S", gmtime())

qcdatafile = args.qcrangejsonfile
variantdb = args.variantdb

# if the case is a positive control, then do not update the variant database
if 'POSITIVE' in caseinfo['exception'].upper():
    variantdb = None

if not qcdatafile:
    qcdatafile = os.path.join(repoLocation,"accessory_files/GWSeq.QCMetrics.json")

if caseinfo['specimen'] == 'BM':
    caseinfo['specimen'] = 'Bone Marrow'
elif caseinfo['specimen'] == 'PB':
    caseinfo['specimen'] = 'Peripheral Blood'


#########################################
#
# Get thresholds and reference values for QC Metrics
#
#########################################

if not Path(qcdatafile).is_file():
    sys.exit("QC json file " + str(qcdatafile) + " not valid.")

qcranges = {}
with open(qcdatafile, 'r') as json_file:
    qcranges = json.load(json_file)

qcranges['qcdatafile'] = qcdatafile
qcranges['variantdb'] = variantdb

# update commandline thresholds/cutoffs if passed -- this overrides the ones in the QC JSON
for arg_name, arg_value in vars(args).items():
    if arg_name in qcranges.keys():
        qcranges[arg_name] = arg_value


minFractionCovered = float(qcranges['Target fraction at coverage'].split(',')[0])
minTargetCoverage = [ int(x) for x in qcranges['Minimum target coverage'].split(',') ]
minCoverage = float(qcranges['COVERAGE SUMMARY: Average alignment coverage over target region'].split(',')[0])

# get annotations classified as nonsynonymous
nonSynon = qcranges['NONSYNONYMOUS_ANNOTATIONS']
nonCodingGeneList = qcranges['NONCODINGVARIANTLIST']

# TMB quartiles
tmbquartiles = [ float(x) for x in qcranges['TMBQUARTILES'].split(',') ]

if len(tmbquartiles)!=4:
    sys.exit("TMB quartiles incorrect! must be 4")

tmbquartiles = dict(zip(tmbquartiles,['0-25','25-50','50-75','75-100']))

#########################################
#
# Get files from the case directory 
#
#########################################

vcffile = list(Path(caseinfo['casedir']).rglob('*.annotated.vcf.gz'))[0]
if not vcffile.is_file():
    sys.exit("VCF file " + str(vcffile) + " not valid.")

svvcffile = list(Path(caseinfo['casedir']).rglob('*.sv_annotated.vcf.gz'))[0]
if not svvcffile.is_file():
    sys.exit("SV VCF file " + str(svvcffile) + " not valid.")

cnvsegoutput = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.seg.called.merged'))[0]
if not cnvsegoutput.is_file():
    sys.exit("CNV file " + str(cnvsegoutput) + " not valid.")

haplotect = list(Path(caseinfo['casedir']).rglob('*.haplotect.txt'))[0]
haplotectloci = list(Path(caseinfo['casedir']).rglob('*.haplotectloci.txt'))[0]
if not haplotect.is_file() or not haplotectloci.is_file():
    sys.exit("Haplotect output " + str(haplotect) + " not valid.")

geneFusionList = os.path.join(repoLocation,"accessory_files/GWSeq.gene_fusions.txt")
if not Path(geneFusionList).is_file():
    sys.exit("Gene fusion list " + str(geneFusionList) + " not valid.")

mappingmetrics = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.mapping_metrics.csv'))[0]
targetmetrics = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.target_bed_coverage_metrics.csv'))[0]
umimetrics = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.umi_metrics.csv'))[0]

# the below command is if the targetbed is from the dragen. we changed the file so in order to not
# rerun the dragen during the validation we will just pass the original BED file for now
#covBedPr = pr.PyRanges(pd.read_csv(genetargetbed, skiprows=1, header=None, names="Chromosome Start End Gene Region ReadCov R1Cov R2Cov".split(), sep="\t"))
#genetargetbed = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.target_bed_read_cov_report.bed'))[0]
#alltargetbed = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.qc-coverage-region-1_read_cov_report.bed'))[0]

alltargetbed = Path(os.path.join(repoLocation,"accessory_files/GWSeq.all.hg38.bed"))

coveragebed = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.qc-coverage-region-1_full_res.bed'))[0]
MSIjson = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.microsat_output.json'))[0]
TMBcsv = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.tmb.metrics.csv'))[0]

if not mappingmetrics.is_file() or not targetmetrics.is_file() or not umimetrics.is_file() or not alltargetbed.is_file() or not coveragebed.is_file() or not MSIjson.is_file() or not TMBcsv.is_file():
    sys.exit("DRAGEN metrics files not found.")

#########################################
#
# Get run info
#
#########################################

if caseinfo['run_info_str'] == 'NONE':
    caseinfo['runid'] = 'NONE'
    caseinfo['instrument'] = 'NONE'
    caseinfo['spec'] = 'NONE'
    caseinfo['flowcell'] = 'NONE'
else:
    run_info = caseinfo['run_info_str'].split(',')
    caseinfo['runid'] = run_info[0]
    caseinfo['instrument'] = ' '.join((run_info[1], 'Side', run_info[2]))
    caseinfo['spec'] = 'x'.join(run_info[-4:])
    caseinfo['flowcell'] = ' '.join((run_info[3], run_info[4], caseinfo['spec']))
    

#########################################
#
# Set up dataframes for all data
#
#########################################

# this is for dragen/umi metrics
qcdf = pd.DataFrame(columns=['metric','value','qcmetric'])

# this is for exon/gene coverage metrics
covqcdf = pd.DataFrame(columns=['Gene','Type','Region','Mean','Covered1','Covered2'])

# dataframe with small variants
variants = pd.DataFrame(columns=['type','filter','chrom','pos','ref','alt','gene','transcript','consequence','csyntax','psyntax','exon','popaf','annotations','coverage','altreads','vaf','dblookup'])

# dataframe with SVs
svs = pd.DataFrame(columns=['type','psyntax','csyntax','filter','chr1','pos1','strand1','gene1','region1','chr2','pos2','strand2','gene2','region2','id','abundance','info','dblookup'])

# dataframe with CNVs
cnvs = pd.DataFrame(columns=['type','gene','chrom','start','end','copy_ratio','qual','filter','copynumber','ploidy','dblookup'])

# this is for tmb and msi
biomarkers = pd.DataFrame(columns=['metric','value'])

# df for haplotect loci
haplotectlocidf = pd.DataFrame(columns=['chr','SNP1','SNP2','all11','all12','all21','all22','popn_counts','distance','total_count','sample_counts'])

#
# Setup/open variant db
#

dbcon = None

if variantdb is not None:
     # create database, if it doesnt exist
     if not os.path.exists(variantdb) or os.path.getsize(variantdb) == 0:
        print(f"Creating variant database {variantdb}")
        dbcon = sqlite3.connect(variantdb)
        pd.concat([pd.DataFrame(columns=caseinfo.keys()),variants]).to_sql('variants', dbcon, if_exists='replace',index=False)
        pd.concat([pd.DataFrame(columns=caseinfo.keys()),svs]).to_sql('svs', dbcon, if_exists='replace',index=False)
        pd.concat([pd.DataFrame(columns=caseinfo.keys()),cnvs]).to_sql('cnvs', dbcon, if_exists='replace',index=False)
        pd.concat([pd.DataFrame(columns=caseinfo.keys()),covqcdf]).to_sql('coverage', dbcon, if_exists='replace',index=False)
        pd.concat([pd.DataFrame(columns=caseinfo.keys()),qcdf]).to_sql('qcdata', dbcon, if_exists='replace',index=False)
        pd.concat([pd.DataFrame(columns=caseinfo.keys()),biomarkers]).to_sql('biomarkers', dbcon, if_exists='replace',index=False)
        pd.concat([pd.DataFrame(columns=caseinfo.keys()),haplotectlocidf]).to_sql('haplotectloci', dbcon, if_exists='replace',index=False)
     else:
        dbcon = sqlite3.connect(variantdb)
         

#########################################
#
# Collect QC metrics
#
#########################################

print("Collecting DRAGEN qc metrics...",file=sys.stderr)

# read in mapping metrics
df = pd.read_csv(mappingmetrics,sep=',',names=['group','readgroup','metric','value','percent'])
df = df[df['group']=='MAPPING/ALIGNING SUMMARY'].drop(columns='readgroup').copy()
df['metric'] = df['group'] + ': ' + df['metric']
df = df.drop(columns='group')
dfpct = df[df['percent'].notna()].copy()
df = df.drop(columns='percent')
dfpct['metric'] = dfpct['metric'].apply(lambda x: x + ' (%)')
dfpct['value'] = dfpct['percent']
dfpct = dfpct.drop(columns='percent')

df['qcmetric'] = 0
df.loc[df['metric'].isin(qcranges.keys()),'qcmetric'] = 1
qcdf = pd.concat([qcdf,df])
dfpct['qcmetric'] = 0
dfpct.loc[dfpct['metric'].isin(qcranges.keys()),'qcmetric'] = 1
qcdf = pd.concat([qcdf,dfpct])

# read in umi metrics
df = pd.read_csv(umimetrics,sep=',',names=['group','readgroup','metric','value','percent'])
df['group'] = 'UMI SUMMARY'
df = df.drop(columns='readgroup').copy()
df['metric'] = df['group'] + ': ' + df['metric']
df = df.drop(columns='group')
dfpct = df[df['percent'].notna()].copy()
df = df.drop(columns='percent')
dfpct['metric'] = dfpct['metric'].apply(lambda x: x + ' (%)')
dfpct['value'] = dfpct['percent']
dfpct = dfpct.drop(columns='percent')
consensusReads = round(df.loc[df.metric=='UMI SUMMARY: Consensus pairs emitted','value'].astype(int).tolist()[0] * 2 / df.loc[df.metric=='UMI SUMMARY: Number of reads','value'].astype(int).tolist()[0] * 100,1)
duplicateReads = round(100-consensusReads,1)

dfpct = pd.concat([dfpct,pd.DataFrame.from_dict({0:['UMI SUMMARY: Consensus reads (%)',consensusReads],1:['UMI SUMMARY: Duplicate reads (%)',duplicateReads]},orient='index',columns=['metric','value'])],axis=0)

df['qcmetric'] = 0
df.loc[df['metric'].isin(qcranges.keys()),'qcmetric'] = 1
qcdf = pd.concat([qcdf,df])
dfpct['qcmetric'] = 0
dfpct.loc[dfpct['metric'].isin(qcranges.keys()),'qcmetric'] = 1
qcdf = pd.concat([qcdf,dfpct])
    
# read in target metrics
df = pd.read_csv(targetmetrics,sep=',',names=['group','readgroup','metric','value','percent'])
df = df.drop(columns='readgroup').copy()
df['metric'] = df['group'] + ': ' + df['metric']
df = df.drop(columns='group').copy()
dfpct = df[df['percent'].notna()].copy()
df = df.drop(columns='percent')
dfpct['metric'] = dfpct['metric'].apply(lambda x: x + ' (%)')
dfpct['value'] = dfpct['percent']
dfpct = dfpct.drop(columns='percent')

df['qcmetric'] = 0
df.loc[df['metric'].isin(qcranges.keys()),'qcmetric'] = 1
qcdf = pd.concat([qcdf,df])
dfpct['qcmetric'] = 0
dfpct.loc[dfpct['metric'].isin(qcranges.keys()),'qcmetric'] = 1
qcdf = pd.concat([qcdf,dfpct])

# get coverage file
fullResCovPr = pr.PyRanges(pd.read_csv(coveragebed, header=None, names="Chromosome Start End Coverage".split(), sep="\t"))

ntAtCovlevel1 = sum(fullResCovPr[fullResCovPr.Coverage>minTargetCoverage[0]].lengths())
ntAtCovlevel2 = sum(fullResCovPr[fullResCovPr.Coverage>minTargetCoverage[1]].lengths())

qcdf.loc[len(qcdf.index)] = ['COVERAGE SUMMARY: Target at ' + str(minTargetCoverage[0]) + 'x (%)', ntAtCovlevel1/sum(fullResCovPr.lengths())*100,1]
qcdf.loc[len(qcdf.index)] = ['COVERAGE SUMMARY: Target at ' + str(minTargetCoverage[1]) + 'x (%)', ntAtCovlevel2/sum(fullResCovPr.lengths())*100,1]

####################
#
# Haplotect
#
####################

haplotectdf = pd.read_csv(haplotect,sep='\t')
haplotectdf.columns = haplotectdf.columns.str.replace('#', '')
haplotectdf = haplotectdf.iloc[:, :-2]
haplotectdf.fillna(0, inplace=True)

haplotectlocidf = pd.read_csv(haplotectloci,sep='\t',skiprows=2)
haplotectlocidf.columns = haplotectlocidf.columns.str.replace('#', '')
haplotectlocidf = haplotectlocidf.iloc[:, :-1]

# update database w/ QC and haplotect output
if variantdb is not None:
    df = haplotectdf.transpose().reset_index().drop(index=0)
    df.columns = ['metric','value']
    df['qcmetric'] = 1
    dbUpdateVariants(dbcon,'qcdata',caseinfo,pd.concat([qcdf,df]))
    if args.noupdate is False:
        dbUpdateVariants(dbcon,'haplotectloci',caseinfo,haplotectlocidf)


# get gene coverage info

# intersect full res coverage bedfile w/ coverage QC bed file to calculate coverage

# the below command is if the targetbed is from the dragen. we changed the file so in order to not
# rerun the dragen during the validation we will just pass the original BED file for now
#covBedPr = pr.PyRanges(pd.read_csv(alltargetbed, skiprows=1, header=None, names="Chromosome Start End Gene Region ReadCov R1Cov R2Cov".split(), sep="\t"))

covBedPr = pr.PyRanges(pd.read_csv(alltargetbed, skiprows=1, header=None, names="Chromosome Start End Gene Region".split(), sep="\t"))

# note, this file has all the transcript IDs, so get those now
transcripts = covBedPr.df.Region.str.split('|',expand=True).replace('\s\+\s\S+','',regex=True).drop_duplicates().loc[:,2:]
transcripts.columns = ['id','cdsStart','cdsEnd','strand']
transcripts = transcripts[transcripts.id!='.'].copy()
transcripts['cdsStart'] = transcripts['cdsStart'].astype(int)
transcripts['cdsEnd'] = transcripts['cdsEnd'].astype(int)

df = covBedPr.join(fullResCovPr).df
df['nt'] = df[['End','End_b']].min(axis=1) - df[['Start','Start_b']].max(axis=1)
df['tcov'] = df['nt'] * df['Coverage']

# mean gene coverage and fraction of gene targets at minTargetCov or higher 
genecovdf = df[df.Region.str.contains(r'GOAL_genes',na=False)].groupby('Gene')[['nt','tcov']].sum().reset_index()
genecovdf['Type'] = 'Gene'
genecovdf['Region'] = genecovdf['Gene']
genecovdf['Mean'] = genecovdf.apply(lambda v: v['tcov']/v['nt'],axis=1)
genecovdf = pd.merge(genecovdf,df[df.Region.str.contains(r'GOAL_genes',na=False) & (df.Coverage>=minTargetCoverage[0])].groupby('Gene')[['nt']].sum().reset_index().rename(columns={'nt':'Covered1','Gene':'Region'}),on='Region')
genecovdf = pd.merge(genecovdf,df[df.Region.str.contains(r'GOAL_genes',na=False) & (df.Coverage>=minTargetCoverage[1])].groupby('Gene')[['nt']].sum().reset_index().rename(columns={'nt':'Covered2','Gene':'Region'}),on='Region')
genecovdf['Covered1'] = genecovdf.Covered1/genecovdf.nt*100
genecovdf['Covered2'] = genecovdf.Covered2/genecovdf.nt*100

df = covBedPr.join(fullResCovPr).df

svgenes = df[df.Region.str.contains(r'GOAL_rearrangements',na=False)]['Gene'].drop_duplicates().tolist()

df['nt'] = df[['End','End_b']].min(axis=1) - df[['Start','Start_b']].max(axis=1)
df['tcov'] = df['nt'] * df['Coverage']

# mean gene coverage and fraction of gene targets at minTargetCov or higher 
df2 = df[df.Region.str.contains(r'GOAL_rearrangements',na=False)].copy()
df2['Gene'] = df2['Gene'].str.extract(r'(\S+_Fusion)')
svcovdf = df2.groupby('Gene')[['nt','tcov']].sum().reset_index()
svcovdf['Type'] = 'Fusion'
svcovdf['Region'] = svcovdf['Gene']
svcovdf['Mean'] = svcovdf.apply(lambda v: v['tcov']/v['nt'],axis=1)
svcovdf = pd.merge(svcovdf,df2[df2.Region.str.contains(r'GOAL_rearrangements',na=False) & (df2.Coverage>=minTargetCoverage[0])].groupby('Gene')[['nt']].sum().reset_index().rename(columns={'nt':'Covered1'}),on='Gene')
svcovdf = pd.merge(svcovdf,df2[df2.Region.str.contains(r'GOAL_rearrangements',na=False) & (df2.Coverage>=minTargetCoverage[1])].groupby('Gene')[['nt']].sum().reset_index().rename(columns={'nt':'Covered2'}),on='Gene')
svcovdf['Covered1'] = svcovdf.Covered1/svcovdf.nt*100
svcovdf['Covered2'] = svcovdf.Covered2/svcovdf.nt*100

# mean MSI loci coverage and fraction of targets at minTargetCov or higher 
othercovdf = df[df.Region.str.contains(r'_MSI|haplotect',na=False)].copy()
othercovdf['Gene'] = 'MSI'
othercovdf.loc[othercovdf.Region.str.contains(r'haplotect',na=False),'Gene'] = 'HAPLOTECT'
x = pd.merge(othercovdf.groupby(['Gene'])[['nt','tcov']].sum().reset_index(),othercovdf[othercovdf.Coverage>=minTargetCoverage[0]].groupby(['Gene'])[['nt']].sum().reset_index().rename(columns={'nt':'Covered1'}),on=['Gene'])
othercovdf = pd.merge(x,othercovdf[othercovdf.Coverage>=minTargetCoverage[1]].groupby(['Gene'])[['nt']].sum().reset_index().rename(columns={'nt':'Covered2'}),on=['Gene'])
othercovdf['Mean'] = othercovdf.tcov/othercovdf.nt
othercovdf['Covered1'] = othercovdf.Covered1/othercovdf.nt*100
othercovdf['Covered2'] = othercovdf.Covered2/othercovdf.nt*100
othercovdf['Region'] = othercovdf['Gene']
othercovdf['Type'] = 'MSI/HAPLOTECT'

covqcdf = pd.concat([genecovdf[['Gene','Type','Region','Mean','Covered1','Covered2']],svcovdf[['Gene','Type','Region','Mean','Covered1','Covered2']],othercovdf[['Gene','Type','Region','Mean','Covered1','Covered2']]])

# update database with new results, if necessary
if variantdb is not None and args.noupdate is False:
    dbUpdateVariants(dbcon,'coverage',caseinfo,covqcdf)


#########################################
#
# Get small variants
#
#########################################

print("Gathering gene variants...",file=sys.stderr)

vcf = VCF(vcffile)

# get VEP fields
vep = {}
i = 0
for j in vcf.get_header_type('CSQ')['Description'].strip('"').split("|"):
    vep[j] = i
    i+=1

# get variants
for variant in vcf:

    vartype = ''
    if len(variant.REF) == len(variant.ALT[0]):
        vartype = 'SNV'
    else:
        vartype = 'INDEL'

    varfilter = 'PASS'
    if variant.FILTER is not None:
        varfilter = variant.FILTER
 
    abundance = round(variant.format('AF')[0][0] * 100,2)

    # dont even show variants below this VAF (%)
    if abundance < qcranges['ExcludeVAF']:
        continue

    if abundance < qcranges['MinVAF'] or variant.format('AD')[0][1] < qcranges['MinReads'] and varfilter=='PASS':
        varfilter = 'LowReads'
        
    gt = [variant.REF] + variant.ALT
    genotype = gt[variant.genotypes[0][0]] + '/' + gt[variant.genotypes[0][1]]
            
    # get VEP annotation
    csq = variant.INFO['CSQ']
    
    if csq is None:
        sys.exit("No VEP fields")
    
    gene='NA'
    transcript='NA'
    csyntax='NA'
    psyntax='NA'
    consequence='NA'
    exon='NA'
    intron='NA'
    popmaf = 'NA'
    customannotation = 'NA'
    for i in variant.INFO['CSQ'].split(','):
        csq = i.split("|")

       # get pop allele frequency. This is present for each transcript annotation, but is always the same 
        if csq[vep['GWSEQ_BLACKLIST']] != '':
            varfilter = 'BLACKLIST'

        # get pop allele frequency. This is present for each transcript annotation, but is always the same 
        if csq[vep['MAX_AF']] != '':
            popmaf = float(csq[vep['MAX_AF']])

        # check if this is in the list of transcripts. only variants annotated with a known transcript will be reported
        if transcripts.id.str.contains(csq[vep['Feature']]).any():
            transcript = csq[vep['Feature']]
            gene = csq[vep['SYMBOL']]
            consequence = csq[vep['Consequence']].split("&")[0]

            csyntax = csq[vep['HGVSc']].split(":")
            if len(csyntax) > 1:
                csyntax = csyntax[1]
            else:
                if csq[vep['STRAND']]!=transcripts[transcripts.id==transcript]['strand'].tolist()[0]:
                    sys.exit('strands dont match')

                if 'upstream' in consequence:
                    if csq[vep['STRAND']]==1:
                        # *  ---->                        
                        distance = transcripts[transcripts.id==transcript]['cdsStart'].min() - variant.POS
                        csyntax = "c.-"+str(distance)+variant.REF+'>'+csq[0]
                    else:
                        # <---- *
                        distance = variant.POS - transcripts[transcripts.id==transcript]['cdsEnd'].max()
                        csyntax = "c.-"+str(distance)+revcomp(variant.REF)+'>'+revcomp(csq[0])

                elif 'downstream' in consequence:
                    if csq[vep['STRAND']]==1:
                        # ---->  *                        
                        distance = variant.POS - transcripts[transcripts.id==transcript]['cdsEnd'].max()
                        csyntax = "c.+"+str(csq[vep['DISTANCE']])+variant.REF+'>'+csq[0]
                    else:
                        # *  <----
                        distance = transcripts[transcripts.id==transcript]['cdsEnd'].min() - variant.POS
                        csyntax = "c.+"+str(csq[vep['DISTANCE']])+revcomp(variant.REF)+'>'+revcomp(csq[0])
                else:
                    csyntax = 'noncoding'

            psyntax = csq[vep['HGVSp']].split(":")
            if len(psyntax) > 1:
                psyntax = convert_aa(psyntax[1])
                psyntax = re.sub("\%3D","=",psyntax)
            else:
                psyntax = csyntax
        
            impact = csq[vep['IMPACT']]
            exon = csq[vep['EXON']] or 'NA'
            intron = csq[vep['INTRON']] or 'NA'
            customannotation = csq[vep['Existing_variation']] or 'NA'            

    # convert pop maf to percent
    if popmaf!='NA':
        popmaf = float(popmaf)*100

    # skip all variants >1% PAF
    if popmaf!='NA' and popmaf >= qcranges['MaxHotspotAF']:
        continue

    # skip all variants > 0.1% and not in a hotspot 
    if popmaf!='NA' and popmaf >= float(qcranges['MaxAF']) and variant.INFO.get('hotspot') is None:
        continue

    # only include all variants <=0.1% and ns or specific noncoding variants 
    if consequence in nonSynon or (gene in nonCodingGeneList.keys() and consequence in nonCodingGeneList[gene]):        
        variants = pd.concat([variants,pd.DataFrame([dict(zip(variants.columns,[vartype,varfilter,str(variant.CHROM),str(variant.POS),variant.REF,variant.ALT[0],gene,transcript,consequence,csyntax,psyntax,exon,str(popmaf) + '%',customannotation,str(variant.format("DP")[0][0]),str(variant.format("AD")[0][1]),str(abundance)+"%",'NA']))])])

# now query database to get counts of variants like this one.
if variantdb is not None and variants.shape[0] > 0:
    queryPass = '''SELECT COUNT(*) FROM variants WHERE filter = "PASS" AND name != ? AND chrom = ? AND pos = ? AND ref = ? AND alt = ?'''
    queryFiltered = '''SELECT COUNT(*) FROM variants WHERE filter != "PASS" AND name != ? AND chrom = ? AND pos = ? AND ref = ? AND alt = ?'''
    variants['dblookup'] = variants.apply(lambda x: dbAnnotateVariants(dbcon,queryPass,queryFiltered,(caseinfo['name'],x['chrom'],x['pos'],x['ref'],x['alt'])),axis=1)

    if args.noupdate is False:
        dbUpdateVariants(dbcon,'variants',caseinfo,variants)

########################
#
# Get CNVs
#
########################

cnvratios = [ float(x) for x in qcranges['CNVRATIOS'].split(',') ]

cnv = pd.read_csv(cnvsegoutput,delimiter='\t')
cnvcalls = pd.concat([cnv[cnv['Segment_Name'].isin(qcranges['CNVTARGETS']['loss'])],cnv[cnv['Segment_Name'].isin(qcranges['CNVTARGETS']['gain'])]])

cnvcalls['Filter'] = cnvcalls.apply(lambda r: testcnv(r,qcranges['CNVTARGETS']['loss'],qcranges['CNVTARGETS']['gain'],cnvratios), axis=1)
cnvcalls['type'] = ''
cnvcalls.loc[cnvcalls['Segment_Call']=='-','type'] = 'LOSS'
cnvcalls.loc[cnvcalls['Segment_Call']=='+','type'] = 'GAIN'
cnvcalls.loc[cnvcalls['Segment_Call']=='0','type'] = 'REF'
cnvcalls = cnvcalls[['type','Filter','Chromosome','Start','End','Segment_Name','Segment_Mean','Qual','Copy_Number','Ploidy']].copy()
cnvcalls.rename(columns={'Segment_Name':'gene','Chromosome':'chrom','Start':'start','End':'end','Segment_Mean':'copy_ratio','Qual':'qual','Filter':'filter','Copy_Number':'copynumber','Ploidy':'ploidy'},inplace=True)
cnvcalls['dblookup'] = 'NA'

cnvs = pd.concat([cnvs,cnvcalls])

# now query database to get counts of variants like this one and then update the db with these variants.
if variantdb is not None and cnvs.shape[0] > 0:
    queryPass = '''SELECT COUNT(*) FROM cnvs WHERE filter = "PASS" AND name != ? AND gene = ? AND type = ?'''
    queryFiltered = '''SELECT COUNT(*) FROM cnvs WHERE filter != "PASS" AND name != ? AND gene = ? AND type = ?'''
    cnvs['dblookup'] = cnvs.apply(lambda x: dbAnnotateVariants(dbcon,queryPass,queryFiltered,(caseinfo['name'],x['gene'],x['type'])),axis=1)

    if args.noupdate is False:
        dbUpdateVariants(dbcon,'cnvs',caseinfo,cnvs)

########################
#
# Get SVs
# 
########################

svvcf = VCF(svvcffile)

# get known fusion events
knownsvs = pd.read_csv(geneFusionList,sep='\t', header=None)
knownsvs.columns = ['gene1','gene2']

# get VEP fields
vep = {}
i = 0
for j in svvcf.get_header_type('CSQ')['Description'].strip('"').split("|"):
    vep[j] = i
    i+=1

passedvars = {}
alreadydone = set()

for variant in svvcf:
    # if its a del, dup, or ins then process here. (BNDs need to capture both ends and then process)
    if variant.INFO['SVTYPE']!='BND':        

        filter = 'PASS'
        if variant.FILTER is not None:
            filter = variant.FILTER

        # get gene 1 info
        chr1 = str(variant.CHROM)
        chr2 = str(variant.CHROM)
        pos1 = str(variant.POS)
        pos2 = str(variant.INFO['END'])
        strand1 = 'NA'
        strand2 = 'NA'
        region1 = 'NA'
        region2 = 'NA'
        gene='NA'
        transcript='NA'
        csyntax = ''
        psyntax = ''

        for i in variant.INFO['CSQ'].split(','):
            csq = i.split("|")
            if transcripts.id.str.contains(csq[vep['Feature']]).any():
                transcript = csq[vep['Feature']]
                gene = csq[vep['SYMBOL']]
                consequence = csq[vep['Consequence']].split("&")[0]
                csyntax = csq[vep['HGVSc']]
                psyntax = csq[vep['HGVSp']]

                if csq[vep['STRAND']]==-1:
                   strand1 = strand2 = '-'
                else:
                    strand1 = strand2 = '+'

                if csq[vep['EXON']]!='':
                   exons = csq[vep['EXON']].replace('/','-').split('-')
                   if exons[0]==1 and exons[1]==exons[2]:
                       region1 = region2 = 'complete'
                   elif len(exons)==1:
                    region1 = 'exon' + str(exons[0])
                    region2 = 'exon' + str(exons[0])
                   else:
                    region1 = 'exon' + str(exons[0])
                    region2 = 'exon' + str(exons[1])

                elif csq[vep['INTRON']]!='':
                   introns = csq[vep['INTRON']].replace('/','-').split('-')
                   if introns[0]==1 and introns[1]==introns[2]:
                       region1 = region2 = 'complete'
                   elif len(introns)==1:
                    region1 = 'intron' + str(introns[0])
                    region2 = 'intron' + str(introns[0])
                   else:
                    region1 = 'intron' + str(introns[0])
                    region2 = 'intron' + str(introns[1])

        abundance = 0.0
        pr = (0,0)
        sr = (0,0)            
        if variant.format("SR") is not None:
            sr = variant.format("SR")[0]
            
        if variant.format("PR")[0] is not None:                
            pr =  variant.format("PR")[0]

        abundance = (sr[1] + pr[1]) / (pr[0] + pr[1] + sr[0] + sr[1]) * 100

        if sr[1] <= qcranges['MinSvReads'] or abundance < qcranges['MinSvAbundance']:
            continue

        if csyntax=='':
            csyntax = chr1 + ":g." + str(pos1) + "_" + str(pos2) + str(variant.INFO['SVTYPE']).lower()
            psyntax = str(variant.INFO['SVTYPE']).lower() + '(' + gene + ')'
            if region1!='complete':
                psyntax = psyntax + '(' + region1 + '-' + region2 + ')'

        infostring = 'PR_READS=' + str(pr[1]) + '/' + str(pr[0]+pr[1]) + ';SR_READS=' + str(sr[1]) + '/' + str(sr[0]+sr[1]) + ';CONTIG=' + str(variant.INFO.get('CONTIG'))

        # only report functional variants.
        if gene!='NA' and region1!='NA' and region2!='NA' and variant.INFO['SVTYPE'] in qcranges['SVTARGETS'].keys() and gene in qcranges['SVTARGETS'][variant.INFO['SVTYPE']]:
            svs.loc[len(svs.index)] = [variant.INFO['SVTYPE'],psyntax,csyntax,filter,chr1,str(pos1),strand1,gene,region1,chr2,str(pos2),strand2,gene,region2,str(variant.ID),str(round(abundance,1))+"%",infostring,'NA']

    else:
        passedvars[variant.ID] = variant

for v in passedvars.items():

    variant = v[1]
    mate = ''

    filter = 'PASS'
    if variant.FILTER is not None:
        filter = variant.FILTER

    # get gene 1 info
    gene1='NA'
    transcript1='NA'
    region1 = 'NA'
    strand1 = '+'
    for i in variant.INFO['CSQ'].split(','):
        csq = i.split("|")
        # if this is the list of transcripts to use for annotation or if its not and its the 'picked' one'
        if transcripts.id.str.contains(csq[vep['Feature']]).any() or knownsvs.gene1.isin([csq[vep['SYMBOL']]]).any() or knownsvs.gene2.isin([csq[vep['SYMBOL']]]).any() or (gene1=='NA' and csq[vep['PICK']]=='1'):
            transcript1 = csq[vep['Feature']]
            gene1 = csq[vep['SYMBOL']]
            if csq[vep['INTRON']]!='':
                region1 = 'intron' + str(csq[vep['INTRON']].split('/')[0])
            elif csq[vep['EXON']]!='':
                region1 = 'exon' + str(csq[vep['EXON']].split('/')[0])
            else:
                region1 = 'INTERGENIC'

            if csq[vep['STRAND']]==-1:
                strand1 = '-'


    chr1 = str(variant.CHROM)
    pos1 = variant.POS
                        
    # skip if this is the mate
    if variant.INFO.get('MATEID') in alreadydone or variant.INFO.get('MATEID') not in passedvars:
        continue
            
    # get the mate
    mate = passedvars[variant.INFO.get('MATEID')]

    matecsq = mate.INFO.get('CSQ')

    if matecsq is None:
        sys.exit("No VEP fields")
                
    gene2='NA'
    transcript2='NA'
    region2 = 'NA'
    strand2 = '+'
    for i in mate.INFO['CSQ'].split(','):
        csq = i.split("|")
        # if this is the list of transcripts to use for annotation or if its not and its the 'picked' one'
        if transcripts.id.str.contains(csq[vep['Feature']]).any() or knownsvs.gene1.isin([csq[vep['SYMBOL']]]).any() or knownsvs.gene2.isin([csq[vep['SYMBOL']]]).any() or (gene2=='NA' and csq[vep['PICK']]=='1'):
            transcript2 = csq[vep['Feature']]
            gene2 = csq[vep['SYMBOL']]
            if csq[vep['INTRON']]!='':
                region2 = 'intron' + str(csq[vep['INTRON']].split('/')[0])
            elif csq[vep['EXON']]!='':
                region2 = 'exon' + str(csq[vep['EXON']].split('/')[0])
            else:
                region2 = 'INTERGENIC'

            if csq[vep['STRAND']]==-1:
                strand2 = '-'

    chr2 = mate.CHROM
    pos2 = mate.POS
        
    # one gene must be a targeted one
    if gene1==gene2 or gene1=='NA' or gene2=='NA' or gene1=='' or gene2=='' or (not knownsvs.gene1.isin([gene1,gene2]).any() and not knownsvs.gene2.isin([gene1,gene2]).any()):
        continue

    # check to see if the gene pairs are in our 'database' otherwise list as novel SV involving one of the targeted genes
    if len(knownsvs[(knownsvs['gene1']==gene1) & (knownsvs['gene2']==gene2)])==0 and len(knownsvs[(knownsvs['gene1']==gene2) & (knownsvs['gene2']==gene1)])==0:
        filter = 'NovelSV'

    # abundance
    abundance = 0.0
    pr = (0,0)
    sr = (0,0)            
    if variant.format("SR") is not None:
        sr = variant.format("SR")[0]
            
    if variant.format("PR")[0] is not None:                
        pr =  variant.format("PR")[0]

    abundance = (sr[1] + pr[1]) / (pr[0] + pr[1] + sr[0] + sr[1]) * 100

    if sr[1] < qcranges['MinSvReads'] or abundance < qcranges['MinSvAbundance']:
        continue

    if pr[1] < qcranges['MinSvPairedReads']:
        if filter=='PASS':
            filter=='NoPairedReads'
        else:
            filter = filter + ';' + 'NoPairedReads'


    alt = variant.ALT[0]
    strand = '+'
    if alt.find("[") == 0 or alt.find("]") > 0:
        strand = '-'
            
    csyntax = ''
    if (chr1.find('X') == -1 and chr2.find('X') == -1 and chr1.find('Y') == -1 and chr2.find('Y') == -1 and int(chr1.replace('chr','')) < int(chr2.replace('chr',''))) or chr1.find('X') > -1 or chr1.find('Y') > -1: # this isnt working. Want to list lower chromosome first in these strings. If X is involved, then X first.
        csyntax = chr1 + ":g." + str(pos1) + "(+)::" + chr2 + ":g." + str(pos2) + "(" + strand + ")"
    else:
        csyntax = chr2 + ":g." + str(pos2) + "(+)::" + chr1 + ":g." + str(pos1) + "(" + strand + ")"
        
    psyntax = ''
    if (strand1=="+" and (alt.find("[")>0 or alt.find("]")>0)) or (strand1=="-" and (alt.find("]")==0 or alt.find("[")==0)):
        psyntax = gene1 + '::' + gene2
    else:
        psyntax = gene2 + '::' + gene1

    infostring = 'PR_READS=' + str(pr[1]) + '/' + str(pr[0]+pr[1]) + ';SR_READS=' + str(sr[1]) + '/' + str(sr[0]+sr[1]) + ';CONTIG=' + str(variant.INFO.get('CONTIG'))
    svs.loc[len(svs.index)] = ['BND',psyntax,csyntax,filter,chr1,str(pos1),strand1,gene1,region1,chr2,str(pos2),strand2,gene2,region2,str(variant.ID) + ";" + str(mate.ID),str(round(abundance,1))+"%",infostring,'NA']

    alreadydone.add(variant.ID)


 # now query database to get counts of variants like this one.
if variantdb is not None and svs.shape[0] > 0:
    queryPass = '''SELECT COUNT(*) FROM svs WHERE filter = "PASS" AND name != ? AND gene1 = ? AND region1 = ? AND gene2 = ? AND region2 = ?'''
    queryFiltered = '''SELECT COUNT(*) FROM svs WHERE filter != "PASS" AND name != ? AND gene1 = ? AND region1 = ? AND gene2 = ? AND region2 = ?'''
    svs['dblookup'] = svs.apply(lambda x: dbAnnotateVariants(dbcon,queryPass,queryFiltered,(caseinfo['name'],x['gene1'],x['region1'],x['gene2'],x['region2'])),axis=1)

    if args.noupdate is False:
        dbUpdateVariants(dbcon,'svs',caseinfo,svs)

########################
#
# MSI
#
########################

with open(MSIjson) as msi_json:
    MSI_json = json.load(msi_json)
msi_json.close()
msi_keys = ['TotalMicrosatelliteSitesAssessed', 'TotalMicrosatelliteSitesUnstable', 'PercentageUnstableSites', 'ResultIsValid']
msi = pd.DataFrame(list(MSI_json.items()), columns=['metric', 'value'])
msi = msi[msi.metric.isin(msi_keys)]

########################
#
# TMB
#
########################

tmb_keys = ['Total Input Variant Count in TMB region', 'Filtered Variant Count', 'Filtered Nonsyn Variant Count','Eligible Region (MB)', 'TMB', 'Nonsyn TMB']
tmb = pd.read_csv(TMBcsv,sep=',',names=['group','fill','metric','value'])[['metric','value']]
tmb = tmb[tmb.metric.isin(tmb_keys)]
tmbQuart = tmbquartiles[min((x for x in tmbquartiles.keys() if x <= float(tmb.set_index('metric')['value'].to_dict()['Nonsyn TMB'])), key=lambda x: abs(x - float(tmb.set_index('metric')['value'].to_dict()['Nonsyn TMB'])))]
tmb = pd.concat([tmb,pd.DataFrame([{'metric':'Nonsyn TMB Quartile','value':tmbQuart}])])

#
# update database with new results, if necessary
#
if variantdb is not None and args.noupdate is False:
    dbUpdateVariants(dbcon,'biomarkers',caseinfo,pd.concat([msi,tmb]))

# close database
if variantdb is not None:
    dbcon.close()

########################
#
# Start report
#
########################

print("Starting report...",file=sys.stderr)

# make dict for report and redirect output for text report
jsonout = {'CASEINFO':{},'VARIANTS':{},'CNV':{},'FUSIONS':{},'QC':{},'MSI':{},'TMB':{}}

f = open(caseinfo['name'] + ".report.txt", "w")
sys.stdout = f

print("GatewaySeq Report for " + caseinfo['name'] + " ---- Generated on: " + caseinfo['date'] + "\n")

print("*** GATEWAY CASE INFORMATION ***\n")
print("MRN:\t" + caseinfo['mrn'])
print("ACCESSION:\t" + caseinfo['accession'])
print("SPECIMEN TYPE:\t" + caseinfo['specimen'])
print("DOB:\t" + caseinfo['DOB'])
print("RUNID:\t" + caseinfo['runid'])
print("INSTRUMENT:\t" + caseinfo['instrument'])
print("FLOWCELL:\t" + caseinfo['flowcell'])
if (caseinfo['exception'] != 'NONE'):
    caseinfo['exception'] = caseinfo['exception'] + "\t(!)"
print("EXCEPTIONS:\t" + caseinfo['exception'])

print()

jsonout['CASEINFO'] = caseinfo

print("*** GENE MUTATIONS ***\n")

if variants[variants['filter']=='PASS'].shape[0] > 0:
    print(variants[variants['filter']=='PASS'].to_csv(sep='\t',header=True, index=False))
else:
    print("None Detected\n")

jsonout['VARIANTS']['PASS'] = variants[variants['filter']=='PASS'].to_dict('split')
del jsonout['VARIANTS']['PASS']['index']

print("*** FILTERED GENE MUTATIONS ***\n")

if variants[variants['filter']!='PASS'].shape[0] > 0:
    print(variants[variants['filter']!='PASS'].to_csv(sep='\t',header=True, index=False))
else:
    print("None Detected\n")

jsonout['VARIANTS']['Filtered'] = variants[variants['filter']!='PASS'].to_dict('split')
del jsonout['VARIANTS']['Filtered']['index']

print("*** COPY NUMBER VARIANTS ***\n")

if 'low input' in caseinfo['exception'].lower():
    print("Not performed\n")
    jsonout['CNV'] = False

else:
    jsonout['CNV']['PASS'] = cnvs[cnvs['filter']=='PASS'].to_dict('split')
    del jsonout['CNV']['PASS']['index']

    jsonout['CNV']['Filtered'] = cnvs[cnvs['filter']!='PASS'].to_dict('split')
    del jsonout['CNV']['Filtered']['index']

    if cnvs[cnvs['filter']=='PASS'].shape[0] > 0:
        print(cnvs[cnvs['filter']=='PASS'].to_csv(sep='\t',header=True, index=False))
    else:
        print("None Detected\n")


    print("*** FILTERED COPY NUMBER VARIANTS ***\n")

    if cnvs[cnvs['filter']!='PASS'].shape[0] > 0:
        print(cnvs[cnvs['filter']!='PASS'].to_csv(sep='\t',header=True, index=False))
    else:
        print("None Detected\n")


print("*** GENE FUSION VARIANTS ***\n")

if 'low input' in caseinfo['exception'].lower():
    print("Not performed\n")
    jsonout['FUSIONS'] = False

else:
    jsonout['FUSIONS']['PASS'] = svs[svs['filter']=='PASS'].to_dict('split')
    del jsonout['FUSIONS']['PASS']['index']
    jsonout['FUSIONS']['Filtered'] = svs[svs['filter']!='PASS'].to_dict('split')
    del jsonout['FUSIONS']['Filtered']['index']

    if svs[svs['filter']=='PASS'].shape[0] > 0:
        print(svs[svs['filter']=='PASS'].to_csv(sep='\t',header=True, index=False))    
    else:
        print("None Detected\n")

    print("*** FILTERED GENE FUSION VARIANTS ***\n")

    if svs[svs['filter']!='PASS'].shape[0] > 0:
        print(svs[svs['filter']!='PASS'].to_csv(sep='\t',header=True, index=False))
    else:
        print("None Detected\n")

print("*** MSI ***\n")

if 'low input' in caseinfo['exception'].lower():
    print("Not performed\n")
    jsonout['MSI'] = False
else:
    for msi_key, msi_value in msi.set_index('metric')['value'].to_dict().items():
        print(msi_key + ': ' + msi_value)
        jsonout['MSI'][msi_key] = msi_value

print("\n*** TMB ***\n")

if 'low input' in caseinfo['exception'].lower():
    print("Not performed\n")
    jsonout['TMB'] = False
else:
    for tmb_key, tmb_value in tmb.set_index('metric')['value'].to_dict().items():
        print(f'{tmb_key}: {tmb_value}')
        jsonout['TMB'][tmb_key] = tmb_value

print("\n*** SEQUENCING QC ***\n")

for qc in ['MAPPING/ALIGNING SUMMARY','COVERAGE SUMMARY','UMI SUMMARY']:
    jsonout['QC'][qc] = {}
    qcgroup = {k: v for k, v in qcranges.items() if k.startswith(qc)}
    for m in qcgroup.keys():
        val = qcdf.loc[qcdf['metric']==m,'value'].iat[0]
        val_ranges = qcranges[m].split(',')
        check_str = check_qc_reference_ranges(float(val), val_ranges[0], val_ranges[1], '') 
        print("\t".join((m, str(val), check_str)))
        
        jsonout['QC'][qc][m] = str(val)

print()

print("*** GENE COVERAGE QC ***\n")

xdf = covqcdf[(covqcdf.Type == "Gene")][['Gene','Mean','Covered1','Covered2']]
xdf['QC'] = np.where((xdf['Mean'] < minCoverage) | (xdf['Covered1']<minFractionCovered), '(!)', '')
xdf.rename(columns={'Covered1':'%CoveredAt'+str(int(minTargetCoverage[0]))+'x','Covered2':'%CoveredAt'+str(int(minTargetCoverage[1]))+'x'},inplace=True)
print(xdf.to_csv(sep='\t',header=True, index=False,float_format='%.1f'))

jsonout['QC']['GENE COVERAGE QC'] = xdf.to_dict('split')
jsonout['QC']['GENE COVERAGE QC'].pop('index', None)

print("*** FUSION COVERAGE QC ***\n")

xdf = covqcdf[(covqcdf.Type == "Fusion")][['Gene','Mean','Covered1','Covered2']]
xdf['QC'] = np.where((xdf['Mean'] < minCoverage) | (xdf['Covered1']<minFractionCovered), '(!)', '')
xdf.rename(columns={'Covered1':'%CoveredAt'+str(int(minTargetCoverage[0]))+'x','Covered2':'%CoveredAt'+str(int(minTargetCoverage[1]))+'x'},inplace=True)
print(xdf.to_csv(sep='\t',header=True, index=False,float_format='%.1f'))

jsonout['QC']['FUSION COVERAGE QC'] = xdf.to_dict('split')
jsonout['QC']['FUSION COVERAGE QC'].pop('index', None)

print("*** MSI/HAPLOTECT COVERAGE QC ***\n")

xdf = covqcdf[(covqcdf.Type == "MSI/HAPLOTECT")][['Gene','Mean','Covered1','Covered2']]
xdf['QC'] = np.where((xdf['Mean'] < minCoverage) | (xdf['Covered1']<minFractionCovered), '(!)', '')
xdf.rename(columns={'Covered1':'%CoveredAt'+str(int(minTargetCoverage[0]))+'x','Covered2':'%CoveredAt'+str(int(minTargetCoverage[1]))+'x'},inplace=True)
print(xdf.to_csv(sep='\t',header=True, index=False,float_format='%.1f'))

jsonout['QC']['MSI/HAPLOTECT COVERAGE QC'] = xdf.to_dict('split')
jsonout['QC']['MSI/HAPLOTECT COVERAGE QC'].pop('index', None)

print("*** Haplotect Contamination Estimate ***\n")

print(haplotectdf.iloc[:,1:].to_csv(sep='\t',header=True, index=False))
print(haplotectlocidf.to_csv(sep='\t',header=True, index=False))

jsonout['QC']['HAPLOTECT SUMMARY'] = haplotectdf.iloc[:,1:].to_dict('split')
jsonout['QC']['HAPLOTECT SUMMARY'].pop('index', None)
jsonout['QC']['HAPLOTECT LOCI'] = haplotectlocidf.iloc[:,1:].to_dict('split')
jsonout['QC']['HAPLOTECT LOCI'].pop('index', None)

print("*** GatewaySeq Assay Version " + str(caseinfo["version"]) + " ***\n")

print(qcranges["DISCLAIMER"])

jsonout['QC']['QCINFO'] = qcranges

original_stdout = sys.stdout
f.close()

# dump json
j = open(caseinfo['name'] + ".report.json", "w")
json.dump(jsonout,j,indent=" ")
j.close()
