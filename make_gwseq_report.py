#!/usr/bin/env python

import sys, os, re, tempfile, csv, pysam, json, binascii, argparse
import pandas as pd
import pyranges as pr
import numpy as np
from time import gmtime, strftime
from cyvcf2 import VCF
from pathlib import Path

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

#
# Script
#

parser = argparse.ArgumentParser(description='Make GatewaySeq report')
parser.add_argument('-n','--name',required=True,help='Sample name')
parser.add_argument('-d','--dir',required=True,help='Output directory')
parser.add_argument('-c','--coverageqc',required=True,help='Coverage QC file')
parser.add_argument('-q','--qcrangejsonfile',required=True,help='QCReferenceRanges.json')
parser.add_argument('-m','--mrn',default='NONE',help='Sample MRN number')
parser.add_argument('-a','--accession',default='NONE',help='Sample accession number')
parser.add_argument('-s','--specimen',default='NONE',help='Sample specimen type')
parser.add_argument('-b','--DOB',default='NONE',help='Date of birth')
parser.add_argument('-e','--exception',default='NONE',help='Exception')
parser.add_argument('-f','--minvaf',default=2.0,help='Minimum VAF for discovery analysis')
parser.add_argument('-i','--runinfostr',default='NONE',help='Illumina Run Information String')
parser.add_argument('-p','--maxaf',default=0.02,help='Maximum population allele frequency for potential somatic variants')

args = parser.parse_args()

caseinfo = {}
caseinfo['name'] = args.name
caseinfo['mrn'] = args.mrn
caseinfo['DOB'] = args.DOB
caseinfo['accession'] = args.accession
caseinfo['specimen'] = args.specimen
caseinfo['casedir'] = args.dir
caseinfo['covqcbedfile'] = args.coverageqc
caseinfo['maxaf'] = args.maxaf
caseinfo['exception'] = args.exception
caseinfo['run_info_str'] = args.runinfostr
caseinfo['qcrange_file'] = args.qcrangejsonfile
caseinfo['mindiscoveryvaf'] = float(args.minvaf)

if caseinfo['specimen'] == 'BM':
    caseinfo['specimen'] = 'Bone Marrow'
elif caseinfo['specimen'] == 'PB':
    caseinfo['specimen'] = 'Peripheral Blood'


#########################################
#
# Get Reference values for QC Metrics
#
#########################################

qcranges = {}
with open(caseinfo['qcrange_file'], 'r') as json_file:
    qcranges = json.load(json_file)

covLevel1 = int(qcranges['Coverage levels'][0])
covLevel2 = int(qcranges['Coverage levels'][1])
minTargetCov = float(qcranges['Target fraction at coverage'])


#########################################
#
# Get files from the case directory 
#
#########################################

vcffile = list(Path(caseinfo['casedir']).rglob('*.annotated.vcf.gz'))[0]
if not vcffile.is_file():
    sys.exit("VCF file " + str(vcffile) + " not valid.")

#haplotect = list(Path(caseinfo['casedir']).rglob('*.haplotect.txt'))[0]
#haplotectloci = list(Path(caseinfo['casedir']).rglob('*.haplotectloci.txt'))[0]
#if not haplotect.is_file() or not haplotectloci.is_file():
#    sys.exit("Haplotect output " + str(haplotect) + " not valid.")
    
mappingmetrics = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.mapping_metrics.csv'))[0]
targetmetrics = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.target_bed_coverage_metrics.csv'))[0]
coveragemetrics = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*_full_res.bed'))[0]
if not mappingmetrics.is_file() or not targetmetrics.is_file() or not coveragemetrics.is_file():
    sys.exit("DRAGEN metrics files not found.")

umimetrics = None
#umimetrics = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.umi_metrics.csv'))[0]

    
#########################################
#
# Set up dataframes for all data
#
#########################################

# this is for dragen/umi metrics
qcdf = pd.DataFrame(columns=['metric','value'])

# this is for exon/gene coverage metrics
covqcdf = pd.DataFrame(columns=['Gene','Type','Region','Mean','covLevel1','covLevel2'])

# dataframe with all variants
variants = pd.DataFrame(columns=['type','chrom','pos','ref','alt','genotype','gene','transcript','consequence','csyntax','psyntax','exon','popaf','annotations','coverage','altreads','vaf'])


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
# Collect QC metrics
#
#########################################

print("Collecting DRAGEN qc metrics...",file=sys.stderr)

# read in mapping metrics
df = pd.read_csv(mappingmetrics,sep=',',names=['group','readgroup','metric','value','percent'])
df = df[df['group']=='MAPPING/ALIGNING SUMMARY'].drop(columns='readgroup')
df['metric'] = df['group'] + ': ' + df['metric']
df = df.drop(columns='group')
dfpct = df[df['percent']==df['percent']].copy()
dfpct['metric'] = dfpct['metric'].apply(lambda x: x + ' (%)')
dfpct['value'] = dfpct['percent']
dfpct = dfpct.drop(columns='percent')
qcdf = pd.concat([qcdf,df])
qcdf = pd.concat([qcdf,dfpct])

# read in umi metrics
if umimetrics is not None:
    df = pd.read_csv(umimetrics,sep=',',names=['group','readgroup','metric','value','percent'])
    df['group'] = 'UMI SUMMARY'
    df = df.drop(columns='readgroup')
    df['metric'] = df['group'] + ': ' + df['metric']
    df = df.drop(columns='group')
    dfpct = df[df['percent']==df['percent']].copy()
    dfpct['metric'] = dfpct['metric'].apply(lambda x: x + ' (%)')
    dfpct['value'] = dfpct['percent']
    dfpct = dfpct.drop(columns='percent')
    
    qcdf = pd.concat([qcdf,df])
    qcdf = pd.concat([qcdf,dfpct])

    
# read in target metrics
df = pd.read_csv(targetmetrics,sep=',',names=['group','readgroup','metric','value','percent'])
df = df.drop(columns='readgroup')
df['metric'] = df['group'] + ': ' + df['metric']
df = df.drop(columns='group')
dfpct = df[df['percent']==df['percent']].copy()
dfpct['metric'] = dfpct['metric'].apply(lambda x: x + ' (%)')
dfpct['value'] = dfpct['percent']
dfpct = dfpct.drop(columns='percent')
qcdf = pd.concat([qcdf,df])
qcdf = pd.concat([qcdf,dfpct])

# get coverage info

print("Collecting coverage metrics...",file=sys.stderr)

# intersect full res coverage bedfile w/ coverage QC bed file to calculate coverage
covQcBedPr = pr.PyRanges(pd.read_csv(caseinfo['covqcbedfile'], header=None, names="Chromosome Start End Exon Strand Gene len geneId transcriptId".split(), sep="\t"))
fullResCovPr = pr.PyRanges(pd.read_csv(coveragemetrics, header=None, names="Chromosome Start End cov".split(), sep="\t"))
df = covQcBedPr.join(fullResCovPr).df
df['nt'] = df[['End','End_b']].min(axis=1) - df[['Start','Start_b']].max(axis=1)
df['tcov'] = df['nt'] * df['cov']

# get hotspot qc
hotspotdf = df[(df.Exon=='HOTSPOTQC')].copy()
hotspotdf['Region'] = 'codon_' + hotspotdf['len']
hotspotdf['Gene'] = hotspotdf['Strand']
x = hotspotdf[['Gene','Region']].join(hotspotdf.groupby(['Gene','Region'])[['cov']].min(),on=['Gene','Region'])
x['Mean'] = x['cov']
x['covLevel1'] = x['cov'] > covLevel1
x['covLevel1'] = x['covLevel1'].replace({True:1,False:0})
x['covLevel2'] = x['cov'] > covLevel2
x['covLevel2'] = x['covLevel2'].replace({True:1,False:0})
x['Type'] = 'hotspot'
x = x.drop_duplicates()

covqcdf = pd.concat([covqcdf,x[['Gene','Type','Region','Mean','covLevel1','covLevel2']]])

# get Gene and exon qc
df = df[df['Exon']!='HOTSPOTQC'].copy()
df['len'] = df['len'].astype(int)

# get list of transcripts
assaygenelist = df[['Gene','geneId','transcriptId']].drop_duplicates()

# exon qc
exonqcdf = df.groupby(['Gene','Exon','len']).sum()[['tcov']].reset_index()
exonqcdf['Mean'] = exonqcdf['tcov']/exonqcdf['len']
x = df[(df['cov']>=covLevel1)].groupby(['Gene','Exon','len']).sum().reset_index()
x['covLevel1']= x.nt/x.len*100
exonqcdf = pd.merge(exonqcdf,x[['Gene','Exon','len','covLevel1']],on=['Gene','Exon','len'],how='left')
x = df[(df['cov']>=covLevel2)].groupby(['Gene','Exon','len']).sum().reset_index()
x['covLevel2']= x.nt/x.len*100
exonqcdf = pd.merge(exonqcdf,x[['Gene','Exon','len','covLevel2']],on=['Gene','Exon','len'],how='left')
exonqcdf['Type'] = 'Exon'
exonqcdf['Region'] = exonqcdf['Exon']
exonqcdf.fillna(0, inplace=True)

covqcdf = pd.concat([covqcdf,exonqcdf[['Gene','Type','Region','Mean','covLevel1','covLevel2']]])

# get Gene qc
#
geneqcdf = df[['Gene','Start','End','len']].drop_duplicates().groupby(['Gene']).sum()[['len']].reset_index()
geneqcdf = pd.merge(geneqcdf,df.groupby(['Gene']).sum()[['tcov']].reset_index(),on=['Gene'],how='left')
geneqcdf['Mean'] = geneqcdf['tcov']/geneqcdf['len']
geneqcdf = pd.merge(geneqcdf,df[(df['cov']>=covLevel1)].groupby(['Gene']).sum().reset_index(),on=['Gene'],how='left',suffixes=['','_y'])
geneqcdf['covLevel1'] = geneqcdf['nt']/geneqcdf['len']*100
geneqcdf = geneqcdf[['Gene','len','Mean','covLevel1']]
geneqcdf = pd.merge(geneqcdf,df[(df['cov']>=covLevel2)].groupby(['Gene']).sum().reset_index(),on=['Gene'],how='left',suffixes=['','_y'])
geneqcdf['covLevel2'] = geneqcdf['nt']/geneqcdf['len']*100
geneqcdf = geneqcdf[['Gene','len','Mean','covLevel1','covLevel2']]
geneqcdf['Type'] = 'Gene'
geneqcdf['Region'] = 'Gene'
geneqcdf.fillna(0, inplace=True)

covqcdf = pd.concat([covqcdf,geneqcdf[['Gene','Type','Region','Mean','covLevel1','covLevel2']]])

# low cov codons
lowcovdf = df[(df['cov']<covLevel1)].sort_values(['Gene','transcriptPos']).reset_index()
lowcovdf['Region'] = lowcovdf.apply(lambda row: [ pos2codon(int(row['CdsStart'])+1,int(row['CdsEnd']),int(row['transcriptPos']),x,row['Strand']) for x in range(max(int(row['
CdsStart']),int(row['Start_b']))+1,min(int(row['CdsEnd']),int(row['End_b']))+1)],axis=1)
lowcovdf = lowcovdf.explode('Region')[['Gene','Region','nt','tcov']]
lowcovdf = lowcovdf[(lowcovdf.Region.notnull())]
lowcovdf = lowcovdf.groupby(['Gene','Region']).sum().reset_index()
lowcovdf['Mean'] = lowcovdf['tcov']/lowcovdf['nt']
lowcovdf = pd.merge(lowcovdf.groupby(['Gene'])[['Mean']].mean().reset_index(),lowcovdf.groupby(['Gene'])['Region'].apply(list).apply(lambda x: make_ranges(x)).reset_index(),on='Gene')
lowcovdf['covLevel1'] = 0
lowcovdf['covLevel2'] = 0
lowcovdf['Type'] = 'Codon'

covqcdf = pd.concat([covqcdf,lowcovdf[['Gene','Type','Region','Mean','covLevel1','covLevel2']]])


covqcdf["Mean"] = covqcdf["Mean"].map(lambda x: float(x))
covqcdf["covLevel1"] = covqcdf["covLevel1"].map(lambda x: float(x))
covqcdf["covLevel2"] = covqcdf["covLevel2"].map(lambda x: float(x))
covqcdf.fillna(0)
    
# get haplotect output
#haplotectdf = pd.read_csv(haplotect,sep='\t')
#haplotectdf = haplotectdf.iloc[:, :-2]
#haplotectdf.fillna(0, inplace=True)

#haplotectlocidf = pd.read_csv(haplotectloci,sep='\t',skiprows=2)
#haplotectlocidf = haplotectlocidf.iloc[:, :-1]


#########################################
#
# Get variants
#
#########################################

print("Gathering variants...",file=sys.stderr)

vcf = VCF(vcffile)

# get VEP fields
vep = {}
i = 0
for j in vcf.get_header_type('CSQ')['Description'].split("|"):
    vep[j] = i
    i+=1

# get variants
for variant in vcf:
    
    vartype = ''
    if len(variant.REF) == len(variant.ALT[0]):
        vartype = 'SNV'
    else:
        vartype = 'INDEL'

    if variant.FILTER is not None and variant.FILTER!='PASS':
        continue

    varfilter = variant.FILTER
    
    abundance = round(variant.format('AF')[0][0] * 100,2)

    gt = [variant.REF] + variant.ALT
    genotype = gt[variant.genotypes[0][0]] + '/' + gt[variant.genotypes[0][1]]
            
    # get VEP annotation
    csq = variant.INFO['CSQ']
    
    if csq is None:
        sys.exit("No VEP fields")
    
    gene=''
    transcript=''
    csyntax='NA'
    psyntax='NA'
    consequence='NA'
    exon='NA'
    intron='NA'
    popmaf = 'NA'
    customannotation = 'NA'
    for i in variant.INFO['CSQ'].split(','):
        csq = i.split("|")
        # if this is the list of transcripts to use for annotation or if its not and its the 'picked' one'
        if csq[vep['Feature']] in list(assaygenelist['transcriptId']): 
            transcript = csq[vep['Feature']]
            gene = csq[vep['SYMBOL']]
            csyntax = csq[vep['HGVSc']].split(":")
            if len(csyntax) > 1:
                csyntax = csyntax[1]
            else:
                csyntax = 'noncoding'
    
            psyntax = csq[vep['HGVSp']].split(":")
            if len(psyntax) > 1:
                psyntax = convert_aa(psyntax[1])
                psyntax = re.sub("\%3D","=",psyntax)
            else:
                psyntax = csyntax

            consequence = csq[vep['Consequence']].split("&")[0]
            impact = csq[vep['IMPACT']]
            exon = csq[vep['EXON']] or 'NA'
            intron = csq[vep['INTRON']] or 'NA'
            customannotation = csq[vep['Existing_variation']] or 'NA'
            
            popmaf = 'NA'
            if csq[vep['MAX_AF']] != '':
                popmaf = float(csq[vep['MAX_AF']])


    # do final categorization

    variants = pd.concat([variants,pd.DataFrame([dict(zip(variants.columns,[vartype,str(variant.CHROM),str(variant.POS),variant.REF,variant.ALT[0],genotype,gene,transcript,consequence,csyntax,psyntax,exon,str(popmaf) + '%',customannotation,str(variant.format("DP")[0][0]),str(variant.format("AD")[0][0]),str(variant.format("AD")[0][1]),str(abundance)+"%"]))])])

print("Starting report...",file=sys.stderr)
    
#
# Start report
#

# make dict for report and redirect output for text report
jsonout = {'CASEINFO':{},'VARIANTS':{},'QC':{}}

f = open(caseinfo['name'] + ".report.txt", "w")
sys.stdout = f

dt = strftime("%Y-%m-%d %H:%M:%S", gmtime())

caseinfo['date'] = dt

print("GatewaySeq Report for " + caseinfo['name'] + " ---- Generated on: " + dt + "\n")

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

jsonout['CASEINFO'] = caseinfo

print("*** REPORTABLE MUTATIONS ***\n")

if variants.shape[0] > 0:
    print(variants.to_csv(sep='\t',header=True, index=False))
    jsonout['VARIANTS'] = variants.to_dict('split')
    del jsonout['VARIANTS']['index']
else:
    print("None Detected\n")

print("*** TIER 4 POLYMORPHISMS ***\n")

if variants[variants['category']=='SNP'].shape[0] > 0:
    print(variants[variants['category']=='SNP'].iloc[:,1:].to_csv(sep='\t',header=True, index=False))
    jsonout['VARIANTS']['SNPS'] = variants[variants['category']=='SNP'].iloc[:,1:].to_dict('split')
    del jsonout['VARIANTS']['SNPS']['index']
else:
    print("None Detected\n")

#varcats = variants['category'].value_counts().to_dict()

#print("*** VARIANT COUNTS BY CATEGORY ***\n")
#for l in ['Tier1-3','Low Level','Filtered','SNP','Genotyping','Silent/Not Reported']:
#    if l not in varcats.keys():
#        print(l + ":\t0")
#    else:
#        print(l + ":\t" + str(varcats[l]))

#jsonout['QC']['VARIANTCOUNTS'] = varcats

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

print("*** HOTSPOT QC ***\n")

xdf = covqcdf[(covqcdf.Type == "hotspot")][['Gene','Region','Mean']]
xdf = xdf.rename(columns={"Region":"Hotspot"})
xdf['QC'] = np.where(xdf['Mean'] < minTargetCov, '(!)', '')
print(xdf.to_csv(sep='\t',header=True, index=False,float_format='%.1f'))
jsonout['QC']['HOTSPOT QC'] = xdf.to_dict('split')
jsonout['QC']['HOTSPOT QC'].pop('index', None)

print("*** GENE COVERAGE QC ***\n")

xdf = covqcdf[(covqcdf.Type == "Gene")][['Gene','Mean','covLevel1','covLevel2']]
xdf['QC'] = np.where((xdf['Mean'] < covLevel2) | (xdf['covLevel1']<minTargetCov), '(!)', '')
xdf = xdf.rename(columns={"covLevel1": str(covLevel1)+"x", "covLevel2": str(covLevel2)+"x"})
print(xdf.to_csv(sep='\t',header=True, index=False,float_format='%.1f'))

jsonout['QC']['GENE COVERAGE QC'] = xdf.to_dict('split')
jsonout['QC']['GENE COVERAGE QC'].pop('index', None)

jsonout['QC']['FAILED GENES'] = ','.join(xdf[(xdf.QC!='')]['Gene'].tolist())
jsonout['QC']['FAILED GENE COUNT'] = xdf[(xdf.QC!='')].shape[0]

print("*** FAILED EXONS ***\n")

xdf = covqcdf[(covqcdf.Type == "Exon")][['Region','Mean','covLevel1','covLevel2']]
xdf['QC'] = np.where((xdf['Mean'] < covLevel2) | (xdf['covLevel1']<minTargetCov), '(!)', '')
xdf = xdf.rename(columns={"covLevel1": str(covLevel1)+"x", "covLevel2": str(covLevel2)+"x"})

jsonout['QC']['EXON COVERAGE QC'] = xdf.to_dict('split')
jsonout['QC']['EXON COVERAGE QC'].pop('index', None)

jsonout['QC']['FAILED EXONS'] = ','.join(xdf[(xdf.QC!='')]['Region'].tolist())
jsonout['QC']['FAILED EXON COUNT'] = xdf[(xdf.QC!='')]['Region'].shape[0]

if xdf[(xdf.QC!='')].shape[0] > 0:
    print(xdf[(xdf.QC!='')].to_csv(sep='\t',header=True, index=False,float_format='%.1f'))
else:
    print("NONE\n")

#print("*** Haplotect Contamination Estimate ***\n")

#print(haplotectdf.iloc[:,1:].to_csv(sep='\t',header=True, index=False))
#print(haplotectlocidf.to_csv(sep='\t',header=True, index=False))

#jsonout['QC']['HAPLOTECT SUMMARY'] = haplotectdf.iloc[:,1:].to_dict('split')
#jsonout['QC']['HAPLOTECT SUMMARY'].pop('index', None)
#jsonout['QC']['HAPLOTECT LOCI'] = haplotectlocidf.iloc[:,1:].to_dict('split')
#jsonout['QC']['HAPLOTECT LOCI'].pop('index', None)

print("*** GatewaySeq Assay Version " + str(qcranges["ASSAY VERSION"]) + " ***\n")

print(qcranges["DISCLAIMER"])

jsonout['QC']['QCINFO'] = qcranges

original_stdout = sys.stdout
f.close()

# dump json
j = open(caseinfo['name'] + ".report.json", "w")
json.dump(jsonout,j,indent=" ")
j.close()
