#!/usr/bin/env python

import sys, os, re, tempfile, csv, pysam, json, binascii, argparse
import pandas as pd
import pyranges as pr
import numpy as np
from time import gmtime, strftime
from cyvcf2 import VCF
from pathlib import Path

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
def testinrange(row,loss,gain):
    if row['Segment_Call']!='0' and ((row['Segment_Name'] in loss and row['Segment_Call']!='-') or (row['Segment_Name'] in gain and row['Segment_Call']!='+')):
        return ';'.join([row['Filter'],'WrongSign'])
    else:
        return row['Filter']

#
# Script
#

parser = argparse.ArgumentParser(description='Make GatewaySeq report')
parser.add_argument('-n','--name',required=True,help='Sample name')
parser.add_argument('-d','--dir',required=True,help='Output directory')
parser.add_argument('-l','--fusionlist',required=True,help='Gene pairs to report for fusions')
parser.add_argument('-q','--qcrangejsonfile',required=True,help='QCReferenceRanges.json')
parser.add_argument('-m','--mrn',default='NONE',help='Sample MRN number')
parser.add_argument('-a','--accession',default='NONE',help='Sample accession number')
parser.add_argument('-s','--specimen',default='NONE',help='Sample specimen type')
parser.add_argument('-b','--DOB',default='NONE',help='Date of birth')
parser.add_argument('-e','--exception',default='NONE',help='Exception')
parser.add_argument('-f','--minvaf',default=5.0,help='Minimum validated VAF')
parser.add_argument('-r','--minreads',default=5.0,help='Minimum alt reads to report variant')
parser.add_argument('-x','--excludevaf',default=1.0,help='Minimum VAF for reporting, below this will not be shown')
parser.add_argument('-i','--runinfostr',default='NONE',help='Illumina Run Information String')
parser.add_argument('-p','--maxaf',default=0.001,help='Maximum population allele frequency for potential somatic variants')

args = parser.parse_args()

caseinfo = {}
caseinfo['name'] = args.name
caseinfo['mrn'] = args.mrn
caseinfo['DOB'] = args.DOB
caseinfo['accession'] = args.accession
caseinfo['specimen'] = args.specimen
caseinfo['casedir'] = args.dir
caseinfo['maxaf'] = float(args.maxaf)
caseinfo['exception'] = args.exception
caseinfo['run_info_str'] = args.runinfostr
caseinfo['qcrange_file'] = args.qcrangejsonfile
caseinfo['minvaf'] = float(args.minvaf)
caseinfo['minreads'] = float(args.minreads)
caseinfo['genefusions'] = args.fusionlist

excludevaf = float(args.excludevaf)

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

minFractionCovered = float(qcranges['Target fraction at coverage'].split(',')[0])
minTargetCoverage = [ int(x) for x in qcranges['Minimum target coverage'].split(',') ]
minCoverage = float(qcranges['COVERAGE SUMMARY: Average alignment coverage over target region'].split(',')[0])

nonSynon = ['splice_acceptor_variant','splice_donor_variant','stop_gained','frameshift_variant','stop_lost','start_lost','transcript_amplification','inframe_insertion','inframe_deletion','missense_variant','protein_altering_variant']
nonCodingList = qcranges['NONCODINGVARIANTLIST']

#########################################
#
# Get files from the case directory 
#
#########################################

vcffile = list(Path(caseinfo['casedir']).rglob('*.annotated.vcf.gz'))[0]
if not vcffile.is_file():
    sys.exit("VCF file " + str(vcffile) + " not valid.")

#svvcffile = list(Path("/storage1/fs1/duncavagee/Active/SEQ/GatewaySeq/process/hackathon/svs/").rglob('*-Gateway-Seq-S15-27904-lib1.sv.annotated2.vcf.gz'))[0] #list(Path(caseinfo['casedir']).rglob('*.sv_annotated.vcf.gz'))[0]
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

mappingmetrics = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.mapping_metrics.csv'))[0]
targetmetrics = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.target_bed_coverage_metrics.csv'))[0]
umimetrics = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.umi_metrics.csv'))[0]
genetargetbed = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.target_bed_read_cov_report.bed'))[0]
alltargetbed = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.qc-coverage-region-1_read_cov_report.bed'))[0]
coveragebed = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.qc-coverage-region-1_full_res.bed'))[0]
MSIjson = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.microsat_output.json'))[0]
TMBcsv = list(Path(os.path.join(caseinfo['casedir'],"dragen")).rglob('*.tmb.metrics.csv'))[0]

if not mappingmetrics.is_file() or not targetmetrics.is_file() or not umimetrics.is_file() or not genetargetbed.is_file() or not alltargetbed.is_file() or not coveragebed.is_file() or not MSIjson.is_file() or not TMBcsv.is_file():
    sys.exit("DRAGEN metrics files not found.")

#########################################
#
# Set up dataframes for all data
#
#########################################

# this is for dragen/umi metrics
qcdf = pd.DataFrame(columns=['metric','value'])

# this is for exon/gene coverage metrics
covqcdf = pd.DataFrame(columns=['Gene','Type','Region','Mean','Covered1','Covered2'])

# dataframe with all variants
variants = pd.DataFrame(columns=['category','type','filter','chrom','pos','ref','alt','gene','transcript','consequence','csyntax','psyntax','exon','popaf','annotations','coverage','altreads','vaf'])

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
df = pd.read_csv(umimetrics,sep=',',names=['group','readgroup','metric','value','percent'])
df['group'] = 'UMI SUMMARY'
df = df.drop(columns='readgroup')
df['metric'] = df['group'] + ': ' + df['metric']
df = df.drop(columns='group')
dfpct = df[df['percent']==df['percent']].copy()
dfpct['metric'] = dfpct['metric'].apply(lambda x: x + ' (%)')
dfpct['value'] = dfpct['percent']
dfpct = dfpct.drop(columns='percent')
consensusReads = round(df.loc[df.metric=='UMI SUMMARY: Consensus pairs emitted','value'].astype(int).tolist()[0] * 2 / df.loc[df.metric=='UMI SUMMARY: Number of reads','value'].astype(int).tolist()[0] * 100,1)
duplicateReads = round(100-consensusReads,1)

dfpct = pd.concat([dfpct,pd.DataFrame.from_dict({0:['UMI SUMMARY: Consensus reads (%)',consensusReads],1:['UMI SUMMARY: Duplicate reads (%)',duplicateReads]},orient='index',columns=['metric','value'])],axis=0)

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

# get coverage file
fullResCovPr = pr.PyRanges(pd.read_csv(coveragebed, header=None, names="Chromosome Start End Coverage".split(), sep="\t"))

ntAtCovlevel1 = sum(fullResCovPr[fullResCovPr.Coverage>minTargetCoverage[0]].lengths())
ntAtCovlevel2 = sum(fullResCovPr[fullResCovPr.Coverage>minTargetCoverage[1]].lengths())

qcdf.loc[len(qcdf.index)] = ['COVERAGE SUMMARY: Target at ' + str(minTargetCoverage[0]) + 'x (%)', ntAtCovlevel1/sum(fullResCovPr.lengths())*100, '']
qcdf.loc[len(qcdf.index)] = ['COVERAGE SUMMARY: Target at ' + str(minTargetCoverage[1]) + 'x (%)', ntAtCovlevel2/sum(fullResCovPr.lengths())*100, '']

# get gene coverage info

# intersect full res coverage bedfile w/ coverage QC bed file to calculate coverage
covBedPr = pr.PyRanges(pd.read_csv(genetargetbed, skiprows=1, header=None, names="Chromosome Start End Gene Region ReadCov R1Cov R2Cov".split(), sep="\t"))
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

# intersect full res coverage bedfile w/ coverage QC bed file to calculate coverage
covBedPr = pr.PyRanges(pd.read_csv(alltargetbed, skiprows=1, header=None, names="Chromosome Start End Gene Region ReadCov R1Cov R2Cov".split(), sep="\t"))

# note, this file has all the transcript IDs, so get those now
transcripts = covBedPr.df.Region.str.split('|',expand=True)[[2]].replace('\s\+\s\S+','',regex=True).drop_duplicates()[2].tolist()

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

# get haplotect output
haplotectdf = pd.read_csv(haplotect,sep='\t')
haplotectdf = haplotectdf.iloc[:, :-2]
haplotectdf.fillna(0, inplace=True)

haplotectlocidf = pd.read_csv(haplotectloci,sep='\t',skiprows=2)
haplotectlocidf = haplotectlocidf.iloc[:, :-1]

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
for j in vcf.get_header_type('CSQ')['Description'].split("|"):
    vep[j] = i
    i+=1

# get variants
for variant in vcf:

    cat = 'Tier1-3'
    
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
    if abundance < excludevaf:
        continue

    if abundance < caseinfo['minvaf'] or variant.format('AD')[0][1] < caseinfo['minreads'] and varfilter=='PASS':
        varfilter = 'LowReads'
        
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
        consequence = csq[vep['Consequence']].split("&")[0]
        # if this is the list of transcripts to use for annotation or if its not and its the 'picked' one'
        if csq[vep['Feature']] in transcripts or (transcript=='' and csq[vep['PICK']] == '1'):
            transcript = csq[vep['Feature']]
            gene = csq[vep['SYMBOL']]

            if consequence in nonSynon or (gene in nonCodingList.keys() and consequence in nonCodingList[gene]):

                csyntax = csq[vep['HGVSc']].split(":")
                if len(csyntax) > 1:
                    csyntax = csyntax[1]
                else:
                    if 'upstream' in consequence:
                        if csq[vep['STRAND']]==1:
                            csyntax = "c.-"+str(csq[vep['DISTANCE']])+variant.REF+'>'+csq[0]
                        else:
                            csyntax = "c.-"+str(csq[vep['DISTANCE']])+revcomp(variant.REF)+'>'+revcomp(csq[0])

                    elif 'downstream' in consequence:
                        if csq[vep['STRAND']]==1:
                            csyntax = "c.+"+str(csq[vep['DISTANCE']])+variant.REF+'>'+csq[0]
                        else:
                            csyntax = "c.+"+str(csq[vep['DISTANCE']])+revcomp(variant.REF)+'>'+revcomp(csq[0])
                    else:
                        csyntax = 'noncoding'

                psyntax = csq[vep['HGVSp']].split(":")
                if len(psyntax) > 1:
                    psyntax = convert_aa(psyntax[1])
                    psyntax = re.sub("\%3D","=",psyntax)
                else:
                    psyntax = consequence
            
                impact = csq[vep['IMPACT']]
                exon = csq[vep['EXON']] or 'NA'
                intron = csq[vep['INTRON']] or 'NA'
                customannotation = csq[vep['Existing_variation']] or 'NA'            

                popmaf = 'NA'
                if csq[vep['MAX_AF']] != '':
                    popmaf = float(csq[vep['MAX_AF']])

                if popmaf!='NA' and popmaf >= caseinfo['maxaf']:
                    cat = 'SNP'

                elif varfilter!='PASS':
                    cat = 'Filtered'
                
                else:
                    cat = 'Tier1-3'

                variants = pd.concat([variants,pd.DataFrame([dict(zip(variants.columns,[cat,vartype,varfilter,str(variant.CHROM),str(variant.POS),variant.REF,variant.ALT[0],gene,transcript,consequence,csyntax,psyntax,exon,str(popmaf) + '%',customannotation,str(variant.format("DP")[0][0]),str(variant.format("AD")[0][1]),str(abundance)+"%"]))])])


#
# Get CNVs
#

cnv = pd.read_csv(cnvsegoutput,delimiter='\t')
cnvcalls = pd.concat([cnv[cnv['Segment_Name'].isin(qcranges['CNVTARGETS']['loss'])],cnv[cnv['Segment_Name'].isin(qcranges['CNVTARGETS']['gain'])]])

cnvcalls['Filter'] = cnvcalls.apply(lambda r: testinrange(r,qcranges['CNVTARGETS']['loss'],qcranges['CNVTARGETS']['gain']), axis=1)
cnvcalls['type'] = ''
cnvcalls.loc[cnvcalls['Segment_Call']=='-','type'] = 'LOSS'
cnvcalls.loc[cnvcalls['Segment_Call']=='+','type'] = 'GAIN'
cnvcalls.loc[cnvcalls['Segment_Call']=='0','type'] = 'REF'
cnvcalls = cnvcalls[['type','Segment_Name','Chromosome','Start','End','Segment_Mean','Qual','Filter','Copy_Number','Ploidy']].copy()
cnvcalls.rename(columns={'Segment_Name':'gene','Chromosome':'chrom','Start':'start','End':'end','Segment_Mean':'copy_ratio','Qual':'qual','Filter':'filter','Copy_Number':'copynumber','Ploidy':'ploidy'},inplace=True)
print("Starting report...",file=sys.stderr)


#
# Get SVs
# 

svvcf = VCF(svvcffile)

# get known fusion events
knownsvs = pd.read_csv(caseinfo['genefusions'],sep='\t', header=None)
knownsvs.columns = ['gene1','gene2']

svs = pd.DataFrame(columns=['type','psyntax','csyntax','filter','chr1','pos1','chr2','pos2','gene1','intron1','strand1','gene2','intron2','strand2','id','abundance','info'])

# get VEP fields
vep = {}
i = 0
for j in svvcf.get_header_type('CSQ')['Description'].split("|"):
    vep[j] = i
    i+=1

passedvars = {}
alreadydone = set()

for variant in svvcf:
    if variant.INFO['SVTYPE']!='BND':
        continue

    passedvars[variant.ID] = variant

for v in passedvars.items():

    variant = v[1]
    isknown = 0
    knowngenes = []
    mate = ''

    filter = 'PASS'
    if variant.FILTER is not None:
        filter = variant.FILTER

    # get gene 1 info
    gene1='NA'
    transcript1='NA'
    intron1 = 'NA'
    strand1 = '+'
    for i in variant.INFO['CSQ'].split(','):
        csq = i.split("|")
        consequence = csq[vep['Consequence']].split("&")[0]
        # if this is the list of transcripts to use for annotation or if its not and its the 'picked' one'
        if csq[vep['Feature']] in transcripts: # or (transcript1=='' and csq[vep['PICK']] == '1'):
            transcript1 = csq[vep['Feature']]
            gene1 = csq[vep['SYMBOL']]
            intron1 = csq[vep['INTRON']].split('/')[0]
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
    intron2 = 'NA'
    strand2 = '+'
    for i in mate.INFO['CSQ'].split(','):
        csq = i.split("|")
        consequence = csq[vep['Consequence']].split("&")[0]
        # if this is the list of transcripts to use for annotation or if its not and its the 'picked' one'
        if csq[vep['Feature']] in transcripts: # or (transcript2=='' and csq[vep['PICK']] == '1'):
            transcript2 = csq[vep['Feature']]
            gene2 = csq[vep['SYMBOL']]
            intron2 = csq[vep['INTRON']].split('/')[0]
            if csq[vep['STRAND']]==-1:
                strand2 = '-'

    chr2 = mate.CHROM
    pos2 = mate.POS
        
    if gene1=='NA':
        gene1 = 'INTERGENIC'
        
    if gene2=='NA':
        gene2 = 'INTERGENIC'
    
    # one gene must be a targeted one
    if any([gene1 in x for x in svgenes]) is False and any([gene2 in x for x in svgenes]) is False:
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

    if sr[1] <= 10 or abundance < 5.0:
        continue

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
    svs.loc[len(svs.index)] = ['BND',psyntax,csyntax,filter,chr1,str(pos1),chr2,str(pos2),gene1,intron1,strand1,gene2,intron2,strand2,str(variant.ID) + ";" + str(mate.ID),str(round(abundance,1))+"%",infostring]

    alreadydone.add(variant.ID)



#
# Start report
#

# make dict for report and redirect output for text report
jsonout = {'CASEINFO':{},'VARIANTS':{},'CNV':{},'FUSIONS':{},'QC':{},'MSI':{},'TMB':{}}

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

print()

jsonout['CASEINFO'] = caseinfo

print("*** GENE MUTATIONS ***\n")

if variants[variants['category']=='Tier1-3'].shape[0] > 0:
    print(variants[variants['category']=='Tier1-3'].iloc[:,1:].to_csv(sep='\t',header=True, index=False))
    jsonout['VARIANTS']['Tier1-3'] = variants[variants['category']=='Tier1-3'].iloc[:,1:].to_dict('split')
    del jsonout['VARIANTS']['Tier1-3']['index']
else:
    print("None Detected\n")

print("*** FILTERED GENE MUTATIONS ***\n")

if variants[variants['category']=='Filtered'].shape[0] > 0:
    print(variants[variants['category']=='Filtered'].iloc[:,1:].to_csv(sep='\t',header=True, index=False))
    jsonout['VARIANTS']['Filtered'] = variants[variants['category']=='Filtered'].iloc[:,1:].to_dict('split')
    del jsonout['VARIANTS']['Filtered']['index']
else:
    print("None Detected\n")

print("*** COPY NUMBER VARIANTS ***\n")

if cnvcalls.shape[0] > 0:
    print(cnvcalls[cnvcalls['filter']=='PASS'].to_csv(sep='\t',header=True, index=False))
    jsonout['CNV']['PASS'] = cnvcalls[cnvcalls['filter']=='PASS'].to_dict('split')
    del jsonout['CNV']['PASS']['index']
else:
    print("None Detected\n")

print("*** FILTERED COPY NUMBER VARIANTS ***\n")

if cnvcalls.shape[0] > 0:
    print(cnvcalls[cnvcalls['filter']!='PASS'].to_csv(sep='\t',header=True, index=False))
    jsonout['CNV']['Filtered'] = cnvcalls[cnvcalls['filter']!='PASS'].to_dict('split')
    del jsonout['CNV']['Filtered']['index']
else:
    print("None Detected\n")

print("*** GENE FUSION VARIANTS ***\n")

if svs[svs['filter']=='PASS'].shape[0] > 0:
    print(svs[svs['filter']=='PASS'].to_csv(sep='\t',header=True, index=False))
    jsonout['FUSIONS']['PASS'] = svs[svs['filter']=='PASS'].to_dict('split')
    del jsonout['FUSIONS']['PASS']['index']
else:
    print("None Detected\n")

print("*** FILTERED GENE FUSION VARIANTS ***\n")

if svs[svs['filter']!='PASS'].shape[0] > 0:
    print(svs[svs['filter']!='PASS'].to_csv(sep='\t',header=True, index=False))
    jsonout['FUSIONS']['Filtered'] = svs[svs['filter']!='PASS'].to_dict('split')
    del jsonout['FUSIONS']['Filtered']['index']
else:
    print("None Detected\n")

print("*** MSI ***\n")
with open(MSIjson) as msi_json:
    MSI_json = json.load(msi_json)
msi_json.close()
print("PercentageUnstableSites: " + MSI_json["PercentageUnstableSites"])
print("ResultIsValid: " + MSI_json["ResultIsValid"])
print("SumDistance: " + MSI_json["SumDistance"])
jsonout['MSI']['PercentageUnstableSites'] = MSI_json['PercentageUnstableSites']
jsonout['MSI']['ResultIsValid'] = MSI_json['ResultIsValid']
jsonout['MSI']['SumDistance'] = MSI_json['SumDistance']

print("\n*** TMB ***\n")
with open(TMBcsv, 'r') as tmb_csv:
    tmb_reader = csv.reader(tmb_csv, delimiter=',')
    for tmb_line in tmb_reader:
        if tmb_line[2] == 'TMB':
            TMB = tmb_line[3]
        if tmb_line[2] == 'Nonsyn TMB':
            NS_TMB = tmb_line[3]
tmb_csv.close()
print('TMB: '+ TMB)
print('Nonsyn TMB: ' + NS_TMB)
jsonout['TMB']['TMB'] = TMB
jsonout['TMB']['Nonsyn TMB'] = NS_TMB

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

print("*** GatewaySeq Assay Version " + str(qcranges["ASSAY VERSION"]) + " ***\n")

print(qcranges["DISCLAIMER"])

jsonout['QC']['QCINFO'] = qcranges

original_stdout = sys.stdout
f.close()

# dump json
j = open(caseinfo['name'] + ".report.json", "w")
json.dump(jsonout,j,indent=" ")
j.close()
