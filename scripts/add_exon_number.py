import argparse, sys
import pybedtools
import pandas as pd

def parse_gtf_attributes(attribute_string):
    attributes = {}
    fields = attribute_string.strip().split(';')
    for field in fields:
        if field.strip() != '':
            key, value = field.strip().split(' ', 1)
            attributes[key] = value.strip('\"')
    return attributes

def extract_exon_number(row):
    attr = {}
    if row['attr'] != '.' and row['category']=='GOAL_genes':
        attr = parse_gtf_attributes(row['attr'])

        if 'transcript_id' in attr.keys() and row['transcript_id']==attr['transcript_id'] and row['type']=='exon':
            return(attr['exon_number'])
        else:
            return('.')
    else:
        return('.')


def main():
    parser = argparse.ArgumentParser(description='Intersect BED file with GTF file and extract gene IDs.')
    parser.add_argument('bed_file', help='Path to the BED file')
    parser.add_argument('gtf_file', help='Path to the GTF file')
    args = parser.parse_args()

    bed_file = args.bed_file
    gtf_file = args.gtf_file

    bed = pybedtools.BedTool(bed_file)
    gtf = pybedtools.BedTool(gtf_file)

    intersected = bed.intersect(gtf, wao=True).to_dataframe(names=["chr","start","end","gene","info","chr_gtf","source","type","gtf_start","gtf_end","frame1","gtf_strand","frame2","attr","overlap"])
    infodf = intersected['info'].str.split('|',expand=True)
    infodf.columns = ['category','gene_id','transcript_id','exon_number','cdsstart','cdsend','strand']
    intersected = pd.concat([intersected, infodf], axis=1)

    intersected['exon_number'] = intersected.apply(lambda r: extract_exon_number(r), axis=1)
    intersected = intersected.sort_values(['chr','start','end','gene','info','exon_number']).reindex().drop_duplicates(subset=['chr','start','end','gene','info'],keep='last')
    intersected['exon_number'] = intersected.apply(lambda r: 'region_'+str(r['start'])+'-'+str(r['end']) if r['exon_number']=='.' else 'exon_'+r['exon_number'],axis=1)
    intersected['info'] = intersected[['category','gene_id','transcript_id','exon_number','cdsstart','cdsend','strand']].apply(lambda r: "|".join(r),axis=1)
    x  = intersected[['chr','start','end','gene','info']].drop_duplicates()
    x.to_csv("~/test.bed",sep="\t", index=False,header=False)


if __name__ == '__main__':
    main()
