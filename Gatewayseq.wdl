import "GatewayseqAnalysis.wdl" as subWF

workflow Gatewayseq {

    File SampleSheet
    # sample sheet has this structure:
    # index  name  RG_ID  RG_FLOWCELL  RG_LANE  RG_LIB  RG_SAMPLE [R1] [R2]
    
    File? DemuxSampleSheet

    Array[String] Adapters = ["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC","AGATCGGAAGAGCGTCGTGTAGGGAAA"]

    String? IlluminaDir
    String JobGroup
    String OutputDir

    String Queue
    String DragenQueue = "duncavagee"

    String DragenReference = "/staging/runs/Chromoseq/refdata/dragen_hg38"
    String Reference    = "/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/refdata/hg38/all_sequences.fa"
    String ReferenceDict = "/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/refdata/hg38/all_sequences.dict"

    String VEP          = "/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/VEP_cache"
    String QcMetrics    = "/storage1/fs1/gtac-mgi/Active/CLE/analysis/gatewayseq/git/cle-gatewayseq/accessory_files/GatewaySeqMetrics.txt"
    String Description  = "/storage1/fs1/gtac-mgi/Active/CLE/analysis/gatewayseq/git/cle-gatewayseq/accessory_files/GatewaySeqDescription.txt" 

    String HaplotectBed = "/storage1/fs1/gtac-mgi/Active/CLE/analysis/gatewayseq/git/cle-gatewayseq/accessory_files/myeloseq.haplotect_snppairs_hg38.041718.bed"
    String AmpliconBed  = "/storage1/fs1/gtac-mgi/Active/CLE/analysis/gatewayseq/git/cle-gatewayseq/accessory_files/GatewaySeq66650-1622560509.Amplicons.hg38.110321.bed"
    String CoverageBed  = "/storage1/fs1/gtac-mgi/Active/CLE/analysis/gatewayseq/git/cle-gatewayseq/accessory_files/GatewaySeq66650-1622560509.CoverageQC.hg38.110321.bed"
    String DragenCoverageBed = "/staging/runs/Haloplex/dragen_align_inputs/GatewaySeq66650-1622560509.CoverageQC.hg38.110321.bed"

    String CustomAnnotationVcf   = "/storage1/fs1/gtac-mgi/Active/CLE/analysis/gatewayseq/git/cle-gatewayseq/accessory_files/myeloseq_custom_annotations.annotated.011618.hg38.vcf.gz"
    String CustomAnnotationIndex = "/storage1/fs1/gtac-mgi/Active/CLE/analysis/gatewayseq/git/cle-gatewayseq/accessory_files/myeloseq_custom_annotations.annotated.011618.hg38.vcf.gz.tbi"
    String CustomAnnotationParameters = "MYELOSEQ,vcf,exact,0,TCGA_AC,MDS_AC,MYELOSEQBLACKLIST"

    String QC_pl = "/storage1/fs1/gtac-mgi/Active/CLE/analysis/gatewayseq/git/cle-gatewayseq/QC_metrics.pl"
    String DemuxFastqDir = "/scratch1/fs1/gtac-mgi/CLE/gatewayseq/demux_fastq"

    if (defined(DemuxSampleSheet)){
      call dragen_demux {
        input: Dir=IlluminaDir,
        OutputDir=OutputDir,
        SampleSheet=DemuxSampleSheet,
        queue=DragenQueue,
        jobGroup=JobGroup
      }

      call prepare_samples {
        input: SampleSheet=SampleSheet,
        Fastq1=dragen_demux.read1,
        Fastq2=dragen_demux.read2,
        queue=Queue,
        jobGroup=JobGroup
      }
    }

    Array[Array[String]] inputData = read_tsv(select_first([prepare_samples.sample_sheet,SampleSheet]))

    # the inputdata should be: index  name  RG_ID  RG_FLOWCELL  RG_LANE  RG_LIB  RG_SAMPLE read1path read2path
    scatter (samples in inputData){

        if(!defined(DemuxSampleSheet)){
          call trim_reads {
              input: Read1=samples[7],
              Read2=samples[8],
              Adapters=Adapters,
              Name=samples[1],
              queue=Queue,
              jobGroup=JobGroup
          }
        }

        call dragen_align {
            input: DragenRef=DragenReference,
                   fastq1=select_first([trim_reads.read1,samples[7]]),
                   fastq2=select_first([trim_reads.read2,samples[8]]),
                   Name=samples[1],
                   RG=samples[3] + '.' + samples[4] + '.' + samples[0],
                   SM=samples[6],
                   LB=samples[5] + '.' + samples[0],
                   AmpliconBed=AmpliconBed,
                   CoverageBed=DragenCoverageBed,
                   OutputDir=OutputDir,
                   SubDir=samples[1] + '_' + samples[0],
                   queue=DragenQueue,
                   jobGroup=JobGroup
        }

        call convert_bam {
            input: Bam=dragen_align.bam,
                   BamIndex=dragen_align.bai,
                   AmpliconBed=AmpliconBed,
                   refFasta=Reference,
                   Name=samples[1],
                   OutputDir=OutputDir,
                   SubDir=samples[1] + '_' + samples[0],
                   queue=Queue,
                   jobGroup=JobGroup
        }

        call subWF.GatewayseqAnalysis {
            input: Cram=convert_bam.cram,
                   CramIndex=convert_bam.crai,
                   CoverageBed=CoverageBed,
                   refFasta=Reference,
                   ReferenceDict=ReferenceDict,
                   Name=samples[1],
                   Vepcache=VEP,
                   CustomAnnotationVcf=CustomAnnotationVcf,
                   CustomAnnotationIndex=CustomAnnotationIndex,
                   CustomAnnotationParameters=CustomAnnotationParameters,
                   HaplotectBed=HaplotectBed,
                   QcMetrics=QcMetrics,
                   Description=Description,
                   OutputDir=OutputDir,
                   SubDir=samples[1] + '_' + samples[0],
                   Queue=Queue,
                   JobGroup=JobGroup
        }
    } 
    
    if (defined(DemuxSampleSheet)){
        call move_demux_fastq {
            input: order_by=GatewayseqAnalysis.all_done,
            Batch=basename(OutputDir),
            DemuxFastqDir=DemuxFastqDir,
            queue=DragenQueue,
            jobGroup=JobGroup
        }

        call remove_rundir {
            input: order_by=GatewayseqAnalysis.all_done,
            rundir=IlluminaDir,
            queue=DragenQueue,
            jobGroup=JobGroup
        }
    }

    call batch_qc {
        input: order_by=GatewayseqAnalysis.all_done,
               BatchDir=OutputDir,
               QC_pl=QC_pl,
               queue=Queue,
               jobGroup=JobGroup
    }
}


task dragen_demux {
     String Dir
     String OutputDir
     String SampleSheet
     String jobGroup
     String queue

     String batch = basename(OutputDir)
     String StagingDir = "/staging/runs/Haloplex/"
     String LocalFastqDir = StagingDir + "demux_fastq/" + batch
     String LocalReportDir = LocalFastqDir + "/Reports"
     String LocalSampleSheet = StagingDir + "sample_sheet/" + batch + '.csv'
     String log = StagingDir + "log/" + batch + "_demux.log"
     String DemuxReportDir = OutputDir + "/dragen_demux_reports"

     command <<<
         /bin/cp ${SampleSheet} ${LocalSampleSheet} && \
         /opt/edico/bin/dragen --bcl-conversion-only true --bcl-only-matched-reads true --strict-mode true --sample-sheet ${LocalSampleSheet} --bcl-input-directory ${Dir} --output-directory ${LocalFastqDir} &> ${log} && \
         /bin/ls ${LocalFastqDir}/*_R1_001.fastq.gz > Read1_list.txt && \
         /bin/ls ${LocalFastqDir}/*_R2_001.fastq.gz > Read2_list.txt && \
         /bin/mv ${log} ./ && \
         /bin/rm -f ${LocalSampleSheet} && \
         /bin/cp -r ${LocalReportDir} ${DemuxReportDir}
     >>>

     runtime {
         docker_image: "seqfu/centos7-dragen-3.9.3:latest"
         cpu: "20"
         memory: "200 G"
         queue: queue
#        host: "compute1-dragen-3"
         job_group: jobGroup 
     }
     output {
         File read1 = "Read1_list.txt"
         File read2 = "Read2_list.txt"
     }
}

task prepare_samples {
     File SampleSheet
     String Fastq1
     String Fastq2
     String jobGroup
     String queue

     command <<<
             /bin/cp ${Fastq1} 1.tmp.txt
             /bin/cp ${Fastq2} 2.tmp.txt
             /usr/bin/perl -e 'open(R1,"1.tmp.txt"); @r1 = <R1>; \
                 chomp @r1; close R1;\
                 open(R2,"2.tmp.txt"); @r2 = <R2>; \
                 chomp @r2; close R2; \
                 open(SS,"${SampleSheet}");
                 while(<SS>){
                     chomp;
                     my @l = split("\t",$_);
                     my $s = $l[1].'_';
                     my $r1 = (grep /$s/, @r1)[0];
                     my $r2 = (grep /$s/, @r2)[0];
                     print join("\t",@l,$r1,$r2),"\n";
                 }
                 close SS;' > sample_sheet.txt
     >>>
     runtime {
         docker_image: "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
         cpu: "1"
         memory: "4 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File sample_sheet = "sample_sheet.txt"
     }
}

task trim_reads {
     String Read1
     String Read2
     Array[String] Adapters
     String Name
     String jobGroup
     String queue

     command {
         if [[ $Read1 =~ _R1_001 ]]; then
             /bin/cp ${Read1} ${Name}.1.fastq.gz && \
             /bin/cp ${Read2} ${Name}.2.fastq.gz
         else
             export PYTHONPATH=/opt/cutadapt/lib/python2.7/site-packages/ && \
             /opt/cutadapt/bin/cutadapt -a ${Adapters[0]} -A ${Adapters[1]} -o ${Name}.1.fastq.gz -p ${Name}.2.fastq.gz ${Read1} ${Read2}
         fi
     }

     runtime {
         docker_image: "registry.gsc.wustl.edu/fdu/cutadapt:1"
         cpu: "1"
         memory: "8 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File read1 = "${Name}.1.fastq.gz"
         File read2 = "${Name}.2.fastq.gz"
     }
}

task dragen_align {
     String Name
     String DragenRef
     String fastq1
     String fastq2
     String RG
     String SM
     String LB
     String AmpliconBed
     String CoverageBed
     Int? TrimLen
     String OutputDir
     String SubDir
     String jobGroup
     String queue

     String batch = basename(OutputDir)

     String StagingDir = "/staging/runs/Haloplex/"
     String LocalAlignDir = StagingDir + "align/" + batch
     String LocalSampleDir = LocalAlignDir + "/" + SubDir
     String log = StagingDir + "log/" + Name + "_align.log"

     String outdir = OutputDir + "/" + SubDir
     String dragen_outdir = outdir + "/dragen"

     command {
         if [ ! -d "${LocalAlignDir}" ]; then
             /bin/mkdir ${LocalAlignDir}
         fi

         /bin/mkdir ${LocalSampleDir} && \
         /bin/mkdir ${outdir} && \
         /opt/edico/bin/dragen -r ${DragenRef} --tumor-fastq1 ${fastq1} --tumor-fastq2 ${fastq2} --RGSM-tumor ${SM} --RGID-tumor ${RG} --RGLB-tumor ${LB} --enable-map-align true --enable-sort true --enable-map-align-output true --vc-enable-umi-liquid true --gc-metrics-enable=true --qc-coverage-region-1 ${CoverageBed} --qc-coverage-reports-1 full_res --umi-enable true --umi-library-type=random-simplex --umi-min-supporting-reads 1 --enable-variant-caller=true --vc-target-bed ${CoverageBed} --umi-metrics-interval-file ${CoverageBed} --read-trimmers=fixed-len --trim-r1-5prime=${default=1 TrimLen} --trim-r1-3prime=${default=1 TrimLen} --trim-r2-5prime=${default=1 TrimLen} --trim-r2-3prime=${default=1 TrimLen} --output-dir ${LocalSampleDir} --output-file-prefix ${Name} --output-format BAM &> ${log} && \
         /bin/mv ${log} ./ && \
         /bin/mv ${LocalSampleDir} ${dragen_outdir}
     }

     runtime {
         docker_image: "seqfu/centos7-dragen-3.9.3:latest"
         cpu: "20"
         memory: "200 G"
         queue: queue
#        host: "compute1-dragen-3"
         job_group: jobGroup 
     }

     output {
         File bam = "${dragen_outdir}/${Name}_tumor.bam"
         File bai = "${dragen_outdir}/${Name}_tumor.bam.bai"
     }
}

task convert_bam {
     String Bam
     String BamIndex
     String Name
     String refFasta
     String AmpliconBed
     String OutputDir
     String SubDir
     String jobGroup
     String queue

     String outdir = OutputDir + "/" + SubDir

     command <<<
         /usr/local/bin/tagbam -v ${Bam} ${AmpliconBed} /tmp/tagged.bam > ${Name}.ampinfo.txt && \
         /usr/local/bin/samtools view -T ${refFasta} -C -o "${Name}.cram" /tmp/tagged.bam && \
         /usr/local/bin/samtools index "${Name}.cram" &&
         (cut -f 5 ${Name}.ampinfo.txt && cut -f 4 ${AmpliconBed}) | sort | uniq -c | awk '!/\./ { print $2,$1-1; }' > ${Name}.ampcounts.txt && \
         /bin/cp ${Name}.ampinfo.txt ${outdir} && \
         /bin/cp ${Name}.ampcounts.txt ${outdir} && \
         /bin/cp ${Name}.cram ${outdir} && \
         /bin/cp ${Name}.cram.crai ${outdir}
     >>>

     runtime {
         docker_image: "registry.gsc.wustl.edu/mgi-cle/myeloseqhd:v1"
         cpu: "1"
         memory: "24 G"
         queue: queue
         job_group: jobGroup 
     }

     output {
         File cram  = "${Name}.cram"
         File crai = "${Name}.cram.crai"
         File info = "${Name}.ampinfo.txt"
         File counts = "${Name}.ampcounts.txt"
     }
}

task move_demux_fastq {
     Array[String] order_by
     String Batch
     String DemuxFastqDir
     String queue
     String jobGroup

     String LocalDemuxFastqDir = "/staging/runs/Haloplex/demux_fastq/" + Batch

     command {
         if [ -d "${LocalDemuxFastqDir}" ]; then
             /bin/mv ${LocalDemuxFastqDir} ${DemuxFastqDir}
         fi
     }
     runtime {
         docker_image: "ubuntu:xenial"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}

task batch_qc {
     Array[String] order_by
     String BatchDir
     String QC_pl
     String queue
     String jobGroup

     command {
         /usr/bin/perl ${QC_pl} ${BatchDir}
     }
     runtime {
         docker_image: "registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:compute1-20"
         memory: "4 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}

task remove_rundir {
     Array[String] order_by
     String rundir
     String queue
     String jobGroup

     command {
         if [ -d "${rundir}" ]; then
             /bin/rm -Rf ${rundir}
         fi
     }
     runtime {
         docker_image: "ubuntu:xenial"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}
