version 1.0

import "GatewayseqAnalysis.wdl" as subWF

workflow Gatewayseq {
    input {
        File SampleSheet
        # sample sheet has this structure:
        # index  name  RG_ID  RG_FLOWCELL  RG_LANE  RG_LIB  RG_SAMPLE [R1] [R2]

        File? DemuxSampleSheet
        String? IlluminaDir
        String? DragenEnv

        String CNVNormFile
        String CNVSegBed
        String CoverageBed
        String GeneCoverageBed
        String OtherCoverageBed
        String SVGeneList
        String CivicCachePath
        String JobGroup
        String OutputDir
        String Queue
        String DragenQueue
        String DragenDockerImage

        Int readfamilysize
        Int CNVfilterlength
    }

    String DragenReference = "/staging/runs/Chromoseq/refdata/dragen_hg38"
    String Reference       = "/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/refdata/hg38/all_sequences.fa"
    String ReferenceDict   = "/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/refdata/hg38/all_sequences.dict"

    String VEP          = "/storage1/fs1/gtac-mgi/Active/CLE/reference/VEP_cache"
    String NirvanaDB    = "/storage1/fs1/gtac-mgi/Active/CLE/reference/dragen_align_inputs/hg38/nirvana_annotation_data"
    String QcMetrics    = "/storage1/fs1/duncavagee/Active/SEQ/GatewaySeq/process/git/cle-gatewayseq/accessory_files/GatewaySeqQCMetrics.json"
    String HaplotectBed = "/storage1/fs1/duncavagee/Active/SEQ/GatewaySeq/process/git/cle-gatewayseq/accessory_files/myeloseq.haplotect_snppairs_hg38.041718.bed"

    String SvNoiseFile = "/staging/runs/Chromoseq/dragen_align_inputs/hg38/WGS_v1.0.0_hg38_sv_systematic_noise.bedpe.gz"

    String MSIMicroSatFile = "/storage1/fs1/duncavagee/Active/SEQ/GatewaySeq/process/git/cle-gatewayseq/accessory_files/microsatellite_file.hg38-1.gws"
    String MSIRefNormalDir = "/storage1/fs1/duncavagee/Active/SEQ/GatewaySeq/process/test/MSI/normal_ref_dir"

    String QC_pl = "/storage1/fs1/duncavagee/Active/SEQ/GatewaySeq/process/git/cle-gatewayseq/scripts/QC_metrics.pl"
    String DemuxFastqDir = "/storage1/fs1/gtac-mgi/Active/CLE/assay/gatewayseq/demux_fastq"


    call update_local_civic_cache {
        input: CivicCachePath=CivicCachePath,
               queue=Queue,
               jobGroup=JobGroup
    }

    if (defined(DemuxSampleSheet)){
        call dragen_demux {
            input: Dir=IlluminaDir,
                   OutputDir=OutputDir,
                   SampleSheet=DemuxSampleSheet,
                   DragenEnv=DragenEnv,
                   DragenDockerImage=DragenDockerImage,
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
        call dragen_align {
            input: DragenRef=DragenReference,
                   fastq1=samples[7],
                   fastq2=samples[8],
                   Name=samples[1],
                   RG=samples[3] + '.' + samples[4] + '.' + samples[0],
                   SM=samples[6],
                   LB=samples[5] + '.' + samples[0],
                   readfamilysize=readfamilysize,
                   CNVfilterlength=CNVfilterlength,
                   CNVNormFile=CNVNormFile,
                   CNVSegBed=CNVSegBed,
                   NirvanaDB=NirvanaDB,
                   CoverageBed=CoverageBed,
                   GeneCoverageBed=GeneCoverageBed,
                   OtherCoverageBed=OtherCoverageBed,
                   SvNoiseFile=SvNoiseFile,
                   MSIMicroSatFile=MSIMicroSatFile,
                   MSIRefNormalDir=MSIRefNormalDir,
                   OutputDir=OutputDir,
                   SubDir=samples[1] + '_' + samples[0],
                   DragenEnv=DragenEnv,
                   DragenDockerImage=DragenDockerImage,
                   queue=DragenQueue,
                   jobGroup=JobGroup
        }

        call subWF.GatewayseqAnalysis {
            input: Name=samples[1],
                   refFasta=Reference,
                   ReferenceDict=ReferenceDict,
                   HaplotectBed=HaplotectBed,
                   CivicCachePath=CivicCachePath,
                   Vepcache=VEP,
                   CoverageBed=CoverageBed,
                   SVGeneList=SVGeneList,
                   QcMetrics=QcMetrics,
                   DragenVcf=dragen_align.vcf,
                   DragenVcfIndex=dragen_align.index,
                   DragenSvVcf=dragen_align.svvcf,
                   DragenSvVcfIndex=dragen_align.svindex,
                   DragenCram=dragen_align.cram,
                   DragenCramIndex=dragen_align.crai,
                   OutputDir=OutputDir,
                   SubDir=samples[1] + '_' + samples[0],
                   Queue=Queue,
                   JobGroup=JobGroup
        }


    }

    call batch_qc {
        input: order_by=GatewayseqAnalysis.all_done,
               BatchDir=OutputDir,
               QC_pl=QC_pl,
               queue=Queue,
               jobGroup=JobGroup
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
}


task dragen_demux {
     input {
         String? Dir
         String OutputDir
         File? SampleSheet
         String? DragenEnv
         String DragenDockerImage
         String jobGroup
         String queue
     }

     String batch = basename(OutputDir)
     String StagingDir = "/staging/runs/GatewaySeq/"
     String LocalFastqDir = StagingDir + "demux_fastq/" + batch
     String LocalReportDir = LocalFastqDir + "/Reports"
     String LocalSampleSheet = StagingDir + "sample_sheet/" + batch + '.csv'
     String log = StagingDir + "log/" + batch + "_demux.log"
     String DemuxReportDir = OutputDir + "/dragen_demux_reports"

     command <<<
         /bin/cp ~{SampleSheet} ~{LocalSampleSheet} && \
         /opt/edico/bin/dragen --bcl-conversion-only true --bcl-only-matched-reads true --strict-mode true --sample-sheet ~{LocalSampleSheet} --bcl-input-directory ~{Dir} --output-directory ~{LocalFastqDir} &> ~{log} && \
         /bin/ls ~{LocalFastqDir}/*_R1_001.fastq.gz > Read1_list.txt && \
         /bin/ls ~{LocalFastqDir}/*_R2_001.fastq.gz > Read2_list.txt && \
         /bin/mv ~{log} ./ && \
         /bin/rm -f ~{LocalSampleSheet} && \
         /bin/cp -r ~{LocalReportDir} ~{DemuxReportDir}
     >>>

     runtime {
         docker_image: DragenDockerImage
         dragen_env: DragenEnv
         cpu: "20"
         memory: "200 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File read1 = "Read1_list.txt"
         File read2 = "Read2_list.txt"
     }
}

task prepare_samples {
     input {
         File SampleSheet
         String Fastq1
         String Fastq2
         String jobGroup
         String queue
     }

     command <<<
         /bin/cp ~{Fastq1} 1.tmp.txt
         /bin/cp ~{Fastq2} 2.tmp.txt
         /usr/bin/perl -e 'open(R1,"1.tmp.txt"); @r1 = <R1>; \
             chomp @r1; close R1;\
             open(R2,"2.tmp.txt"); @r2 = <R2>; \
             chomp @r2; close R2; \
             open(SS,"~{SampleSheet}");
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
         docker_image: "docker1(registry.gsc.wustl.edu/genome/lims-compute-xenial:1)"
         cpu: "1"
         memory: "4 G"
         queue: queue
         job_group: jobGroup
     }

     output {
         File sample_sheet = "sample_sheet.txt"
     }
}

task dragen_align {
     input {
         String Name
         String DragenRef
         String fastq1
         String fastq2
         String RG
         String SM
         String LB
         String CNVfilterlength
         String CNVNormFile
         String CNVSegBed
         String NirvanaDB
         String CoverageBed
         String GeneCoverageBed
         String OtherCoverageBed
         String SvNoiseFile
         String MSIMicroSatFile
         String MSIRefNormalDir
         String OutputDir
         String SubDir
         String jobGroup
         String queue
         String DragenDockerImage
         String? DragenEnv

         Int? TrimLen
         Int readfamilysize
     }

     String batch = basename(OutputDir)
     String StagingDir = "/staging/runs/GatewaySeq/"
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
         /opt/edico/bin/dragen -r ${DragenRef} --tumor-fastq1 ${fastq1} --tumor-fastq2 ${fastq2} --RGSM-tumor ${SM} --RGID-tumor ${RG} --RGLB-tumor ${LB} \
         --umi-enable true --umi-library-type=random-simplex --umi-min-supporting-reads ${readfamilysize} --umi-metrics-interval-file ${CoverageBed} \
         --enable-map-align true --enable-sort true --enable-map-align-output true --gc-metrics-enable=true \
         --qc-coverage-ignore-overlaps=true --qc-coverage-region-1 ${CoverageBed} --qc-coverage-reports-1 full_res --qc-coverage-region-2 ${OtherCoverageBed} \
         --enable-variant-caller=true --vc-target-bed ${GeneCoverageBed} --vc-enable-umi-solid true --vc-combine-phased-variants-distance 3 \
         --vc-enable-orientation-bias-filter true --vc-enable-triallelic-filter false \
         --enable-sv true --sv-exome true --sv-output-contigs true --sv-systematic-noise ${SvNoiseFile} --sv-hyper-sensitivity true \
         --sv-min-edge-observations 2 --sv-min-candidate-spanning-count 1 --sv-use-overlap-pair-evidence true \
         --enable-cnv true --cnv-target-bed ${GeneCoverageBed} --cnv-enable-ref-calls false --cnv-filter-length ${CNVfilterlength} \
         --cnv-normals-list ${CNVNormFile} --cnv-segmentation-mode bed --cnv-segmentation-bed ${CNVSegBed} \
         --enable-tmb true --qc-coverage-region-3 ${CoverageBed} --qc-coverage-tag-3=tmb --qc-coverage-reports-3=callability \
         --enable-variant-annotation true --variant-annotation-assembly GRCh38 --variant-annotation-data ${NirvanaDB} \
         --msi-command tumor-only --msi-coverage-threshold 60 --msi-microsatellites-file ${MSIMicroSatFile} --msi-ref-normal-dir ${MSIRefNormalDir} \
         --output-dir ${LocalSampleDir} --output-file-prefix ${Name} --output-format CRAM &> ${log} && \
         /bin/mv ${log} ./ && \
         /bin/mv ${LocalSampleDir} ${dragen_outdir}
     }

     runtime {
         docker_image: DragenDockerImage
         dragen_env: DragenEnv
         cpu: "20"
         memory: "200 G"
         queue: queue
         job_group: jobGroup
     }

     output {
         File cram = "${dragen_outdir}/${Name}_tumor.cram"
         File crai = "${dragen_outdir}/${Name}_tumor.cram.crai"
         File vcf = "${dragen_outdir}/${Name}.hard-filtered.vcf.gz"
         File index = "${dragen_outdir}/${Name}.hard-filtered.vcf.gz.tbi"
         File svvcf = "${dragen_outdir}/${Name}.sv.vcf.gz"
         File svindex = "${dragen_outdir}/${Name}.sv.vcf.gz.tbi"
     }
}

task move_demux_fastq {
     input {
         Array[String] order_by
         String Batch
         String DemuxFastqDir
         String queue
         String jobGroup
     }

     String LocalDemuxFastqDir = "/staging/runs/GatewaySeq/demux_fastq/" + Batch

     command {
         if [ -d "${LocalDemuxFastqDir}" ]; then
             /bin/mv ${LocalDemuxFastqDir} ${DemuxFastqDir}
         fi
     }
     runtime {
         docker_image: "docker1(ubuntu:xenial)"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}

task remove_rundir {
     input {
         Array[String] order_by
         String? rundir
         String queue
         String jobGroup
     }
     command {
         if [ -d "${rundir}" ]; then
             /bin/rm -Rf ${rundir}
         fi
     }
     runtime {
         docker_image: "docker1(ubuntu:xenial)"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}

task update_local_civic_cache {
     input {
         String CivicCachePath
         String jobGroup
         String queue
     }
     command {
         /usr/local/bin/civicpy update --hard --cache-save-path ${CivicCachePath}
     }
     runtime {
         docker_image: "docker1(griffithlab/civicpy:3.0.0)"
         cpu: "1"
         memory: "8 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}

task batch_qc {
     input {
         Array[String] order_by
         String BatchDir
         String QC_pl
         String queue
         String jobGroup
     }

     command {
         if [ -n "$(/bin/ls -d ${BatchDir}/TW*)" ]; then
             /bin/chmod -R 777 ${BatchDir}/TW*
         fi

         /usr/bin/perl ${QC_pl} ${BatchDir}
     }
     runtime {
         docker_image: "docker1(registry.gsc.wustl.edu/mgi-cle/myeloseqhd:v2)"
         memory: "4 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}
