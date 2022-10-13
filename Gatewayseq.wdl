version 1.0

workflow Gatewayseq {
    input {
        File SampleSheet
        # sample sheet has this structure:
        # index  name  RG_ID  RG_FLOWCELL  RG_LANE  RG_LIB  RG_SAMPLE [R1] [R2]

        File? DemuxSampleSheet
        String? IlluminaDir

        String CNVNormFile
        String CoverageBed
        String CivicCachePath
        String JobGroup
        String OutputDir
        String Queue
        String DragenQueue
    }

    String DragenReference = "/staging/runs/Chromoseq/refdata/dragen_hg38"
    String Reference       = "/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/refdata/hg38/all_sequences.fa"
    String ReferenceDict   = "/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/refdata/hg38/all_sequences.dict"

    String VEP          = "/storage1/fs1/gtac-mgi/Active/CLE/reference/VEP_cache"
    String QcMetrics    = "/storage1/fs1/duncavagee/Active/SEQ/GatewaySeq/process/git/cle-gatewayseq/accessory_files/GatewaySeqQCMetrics.json"
    String HaplotectBed = "/storage1/fs1/duncavagee/Active/SEQ/GatewaySeq/process/git/cle-gatewayseq/accessory_files/myeloseq.haplotect_snppairs_hg38.041718.bed"

    String SvNoiseFile = "/staging/runs/Chromoseq/dragen_align_inputs/hg38/WGS_v1.0.0_hg38_sv_systematic_noise.bedpe.gz"

    String QC_pl = "/storage1/fs1/duncavagee/Active/SEQ/GatewaySeq/process/git/cle-gatewayseq/QC_metrics.pl"
    String DemuxFastqDir = "/scratch1/fs1/gtac-mgi/CLE/gatewayseq/demux_fastq"


    Int readfamilysize  = 1
    Int CNVfilterlength = 100000

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
                   CoverageBed=CoverageBed,
                   SvNoiseFile=SvNoiseFile,
                   OutputDir=OutputDir,
                   SubDir=samples[1] + '_' + samples[0],
                   queue=DragenQueue,
                   jobGroup=JobGroup
        }

        call run_civic {
            input: order_by=update_local_civic_cache.done,
                   CivicCachePath=CivicCachePath,
                   Vcf=dragen_align.vcf,
                   VcfIndex=dragen_align.index,
                   Name=samples[1],
                   queue=Queue,
                   jobGroup=JobGroup
        }

        call run_vep {
            input: Vcf=run_civic.vcf,
                   refFasta=Reference,
                   Vepcache=VEP,
                   Name=samples[1],
                   queue=Queue,
                   jobGroup=JobGroup
        }

        call run_haplotect {
            input: refFasta=Reference,
                   refDict=ReferenceDict,
                   Cram=dragen_align.cram,
                   CramIndex=dragen_align.crai,
                   Bed=HaplotectBed,
                   Name=samples[1],
                   queue=Queue,
                   jobGroup=JobGroup
        }

        call gather_files {
            input: OutputFiles=[run_haplotect.out_file,
                   run_haplotect.sites_file,
                   run_vep.vcf],
                   OutputDir=OutputDir,
                   SubDir=samples[1] + '_' + samples[0],
                   queue=Queue,
                   jobGroup=JobGroup
        }

        call make_report {
            input: order_by=gather_files.done,
                   Name=samples[1],
                   CoverageBed=CoverageBed,
                   QcMetrics=QcMetrics,
                   SampleOutDir=OutputDir + '/' + samples[1] + '_' + samples[0],
                   queue=Queue,
                   jobGroup=JobGroup
        }
    }

    if (defined(DemuxSampleSheet)){
        call move_demux_fastq {
            input: order_by=make_report.done,
                   Batch=basename(OutputDir),
                   DemuxFastqDir=DemuxFastqDir,
                   queue=DragenQueue,
                   jobGroup=JobGroup
        }

        call remove_rundir {
             input: order_by=make_report.done,
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
         docker_image: "seqfu/centos7-dragen-4.0.3:latest"
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
         String CoverageBed
         String SvNoiseFile
         String OutputDir
         String SubDir
         String jobGroup
         String queue

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
         /opt/edico/bin/dragen -r ${DragenRef} --tumor-fastq1 ${fastq1} --tumor-fastq2 ${fastq2} --RGSM-tumor ${SM} --RGID-tumor ${RG} --RGLB-tumor ${LB} --enable-map-align true --enable-sort true --enable-map-align-output true --cnv-target-bed ${CoverageBed} --vc-target-bed ${CoverageBed} --sv-call-regions-bed ${CoverageBed} --vc-enable-umi-solid true --vc-combine-phased-variants-distance 3 --vc-enable-orientation-bias-filter true --vc-enable-triallelic-filter false --enable-sv true --sv-exome true --sv-output-contigs true --sv-systematic-noise ${SvNoiseFile} --enable-cnv true --cnv-enable-ref-calls false --cnv-filter-length ${CNVfilterlength} --cnv-normals-file ${CNVNormFile} --gc-metrics-enable=true --qc-coverage-ignore-overlaps=true --qc-coverage-region-1 ${CoverageBed} --qc-coverage-reports-1 full_res --umi-enable true --umi-library-type=random-simplex --umi-min-supporting-reads ${readfamilysize} --enable-variant-caller=true --umi-metrics-interval-file ${CoverageBed} --output-dir ${LocalSampleDir} --output-file-prefix ${Name} --output-format CRAM &> ${log} && \
         /bin/mv ${log} ./ && \
         /bin/mv ${LocalSampleDir} ${dragen_outdir}
     }

     runtime {
         docker_image: "seqfu/centos7-dragen-4.0.3:latest"
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
         docker_image: "ubuntu:xenial"
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
         docker_image: "ubuntu:xenial"
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
         docker_image: "griffithlab/civicpy:1.1.1"
         cpu: "1"
         memory: "8 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}

task run_civic {
     input {
         File Vcf
         File VcfIndex
         String CivicCachePath
         String order_by
         String Name
         String jobGroup
         String queue
     }
     command {
         export CIVICPY_CACHE_FILE=${CivicCachePath} && \
         /usr/local/bin/civicpy annotate-vcf --input-vcf ${Vcf}  --output-vcf ${Name}.civic.vcf.gz --reference GRCh38 -i accepted
     }
     runtime {
         docker_image: "griffithlab/civicpy:1.1.1"
         cpu: "1"
         memory: "10 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File vcf = "${Name}.civic.vcf.gz"
     }
}

task run_vep {
     input {
         File Vcf
         String refFasta
         String Vepcache
         Float? maxAF
         String Name
         String jobGroup
         String queue
     }
     command {
         /usr/local/bin/tabix -p vcf ${Vcf} && \
         /usr/bin/perl -I /opt/vep/lib/perl/VEP/Plugins /opt/vep/src/ensembl-vep/vep --format vcf \
         --vcf --plugin Downstream --fasta ${refFasta} --hgvs --symbol --term SO --flag_pick \
         -i ${Vcf} --offline --cache --max_af --dir ${Vepcache} -o ${Name}.annotated.vcf && \
         /usr/local/bin/bgzip ${Name}.annotated.vcf && /usr/local/bin/tabix -p vcf ${Name}.annotated.vcf.gz
     }
     runtime {
         docker_image: "registry.gsc.wustl.edu/mgi-cle/vep105-htslib1.9:1"
         cpu: "1"
         memory: "10 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File vcf = "${Name}.annotated.vcf.gz"
     }
}

task run_haplotect {
     input {
         String Cram
         String CramIndex
         String Bed
         String Name
         String refDict
         String refFasta
         String queue
         String jobGroup

         Int? MinReads
     }
     command <<<
         /usr/bin/awk -v OFS="\t" '{ $2=$2-1; print; }' ~{Bed} > /tmp/pos.bed && \
         /usr/local/openjdk-8/bin/java -Xmx6g \
         -jar /opt/hall-lab/gatk-package-4.1.8.1-18-ge2f02f1-SNAPSHOT-local.jar Haplotect \
         -I ~{Cram} -R ~{refFasta} --sequence-dictionary ~{refDict} \
         -mmq 20 -mbq 20 -max-depth-per-sample 10000 -gstol 0.001 -mr ~{default=10 MinReads} \
         -htp ~{Bed} -L /tmp/pos.bed -outPrefix ~{Name}
     >>>

     runtime {
         docker_image: "registry.gsc.wustl.edu/mgi-cle/haplotect:0.3"
         cpu: "1"
         memory: "8 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File out_file = "${Name}.haplotect.txt"
         File sites_file = "${Name}.haplotectloci.txt"
     }
}

task gather_files {
     input {
         Array[String] OutputFiles
         String OutputDir
         String? SubDir
         String jobGroup
         String queue
     }
     command {
         if [[ ${SubDir} != "" ]] && [[ ! -e ${OutputDir}/${SubDir} ]]; then
             mkdir ${OutputDir}/${SubDir}
         fi
         /bin/mv -f -t ${OutputDir}/${SubDir} ${sep=" " OutputFiles}
     }
     runtime {
         docker_image: "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}

task make_report {
     input {
         String order_by
         String Name
         String CoverageBed
         String QcMetrics
         String SampleOutDir
         String queue
         String jobGroup
     }
     command {
         /usr/bin/python3 /storage1/fs1/duncavagee/Active/SEQ/GatewaySeq/process/git/cle-gatewayseq/make_gwseq_report.py -n ${Name} -d ${SampleOutDir} -q ${QcMetrics} && \
         /bin/mv ./*.report.txt ./*.report.json ${SampleOutDir}
     }
     runtime {
         docker_image: "registry.gsc.wustl.edu/mgi-cle/myeloseqhd:v2"
         cpu: "1"
         memory: "16 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}
