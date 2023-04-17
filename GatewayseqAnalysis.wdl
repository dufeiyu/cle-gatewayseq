version 1.0

workflow GatewayseqAnalysis {
    input {
        String Name
        String DragenCram
        String DragenCramIndex
        String DragenVcf
        String DragenVcfIndex
        String DragenSVVCF
        
        String refFasta 
        String ReferenceDict
        String Vepcache        

        String QcMetrics
        String CoverageBed
        String HaplotectBed
        String CivicCachePath

        String SubDir
        String OutputDir

        String Queue
        String JobGroup 
    }

    Int MinSVReads = 10
    
    call run_civic {
        input: CivicCachePath=CivicCachePath,
               Vcf=DragenVcf,
               VcfIndex=DragenVcfIndex,
               Name=Name,
               queue=Queue,
               jobGroup=JobGroup
    }

    call run_vep {
        input: Vcf=run_civic.vcf,
               refFasta=refFasta,
               Vepcache=Vepcache,
               Name=Name,
               queue=Queue,
               jobGroup=JobGroup
    }

    call filter_sv {
        input: Vcf=DragenSVVCF,
               Bed=CoverageBed,
               Name=Name,
               MinReads=MinSVReads,
               queue=Queue,
               jobGroup=JobGroup
    }

    call annotate_sv {
        input: Vcf=filter_sv.vcf,
               refFasta=refFasta,
               Vepcache=Vepcache,
               Name=Name,
               queue=Queue,
               jobGroup=JobGroup
    }


    call run_haplotect {
        input: refFasta=refFasta,
               refDict=ReferenceDict,
               Cram=DragenCram,
               CramIndex=DragenCramIndex,
               Bed=HaplotectBed,
               Name=Name,
               queue=Queue,
               jobGroup=JobGroup
    }

    call gather_files {
        input: OutputFiles=[run_haplotect.out_file,
               run_haplotect.sites_file,
               run_vep.vcf,run_vep.index,
               annotate_sv.vcf,
               annotate_sv.index],
               OutputDir=OutputDir,
               SubDir=SubDir,
               queue=Queue,
               jobGroup=JobGroup
    }

    call make_report {
        input: order_by=gather_files.done,
               Name=Name,
               CoverageBed=CoverageBed,
               QcMetrics=QcMetrics,
               OutputDir=OutputDir,
               SubDir=SubDir,
               queue=Queue,
               jobGroup=JobGroup
    }

    output {
        String all_done = make_report.done
    }
}


task run_civic {
     input {
         File Vcf
         File VcfIndex
         String CivicCachePath
         String Name
         String jobGroup
         String queue
     }
     command {
         export CIVICPY_CACHE_FILE=${CivicCachePath} && \
         /usr/local/bin/civicpy annotate-vcf --input-vcf ${Vcf}  --output-vcf ${Name}.civic.vcf.gz --reference GRCh38 -i accepted
     }
     runtime {
         docker_image: "docker1(griffithlab/civicpy:3.0.0)"
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
         docker_image: "docker1(registry.gsc.wustl.edu/mgi-cle/vep105-htslib1.9:1)"
         cpu: "1"
         memory: "10 G"
         queue: queue
         job_group: jobGroup
     }
     output {
        File vcf = "${Name}.annotated.vcf.gz"
        File index = "${Name}.annotated.vcf.gz.tbi"

     }
}

task filter_sv {
    input { 
            
        String Vcf
        String Name
        String Bed
        Int MinReads
    
        String jobGroup
        String queue
    
    }

    command {
        grep Fusion ${Bed} > sv_region.bed && \
        bcftools view -i 'SVTYPE=="BND" && FMT/SR[0:1]>=${MinReads}' -R sv_region.bed -Oz -o ${Name}.sv.filtered.vcf.gz $Vcf && \
        tabix -p vcf $Name".sv.filtered.vcf.gz"
    }

    runtime {
         docker_image: "docker1(ghcr.io/dhslab/docker-baseimage:latest)"
         cpu: "1"
         memory: "10 G"
         queue: queue
         job_group: jobGroup
     }

     output {
        File vcf = "${Name}.sv.filtered.vcf.gz"
        File index = "${Name}.sv.filtered.vcf.gz.tbi"

     }
}

task annotate_sv {
    input { 
            
        String Vcf
        String Name
    
        String refFasta
        String Vepcache
    
        String jobGroup
        String queue
    
    }

    command {
        /usr/bin/perl -I /opt/vep/lib/perl/VEP/Plugins /opt/vep/src/ensembl-vep/vep --format vcf --vcf --fasta ${refFasta} \
        --per_gene --symbol --term SO -o ${Name}.sv.annotated.vcf -i ${Vcf} --offline --cache --dir ${Vepcache} && \
        /usr/local/bin/bgzip ${Name}.sv.annotated.vcf && /usr/local/bin/tabix -p vcf ${Name}.sv.annotated.vcf.gz
    }

    runtime {
         docker_image: "docker1(registry.gsc.wustl.edu/mgi-cle/vep105-htslib1.9:1)"
         cpu: "1"
         memory: "10 G"
         queue: queue
         job_group: jobGroup
     }

     output {
        File vcf = "${Name}.sv.annotated.vcf.gz"
        File index = "${Name}.sv.annotated.vcf.gz.tbi"

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
         docker_image: "docker1(registry.gsc.wustl.edu/mgi-cle/haplotect:0.3)"
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
         docker_image: "docker1(registry.gsc.wustl.edu/genome/lims-compute-xenial:1)"
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
         String OutputDir
         String SubDir
         String queue
         String jobGroup
     }

     String SampleOutDir = OutputDir + "/" + SubDir

     command {
         /usr/bin/python3 /storage1/fs1/duncavagee/Active/SEQ/GatewaySeq/process/git/cle-gatewayseq/make_gwseq_report.py -n ${Name} -d ${SampleOutDir} -q ${QcMetrics} && \
         /bin/mv ./*.report.txt ./*.report.json ${SampleOutDir}
     }
     runtime {
         docker_image: "docker1(registry.gsc.wustl.edu/mgi-cle/myeloseqhd:v2)"
         cpu: "1"
         memory: "16 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}

