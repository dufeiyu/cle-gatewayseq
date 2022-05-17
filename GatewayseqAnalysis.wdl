workflow GatewayseqAnalysis {

    String Cram
    String CramIndex
    String DragenVcf
    String DragenVcfIndex
    String Name

    String refFasta 
    String ReferenceDict
    String Vepcache

    String CoverageBed
#    String HaplotectBed

#    String QcMetrics

#    String mrn
#    String accession
#    String DOB
#    String sex
#    String exception
#    String RunInfoString

    String SubDir
    String OutputDir

    String Queue
    String JobGroup 

    call run_vep {
        input: Vcf=DragenVcf,
	       VcfIndex=DragenVcfIndex,
               refFasta=refFasta,
               Vepcache=Vepcache,
               Name=Name,
               queue=Queue,
               jobGroup=JobGroup
    }

    call run_haplotect {
        input: refFasta=refFasta,
               refDict=ReferenceDict,
               Cram=Cram,
               CramIndex=CramIndex,
               Bed=HaplotectBed,
               Name=Name,
               queue=Queue,
               jobGroup=JobGroup
    }

     call gather_files {
        input: OutputFiles=[run_haplotect.out_file,
               run_haplotect.sites_file,
               run_vep.vcf,
               run_vep.filtered_vcf],
               OutputDir=OutputDir,
               SubDir=SubDir,
               queue=Queue,
               jobGroup=JobGroup
    }

    call make_report {
        input: order_by=gather_files.done,
#               Name=Name,
#               mrn=mrn,
#               accession=accession,
#               DOB=DOB,
#               sex=sex,
#               exception=exception,
#               RunInfoString=RunInfoString,
               CoverageBed=CoverageBed,
#               QcMetrics=QcMetrics,
               OutputDir=OutputDir,
               SubDir=SubDir,
               queue=Queue,
               jobGroup=JobGroup

    output {
        String all_done = haloplex_qc.done
    }
}


task run_vep {
     File Vcf
     File VcfIndex
     String refFasta
     String Vepcache
     Float? maxAF
     String Name
     String jobGroup
     String queue

     command {
             /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl --format vcf \
             --vcf --plugin Downstream --fasta ${refFasta} --hgvs --symbol --term SO --flag_pick \
             -i ${Vcf} --offline --cache --max_af --dir ${Vepcache} -o ${Name}.annotated.vcf && \
             /opt/htslib/bin/bgzip ${Name}.annotated.vcf && /usr/bin/tabix -p vcf ${Name}.annotated.vcf.gz
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
     String Cram
     String CramIndex
     String Bed
     String Name
     String refDict
     String refFasta
     String queue
     String jobGroup

     Int? MinReads

     command <<<
         /usr/bin/awk -v OFS="\t" '{ $2=$2-1; print; }' ${Bed} > /tmp/pos.bed && \
         /usr/local/openjdk-8/bin/java -Xmx6g \
         -jar /opt/hall-lab/gatk-package-4.1.8.1-18-ge2f02f1-SNAPSHOT-local.jar Haplotect \
         -I ${Cram} -R ${refFasta} --sequence-dictionary ${refDict} \
         -mmq 20 -mbq 20 -max-depth-per-sample 10000 -gstol 0.001 -mr ${default=10 MinReads} \
         -htp ${Bed} -L /tmp/pos.bed -outPrefix ${Name}
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
     Array[String] OutputFiles
     String OutputDir
     String? SubDir
     String jobGroup
     String queue

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
