
workflow BwaAndBedWorkflow{

    Array[File] r1_files
    Array[File] r2_files
    File bwa_fa
    File bwa_amb
    File bwa_ann
    File bwa_bwt
    File bwa_fai
    File bwa_pac
    File bwa_sa
    File bwa_dict

    String git_repo_url
    String git_commit_hash

    Array[Pair[File, File]] fastq_pairs = zip(r1_files, r2_files)

    scatter(item in fastq_pairs){

        call trim_reads {
            input:
                r1 = item.left,
                r2 = item.right
        } 

        call align_reads {
            input:
                r1 = trim_reads.trimmed_r1,
                r2 = trim_reads.trimmed_r2,
                bwa_fa = bwa_fa,
                bwa_amb = bwa_amb,
                bwa_ann = bwa_ann,
                bwa_bwt = bwa_bwt,
                bwa_fai = bwa_fai,
                bwa_pac = bwa_pac,
                bwa_sa = bwa_sa,
                bwa_dict = bwa_dict
        }   

        call run_coverage_and_flagstat {
            input:
                bam = align_reads.bam,
                bai = align_reads.bai
        }
    }



    output {
        Array[File] trimmed_r1 = trim_reads.trimmed_r1
        Array[File] trimmed_r2 = trim_reads.trimmed_r2
        Array[File] bams = align_reads.bam
        Array[File] bais = align_reads.bai
        Array[File] beds = run_coverage_and_flagstat.bed
        Array[File] flagstats = run_coverage_and_flagstat.flagstat
    }

    meta {
        workflow_title : "BWA align and coverage"
        workflow_short_description : "Run BWA and bedtools coverage on FASTQ."
        workflow_long_description : "Run BWA and bedtools coverage on FASTQ."
    }

}

task trim_reads {

    File r1
    File r2

    # Extract the samplename from the fastq filename
    String sample_name = basename(r1, "_R1.fastq.gz")

    Int disk_size = 200


    command {
        java -jar /opt/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
            -trimlog ${sample_name}.trim.log \
            -summary ${sample_name}.trim_summary.log \
            ${r1} ${r2} \
            -baseout ${sample_name}.trimmed.fastq.gz \
            ILLUMINACLIP:/opt/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:8:true
    }

    output {
        File trimmed_r1 = "${sample_name}.trimmed_1P.fastq.gz"
        File trimmed_r2 = "${sample_name}.trimmed_2P.fastq.gz"
    }

    runtime {
        docker: "docker.io/blawney/msi_container:v2"
        cpu: 4
        memory: "12 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task align_reads {

    File r1
    File r2
    File bwa_fa
    File bwa_amb
    File bwa_ann
    File bwa_bwt
    File bwa_fai
    File bwa_pac
    File bwa_sa
    File bwa_dict

    Int disk_size = 200

    String sample_name = basename(r1, ".trimmed_1P.fastq.gz")

    command {

        /opt/software/bwa-0.7.17/bwa mem -t 4 ${bwa_fa} ${r1} ${r2} \
        | /opt/software/samtools-1.10/samtools view -bht ${bwa_fa} - \
        | /opt/software/samtools-1.10/samtools sort -o MM23-CT1.bam -;
        /opt/software/samtools-1.10/samtools index ${sample_name}.bam
    }

    output {
        File bam = "${sample_name}.bam"
        File bai = "${sample_name}.bam.bai"
    }

    runtime {
        docker: "docker.io/blawney/msi_container:v2"
        cpu: 8
        memory: "12 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task run_coverage_and_flagstat {

    File bam
    File bai
    String sample_name = basename(bam, ".bam")
    Int disk_size = 100

    command {
        /opt/software/bedtools2/bin/bedtools genomecov -ibam ${bam} -bga > ${sample_name}.bed
        /opt/software/samtools-1.10/samtools flagstat ${bam} > ${sample_name}.flagstat
    }

    output {
        File bed = "${sample_name}.bed" 
        File flagstat = "${sample_name}.flagstat" 
    }

    runtime {
        docker: "docker.io/blawney/msi_container:v2"
        cpu: 4
        memory: "12 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
