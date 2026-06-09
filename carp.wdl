version 1.0

workflow Carp {
    input {
        String family
        String pipelineProject
        String folderPath
        String mode
        String location
        String indexSample
        String? samples
        String? genotypes
        String? options
        String outPrefix = ""
        File cytobands = "dx://project-G729Kkj4fq4Q9X1BPy0807bK:file-GGbXZvQ4fq4pKZZ54kQ3QBKB"
        File iscaRegions = "dx://project-G729Kkj4fq4Q9X1BPy0807bK:file-J6PKpK84fq4gqYJV4ff3XKFz"
        String dockerCarp = "dx://project-G729Kkj4fq4Q9X1BPy0807bK:file-J8Kf8Xj4fq4zfky9V8jX70Qq"
    }

    String defaultPrefix = if outPrefix == "" then family else outPrefix

    call findCarpFiles {
        input:
            family=family,
            pipelineProject=pipelineProject,
            folderPath=folderPath
    }

    call runCarp {
        input:
            arrayCnvs=findCarpFiles.arrCnvs,
            arrayDataFiles=findCarpFiles.arrData,
            familyVcf=findCarpFiles.vcf,
            indexSample=indexSample,
            outPrefix=defaultPrefix,
            mode=mode,
            location=location,
            samples=samples,
            genotypes=genotypes,
            options=options,
            cytobands=cytobands,
            iscaRegions=iscaRegions,
            dockerCarp=dockerCarp
    }

    output {
        Array[File] carp_pdf = runCarp.pdf
        Array[File] carp_png = runCarp.png
    }
}

task findCarpFiles {
    input {
        String family
        String pipelineProject
        String folderPath
    }

    String validFolderPath = sub(folderPath, "/?$", "/")

    command <<<
        set -euo pipefail

        mkdir arrCnvs
        mkdir arrData
        mkdir vcf

        for f in $(dx find data \
            --path ~{pipelineProject}:~{validFolderPath} \
            --name "cnvs_*.20000" \
            --brief); do \
            dx download $f -o arrCnvs/; \
        done

        for f in $(dx find data \
            --path ~{pipelineProject}:~{validFolderPath} \
            --name "*coverageBinner6.20000.data" \
            --brief); do \
            dx download $f -o arrData/; \
        done            

        for f in $(dx find data \
            --path ~{pipelineProject}:~{validFolderPath}vcf \
            --name "*.vcf.gz" \
            --brief); do \
            dx download $f -o vcf/; \
        done

        for f in $(dx find data \
            --path ~{pipelineProject}:~{validFolderPath}vcf \
            --name "*.vcf.gz.tbi" \
            --brief); do \
            dx download $f -o vcf/; \
        done
    >>>

    output {
        Array[File] arrCnvs = glob("arrCnvs/*")
        Array[File] arrData = glob("arrData/*")
        Array[File] vcf = glob("vcf/*")
    }
}

task runCarp {
    input {
        Array[File] arrayCnvs
        Array[File] arrayDataFiles
        Array[File] familyVcf
        String indexSample
        String outPrefix
        String mode
        String location
        String? samples
        String? genotypes
        String? options
        File cytobands
        File iscaRegions
        String dockerCarp
    }

    String cytobandsFile = basename(cytobands)
    String iscaFile = basename(iscaRegions)

    Int disk_gb = ceil(2*size(arrayCnvs, "GiB") + size(arrayDataFiles, "GiB") + size(familyVcf, "GiB")) + 5
    String samplesFlag = if defined(samples) then "--samples '" + select_first([samples]) + "'" else ""
    String genotypesFlag = if defined(genotypes) then "--genotypes '" + select_first([genotypes]) + "'" else ""
    String optionsFlag = if defined(options) then  select_first([options]) else ""

    command <<<
        cd /usr/carp
        inputDir=/usr/carp/input
        mkdir -p ${inputDir}

        cp ~{sep=' ' arrayCnvs} "${inputDir}/."
        cp ~{sep=' ' arrayDataFiles} "${inputDir}/."
        cp ~{sep=' ' familyVcf} "${inputDir}/."
        cp ~{cytobands} .
        cp ~{iscaRegions} .

        python3 -m src.carp \
            --inDir "${inputDir}" \
            --outDir /home/dnanexus/work/out \
            --mode ~{mode} \
            --prefix ~{outPrefix} \
            --location ~{location} \
            --proband_id ~{indexSample} ~{samplesFlag} ~{genotypesFlag} \
            --cyto ~{cytobandsFile} \
            --isca ~{iscaFile} \
            ~{optionsFlag}
    >>>

    output {
        Array[File] pdf = glob("/home/dnanexus/work/out/*.pdf")
        Array[File] png = glob("/home/dnanexus/work/out/*.png")
    }

    runtime {
        docker: "${dockerCarp}"
        gpu: false
        memory: "60 GB"
        cpu: 4
        disks: "local-disk ${disk_gb} SSD"
    }
}