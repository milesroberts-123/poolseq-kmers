rule seqkit_hetmer_key:
    input:
        vcf="slim_results/samples_{SID}.vcf.gz",
        fasta="seqkit_results/samples_across_pop_{SID}.fasta",
    output:
        poskey="seqkit_results/center_kmer_pairs_{SID}.txt",
        snppos=temp("seqkit_results/snp_positions_{SID}.txt"),
    params:
        L=config["L"],
    conda:
        "../envs/seqkit.yaml"
    log:
        "logs/seqkit_hetmer_key/{SID}.log",
    shell:
        """        
        # get list of snp positions
        zcat {input.vcf} | grep -v "^#" | cut -f 2 > {output.snppos}

        # load list
        readarray posarray < {output.snppos}

        # if output already exists (i.e. previous failed run) then delete it
        if [ -e "{output.poskey}" ]; then
            rm {output.poskey}
        fi
        
        # loop over snp positions
        for i in "${{posarray[@]}}"
        do
            echo Extracting k-mers for snp $i &>> {log}

            # calculate bounds of k-mers
            center_start=$(($i-15))
            center_end=$(($i+15))

            # check if bounds of k-mer extend beyond edges of chromosome
            if [ "$center_start" -lt 1 ]; then
                continue
            fi

            if [ "$center_end" -gt {params.L} ]; then
                continue
            fi
            
            # extract k-mers
            cat {input.fasta} | seqkit subseq -r $(echo $center_start):$(echo $center_end) | grep -v "^>" | sort -u | tr '\n' ' ' | echo $i $(cat -) >> {output.poskey}
        done
        """
