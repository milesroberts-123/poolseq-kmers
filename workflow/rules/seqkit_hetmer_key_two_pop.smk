rule seqkit_hetmer_key_two_pop:
    input:
        tempsamplefasta="seqkit_results/samples_{ID}_p1p2.fasta",
        vcffilled="slim_results/samples_filled_{ID}_p1p2.vcf.gz",
    output:
        poskey="seqkit_results/center_kmer_pairs_{ID}_p1p2.txt",
        snppos=temp("seqkit_results/snp_positions_{ID}_p1p2.txt"),
    params:
        L=config["L"],
    conda:
        "../envs/seqkit.yaml"
    log:
        "logs/seqkit/{ID}.log",
    shell:
        """
        # get list of snp positions
        zcat {input.vcffilled} | grep -v "^#" | cut -f 2 > {output.snppos}

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
            cat {input.tempsamplefasta} | seqkit subseq -r $(echo $center_start):$(echo $center_end) | grep -v "^>" | sort -u | tr '\n' ' ' | echo $i $(cat -) >> {output.poskey}
        done
        """
