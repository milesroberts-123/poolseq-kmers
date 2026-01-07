rule kmer_distances:
    input:
        "cbf_table_{SID}.txt",
    output:
        "distances_{SID}.txt",
    conda:
        "../envs/cbf.yaml"
    shell:
        "python scripts/kmer_distances.py --input {input} --output {output}"
