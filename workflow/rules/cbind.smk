# https://unix.stackexchange.com/questions/16443/combine-text-files-column-wise

rule cbind:
    input:
        expand("cbf_results/{{SID}}_{PID}.txt", PID=[0,1])
    output:
        temp("cbf_table_{SID}.txt")
    shell:
        r"""
        paste -d' ' {input} > {output}
        """
