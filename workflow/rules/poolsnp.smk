rule poolsnp:
	input:
	output:
	shell:
	"""
	PoolSNP.sh   \
	mpileup=input.mpileup.gz \
	reference=reference.fasta.gz \
	names=Sample1,Sample2,Sample3,Sample4 \
	max-cov=0.98 \
	min-cov=10 \
	min-count=10 \
	min-freq=0.01 \
	miss-frac=0.2 \
	jobs=22 \
	badsites=1 \
	allsites=0 \
	output=output-file
	"""