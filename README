# Synopsis

`moreCigar(read.cigar, read.opt("MD"))` 

# Goal 

It would enable `read.cigar` in pysam to absorb `MD` tag information. The syntax of pysam could be found [here](http://www.cgat.org/~andreas/documentation/pysam/api.html?highlight=cigar#pysam.AlignedRead.cigar) within which `0` denotes **alignment matched** position. However, though `7` and `8` are claimed to denote equal and mismatch respectively, the **sequence** mismatched and matched base in `bam` file with `MD` tag attached and even with identical base changed to `=` cannot be distinguished and are still mixed to be marked by `0`. 

# Demo

	>>> ## Demo
	>>> demoMD = "31^AT3T14"
	>>> ## Demo read.cigar()
	>>> cigarList = [(0, 31), (2,2), (3, 94), (0, 18)]
	>>> ## cigar list with match/mismatch information: "
	>>> moreCigarList = moreCigar(cigarList, demoMD)
	>>> print moreCigarList
	[(7, 31), (2, 2), (3, 94), (7, 3), (8, 1), (7, 14)]

# Misc

`samtools` could fill `MD` tag if your aligner did not write `MD` tag information by using: `samtools calmd`. See [here](http://samtools.sourceforge.net/samtools.shtml) for its usage. 