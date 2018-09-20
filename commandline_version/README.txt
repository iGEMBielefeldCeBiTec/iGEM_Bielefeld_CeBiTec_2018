Usage:
	python siRNA_commandline.py
	
	required:
		-i <FULL PATH TO INPUT FILE>		gene sequnce
		-m <siRNA | RNAi>			desifn method
		-c					check siRNA
		
	if the check function is used a second input is required:
		-si <siRNA sequence | FULL PATH TO INPUT FILE>		siRNA sequence
	
	optionally:		  
		-v (Bielefeld Method)			if the Bielefeld vector is used
		-o <FULL PATH TO OUTPUT FILE>		output file
	
	All input files should be in FASTA format:
		> header
		gene/siRNA sequence
