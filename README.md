# SV_filtering

This simple workflow filters vcf files based on provided cases, and automatically detects remaining samples as controls.
To run the workflow, you must provide

* your input VCF containing multiple samples
* a cases file
* file names and output directory locations in the config file

## Cases file

to use the tool, you need a file specifying cases.
Please add it in the format "cases/{species}_chr{number}.txt"
the file should be a text file with a new line for each sample, i.e.:

case1
case2
case3

## Filtering modes

Two modes are included in this workflow: 

- optimized for recessive filtering (more strict) 
- optimized for dominant filtering (less strict) 

Please specify the mode in the config file.

If recessive is selected, only variants which are present as 1/1 and 1|1 in *all* cases will be kept. 
All variants with 1/1 or 1|1 in *any* control will be discarded

If dominant is selected, all variants with at least one alternative allele (1/0, 0/1, 1|0, 0|1, 1/1, 1|1) in *all* cases will be kept
The variant will be filtered out if it has anything other than 0/0 or 0|0 in *any* control.
