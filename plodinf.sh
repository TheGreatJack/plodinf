#!/bin/bash

## DEFAULT VALUES



## Help message
plodinf_help() {
    echo "
    PLODINF - PLOidy INFerece for RAD-seq

    Authors:
	    Andres Florian Cruz
    	Juan Picon Cossio
	    Elisa Correa Pineda

    Version: 0.1

    Usage: plodinf.sh [options] FASTQ

    Options:
        -i        Input file. [required].
        -t        Threads. [4].
        -w        Working directory. Path to create the folder which will contain all mitnanex information. [./mitnanex_results].
        -r        Prefix name add to every produced file. [input file name].
        -d        Different output directory. Create a different output directory every run (it uses the date and time). [False].
        *         Help.
    "
    exit 1
}

while getopts 'i:t:m:w:r:d' opt; do
    case $opt in
        i)
        input_file=$OPTARG
        ;;
        t)
        threads=$OPTARG
        ;;
        w)
        wd=$OPTARG
        ;;
        r)
        prefix=$OPTARG
        ;;
        d)
        output_dir="plodinf_results_$(date  "+%Y-%m-%d_%H-%M-%S")/"
        ;;
        *)
        plodinf_help
        ;;
    esac 
done



get_snp_count (){
    ## Get the read count for each SNP
    bcftools mpileup -Q 20 -q 30 --skip-all-unset 3 -a "FORMAT/AD" -Ou --no-reference $bam_file 
}

filter_snp(){
    ## Filter SNPs by coverage and quality
    bcftools filter -Ou -e 'INFO/QS[1] == 1 | INFO/DP <= 10' -
}
    
format_snp_file (){
    ## Change the format for the outputted file
    bcftools query -f '%CHROM\t%POS\t%ALT\t%DP\tQS:%QS[\t%AD]]\n' 
}

compute_dist_snp(){
    ## Compute the distribution of reads to assign a ploidy base on that distribution
    python src/allele_counts_parser.py -m $minimum_ad -d $minimum_depth > $out_path/$prefix\_$pop\_allele_prop_m$minimum_ad\_d$minimum_depth.tsv
}
