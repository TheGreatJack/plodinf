#!/bin/bash

## DEFAULT VALUES
minimum_ad=3
minimum_depth=12
threads=4
output_dir='mitnanex_results/'
wd="./"

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
        -w        Working directory. Path to create the folder which will contain all mitnanex information. [./plodinf_results].
        -r        Prefix name add to every produced file. [input file name].
        -z        Different output directory. Create a different output directory every run (it uses the date and time). [False].
        -m        Minimum number of reads per allele in a site. [3].
        -d        Minimum real depth per site. [12].
        *         Help.
    "
    exit 1
}

while getopts 'i:t:m:w:r:zm:d:' opt; do
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
        z)
        output_dir="plodinf_results_$(date  "+%Y-%m-%d_%H-%M-%S")/"
        ;;
        m)
        minimum_ad=$OPTARG
        ;;
        d)
        minimum_depth=$OPTARG
        ;;
        *)
        plodinf_help
        ;;
    esac 
done

# Check if required arguments are provided
if [ -z "$input_file" ];
then
  echo "Error: Input file is required."
  plodinf_help
fi

## Add a slash if it is absent in the working directory
if [ ${wd: -1} = / ];
then 
    wd=$wd$output_dir
else
    wd=$wd"/"$output_dir
fi

## CREATE WORKING DIRECTORY
create_wd(){

    if [ -d $wd ]
    then
        echo $timestamp": Rewriting directory..."
        echo " "
    else 
        echo $timestamp": Creating directory..."
        echo " "
        mkdir $wd
    fi
}

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
    python src/allele_counts_parser.py -m $minimum_ad -d $minimum_depth > $wd"/_allele_prop_m"$minimum_ad"_d"$minimum_depth".tsv"
}
