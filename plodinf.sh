#!/bin/bash

## DEFAULT VALUES
minimum_ad=3
minimum_depth=12
threads=4
output_dir='mitnanex_results/'
wd="./"
minimum_coverage="10"

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
        -i        Input file. BAM file. [required].
        -t        Threads. [4].
        -w        Working directory. Path to create the folder which will contain all mitnanex information. [./plodinf_results].
        -r        Prefix name add to every produced file. [input file name].
        -z        Different output directory. Create a different output directory every run (it uses the date and time). [False].
        -m        Minimum number of reads per allele in a site. [3].
        -d        Minimum real depth per site. [12].
        -p        Minimum coverage. [10].
        -Q        Minimum mapping quality of an alignment. [20].
        -q        Minimum base quality. [30].
        *         Help.
    "
    exit 1
}

while getopts 'i:t:m:w:r:zm:d:p:Q:q' opt; do
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
        p)
        minimum_coverage=$OPTARG
        ;;
        Q)
        mapping_quality=$OPTARG
        ;;
        q)
        base_quality=$OPTARG
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

## PREFIX name to use for the resulting files
if [ -z $prefix ];
then 
    prefix=$(basename $input_file)
    prefix=${prefix%%.*}
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
    bcftools mpileup -Q 20 -q 30 --skip-all-unset 3 -a "FORMAT/AD" -Ou --no-reference $input_file 
}

filter_snp(){
    ## Filter SNPs by coverage and quality
    bcftools filter -Ou -e "INFO/QS[1] == 1 | INFO/DP < $minimum_coverage" $1
}
    
format_snp_file (){
    ## Change the format for the output file
    bcftools query -f '%CHROM\t%POS\t%ALT\t%DP\tQS:%QS[\t%AD]]\n' $1 -o $wd"/$prefix"".bfc"
}

compute_dist_snp(){
    ## Extract the number of reads supporting SNPs
    python src/allele_counts_parser.py -m $minimum_ad -d $minimum_depth > $wd"/_allele_prop_m"$minimum_ad"_d"$minimum_depth".tsv"
}


## Pipeline

#output=$(get_snp_count) && output=$(filter_snp $output) && format_snp_file $output