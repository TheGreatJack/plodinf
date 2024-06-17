#!/bin/bash

## DEFAULT VALUES



## Help message
plodinf_help() {
    echo "
    PLODINF - PLOidy INFerece for RAD-seq

    Authors:
	Andres Florian Cruz
    	Juan Picon Cossio
	Elisa Correa Perez

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
        output_dir="mitnanex_results_$(date  "+%Y-%m-%d_%H-%M-%S")/"
        ;;
        *)
        plodinf_help
        ;;
    esac 
done
