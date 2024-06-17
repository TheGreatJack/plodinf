#!/bin/bash
#SBATCH --partition=longjobs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=13-00:00:00
#SBATCH --job-name=bwa-alignment_analysis
#SBATCH --cpus-per-task=32
#SBATCH -o result_%N_%j.out      # File to which STDOUT will be written
#SBATCH -e result_%N_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user=afflorianc@eafit.edu.co

#export SBATCH_EXPORT=NONE
#export OMP_NUM_THREADS=1

## Para apolo:
#module load miniconda3-4.10.3-gcc-11.2.0-vcglj27
module load python/3.10_miniconda-23.5.2

## Define variables
DATE=`date +%d-%m-%y_%H-%M`
env="bcftools"
threads=$SLURM_CPUS_PER_TASK
folder=$(pwd)


source activate $env


bam_path=""
#a_fist_pop_path=""
out_path="1"
pop="P1"
max_processes="24"
minimum_ad="15"
minimum_depth="25"

# Check if file path exists
if [ ! -d "$out_path" ]; then
    echo "proper_sorted_reads folder path does not exist."
    exit 1
fi

# Check if file path exists
if [ ! -d "$bam_path" ]; then
    echo "a_cepa_pop_path folder path does not exist."
    exit 1
fi

# Create an empty array to store file names
file_list=()

# Loop through files in the folder and add them to the array
for file in "$bam_path"/*.bam; do
    if [ -f "$file" ]; then
        prefix=$(basename $file _SORTED.bam)
        file_list+=("$prefix")
    fi
done

# Check if the array is empty
if [ ${#file_list[@]} -eq 0 ]; then
    echo "No files found in the bam folder."
    exit 1
fi


# Function to perform the task on a file
process_file() {
    prefix="$1"
    bam_path="$2"
    out_path="$3"
    pop="$4"
    minimum_ad="$5"
    minimum_depth="$6"

    echo "Processing file: $file"
    bam_file=$bam_path/$prefix\_SORTED.bam

    # Replace the following line with the command you want to execute on each file


    echo "START------------------------------"
    echo $bam_file
    echo $out_path/$prefix\_$pop\_allele_prop_m$minimum_ad\_d$minimum_depth.tsv

    bcftools mpileup -Q 20 -q 30 --skip-all-unset 3 -a "FORMAT/AD" -Ou --no-reference $bam_file | bcftools filter -Ou -e 'INFO/QS[1] == 1 | INFO/DP <= 10' - | bcftools query -f '%CHROM\t%POS\t%ALT\t%DP\tQS:%QS[\t%AD]]\n' | python ~/programs/sctacks_utils/allele_counts_parser.py -m $minimum_ad -d $minimum_depth > $out_path/$prefix\_$pop\_allele_prop_m$minimum_ad\_d$minimum_depth.tsv

    echo "Finished processing file: $prefix"
}

# Loop through files
for file in "${file_list[@]}"; do
    # Wait until the number of background processes is less than the maximum allowed
    while (( $(jobs -p | wc -l) >= max_processes )); do
        sleep 1
    done

    # Execute the process_file function in the background
    process_file "$file" "$bam_path" "$out_path" "$pop" "$minimum_ad" "$minimum_depth"  &

    # Optionally, you can add a delay between starting each process
    # to avoid launching all processes simultaneously
    # sleep 1
done

# Wait for all background processes to finish
wait

echo "All processes have completed."
