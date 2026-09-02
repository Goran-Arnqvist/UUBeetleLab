#!/bin/bash -l
#example batch script for calculating depth of coverage from bam files using mosdepth


# mosdepth and parallel msut be in $PATH 

# Function to display usage instructions
usage() {
  echo -e "Script for creating depth of coverage files from bam files using mosdepth"
  echo "  -i  Input directory containing BAM files"
  echo "  -o  Output directory to store results"
  echo "  -n  Naming scheme for output files (default: depth_coverage)"
  echo "  --by <window_size>    Specify window size for mosdepth --by option"
  echo "  -Q <mapq>             Mapping quality threshold (default: 0)"
  echo "  --include-flag <FLAG> Include reads with any of the bits in FLAG set"
  echo "  --use-scratch         Use scratch for temporary files"
  exit 1
}

# Initialize variables
input_dir=""
output_dir=""
flag_scratch=false
threads_available=$(nproc)
output_naming_scheme="depth_of_coverage"
by_window_size=""
mapq_threshold="0"
include_flag=""

# Parse options (flags and arguments that start with -)
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -i)
      input_dir="$2"
      shift 2
      ;;
    -o)
      output_dir="$2"
      shift 2
      ;;
    -n)
      output_naming_scheme="$2"
      shift 2
      ;;
    --by)
      by_window_size="$2"
      shift 2
      ;;
    -Q)
      mapq_threshold="$2"
      shift 2
      ;;
    --include-flag)
      include_flag="$2"
      shift 2
      ;;
    --use-scratch)
      flag_scratch=true
      shift
      ;;
    *)
      echo "Unknown option: $1"
      usage
      ;;
  esac
done

# Validate required arguments
if [[ -z "$input_dir" || -z "$output_dir" ]]; then
  echo "Error: Both input and output directories must be specified."
  usage
fi

# Deal with relative paths for directories
input_dir=$(realpath "$input_dir")
output_dir=$(realpath "$output_dir")

# Debugging: Print out received arguments
echo "Input directory: $input_dir"
echo "Output directory: $output_dir"
echo "Use scratch: $flag_scratch"
echo "Window size (--by): ${by_window_size:-Not specified}"
echo "Mapping quality threshold (-Q): $mapq_threshold"
echo "Include flag (--include-flag): ${include_flag:-Not specified}"

# Create output directory if it does not exist
mkdir -p "$output_dir" || exit 1

# Set working directory
working_dir="$input_dir"

cd "$working_dir" || exit 1

if [ "$flag_scratch" = true ]; then
  echo -e "Using scratch for temporary files.\n"

  mkdir -pv "$SNIC_TMP"/working || exit 1
  mkdir -pv "$SNIC_TMP"/output || exit 1

  working_dir="$SNIC_TMP"/working

  # Match bam files and index files
  rsync -Lv ./*.bam* "$working_dir" || exit 1
  cd "$working_dir" || exit 1
fi

# Build mosdepth options
mosdepth_options=""

if [[ -n "$by_window_size" ]]; then
  mosdepth_options+=" --by $by_window_size"
fi

if [[ -n "$mapq_threshold" ]]; then
  mosdepth_options+=" -Q $mapq_threshold"
fi

if [[ -n "$include_flag" ]]; then
  mosdepth_options+=" -i $include_flag"
fi

# Export variables needed in parallel block
export output_naming_scheme
export mosdepth_options

# Create coverage files
# Parallel if more than one bam file
parallel --jobs "$threads_available" --env mosdepth_options --env output_naming_scheme '
  file={1}

  mosdepth -t 1 $mosdepth_options "$(basename "$file" .bam)_$output_naming_scheme" "$file"

  echo "Coverage calculation for $file is done."
' ::: ./*.bam

echo -e "Coverage calculation for all bam files is done.\n"

# Move output files to output directory
if [ "$flag_scratch" = true ]; then
  rsync -Lv ./*.mosdepth* "$output_dir"
  rsync -Lv ./*bed* "$output_dir"
else
  mv -v ./*.mosdepth* "$output_dir"
  mv -v ./*bed* "$output_dir"
fi

echo -e "Output files moved to $output_dir.\n"

echo -e "All done.\n"
exit 0
