#!/bin/bash
# count_variants.sh
#
# Usage: ./count_variants.sh -i input.bam [-o output.csv]
#
# This script:
#  - Takes a BAM file as input (-i) and uses its base name (without .bam) as the sample name.
#  - Uses samtools to stream SAM records.
#  - Uses an embedded AWK program to count occurrences of specific variant sequences
#    (or their reverse complements) in the read sequences.
#
# The target sequences are:
#   MUT_27dupC     : GGGCTCCACCGCCCCCCCCAGCCCACGGTGTC
#   MUT_27insCCCC  : GGGCTCCACCGCCCCCCCCCCCAGCCCACGGTGTC
#   MUT_26_27insG  : GGGCTCCACCGCCCCCCGCAGCCCACGGTGTC
#   MUT_28dupA     : GGGCTCCACCGCCCCCCCAAGCCCACGGTGT
#   MUT_23delinsAT : GGCTCCACCGCCATCCCCAGCCCACGGTGTC
#   WT             : GGGCTCCACCGCCCCCCCAGCCCACGGTGTC
#   All            : GGGCTCCACCG.*CCCCCCCAGCCCACGGTGTC
#
# If an output file is specified and already exists, the script will:
#  - Overwrite the row for the sample if it is already present.
#  - Otherwise, append the new row.
#
# If no output file is provided, output is written to stdout.

# Function to display usage.
usage() {
    echo "Usage: $0 -i input.bam [-o output.csv]" >&2
    exit 1
}

# Parse command-line options.
while getopts "i:o:" opt; do
  case "$opt" in
    i)
      bamfile="$OPTARG"
      ;;
    o)
      outfile="$OPTARG"
      ;;
    *)
      usage
      ;;
  esac
done

# Check that an input file was provided.
if [ -z "$bamfile" ]; then
    echo "Error: Input BAM file (-i) is required." >&2
    usage
fi

# Check that the BAM file exists.
if [ ! -f "$bamfile" ]; then
    echo "Error: File '$bamfile' not found." >&2
    exit 1
fi

# Derive the sample name by stripping the .bam extension.
sampleName=$(basename "$bamfile" .bam)

# Determine if we are writing to a file.
if [ -n "$outfile" ]; then
    if [ -f "$outfile" ]; then
        # File exists: set flag so that AWK prints only the row (no header).
        appendFlag=1
    else
        # File does not exist: we'll print header + row.
        appendFlag=0
    fi
fi

# Store the AWK program in a variable.
awk_prog='
function revcomp(seq,    i, len, res, c) {
    len = length(seq)
    res = ""
    for (i = len; i > 0; i--) {
        c = substr(seq, i, 1)
        if (c ~ /[Aa]/) {
            res = res "T"
        } else if (c ~ /[Tt]/) {
            res = res "A"
        } else if (c ~ /[Gg]/) {
            res = res "C"
        } else if (c ~ /[Cc]/) {
            res = res "G"
        } else {
            res = res c
        }
    }
    return res
}

BEGIN {
    targets["MUT_27dupC"]     = "GGGCTCCACCGCCCCCCCCAGCCCACGGTGTC";
    targets["MUT_27insCCCC"]  = "GGGCTCCACCGCCCCCCCCCCCAGCCCACGGTGTC";
    targets["MUT_26_27insG"]  = "GGGCTCCACCGCCCCCCGCAGCCCACGGTGTC";
    targets["MUT_28dupA"]     = "GGGCTCCACCGCCCCCCCAAGCCCACGGTGT";
    targets["MUT_23delinsAT"] = "GGCTCCACCGCCATCCCCAGCCCACGGTGTC";
    targets["WT"]             = "GGGCTCCACCGCCCCCCCAGCCCACGGTGTC";
    targets["All"]            = "GGGCTCCACCG.*CCCCCCCAGCCCACGGTGTC";
    for (t in targets) {
         counts[t] = 0;
    }
}
$0 ~ /^@/ { next }
{
    seq = $10;
    seq_rc = revcomp(seq);
    for (name in targets) {
         pattern = targets[name];
         if (seq ~ pattern || seq_rc ~ pattern) {
              counts[name]++;
         }
    }
}
END {
    # If not appending (i.e. new file), print header.
    if (append == 0) {
        printf "Sample"
        for (name in targets) {
             printf "\t%s", name;
        }
        printf "\n";
    }
    # Always print the sample row.
    printf "%s", sample;
    for (name in targets) {
         printf "\t%d", counts[name];
    }
    printf "\n";
}
'

# Function to generate the output (header + row or just row).
generate_output() {
    if [ -n "$outfile" ]; then
        samtools view "$bamfile" | awk -v sample="$sampleName" -v append="$appendFlag" "$awk_prog"
    else
        samtools view "$bamfile" | awk -v sample="$sampleName" -v append=0 "$awk_prog"
    fi
}

if [ -n "$outfile" ]; then
    # Generate the new output (as one or more lines).
    new_output=$(generate_output)
    # Separate header (if any) from the row.
    header=$(echo "$new_output" | head -n 1)
    row=$(echo "$new_output" | tail -n 1)
    
    if [ ! -f "$outfile" ]; then
        # If the file doesn't exist, write header and row.
        echo -e "$header" > "$outfile"
        echo -e "$row" >> "$outfile"
        echo "Output written to $outfile"
    else
        # File exists. Check if a row for this sample already exists.
        if grep -q "^${sampleName}[[:space:]]" "$outfile"; then
            # Replace the existing row.
            # Using sed to replace the line starting with sampleName.
            sed -i.bak "s/^${sampleName}[[:space:]].*/${row}/" "$outfile"
            echo "Row for sample '${sampleName}' updated in $outfile"
        else
            # Append the new row.
            echo -e "$row" >> "$outfile"
            echo "Row for sample '${sampleName}' appended to $outfile"
        fi
    fi
else
    # No output file provided: print to stdout.
    generate_output
fi
