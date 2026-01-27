#!/bin/bash

# param.sh

# Assumes file structure like:
#   -data
#       --202510**
#           --pcap
#               **.root.array

data_dir=$1

# Check if directory argument is provided
if [ -z "$data_dir" ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

# Function to get corrections for each date
get_corrections() {
    # extract date from folder name
    local date_dir="$1"
    local date_str=$(basename "$date_dir")
    
    # Extract time from file name
    local filename=$(basename "$file")
    local time_str=$(echo "$filename" | sed 's/.*_\([0-9]\{6\}\)\..*$/\1/')
    local time_int=$((10#$hour_min))

    case "$date_str" in
        20251017)
            if [ $time_int -lt 115200 ]; then
                echo "setCorrections(3,0,-3,-1.106,0,0,0);"
            else
                echo "setCorrections(3,0,3,1.106,0,0,180);"
            fi
            ;;
        20251018)
            if [ $time_int -lt 121034 ]; then
                echo "setCorrections(3,0,-3,-1.106,0,0,0);"
            else
                echo "setCorrections(3,0,3,1.106,0,0,180);"
            fi
            ;;
        20251019)
            if [ $time_int -lt 113040 ]; then
                echo "setCorrections(3,0,-2,-0.106,0,0,0);"
            else
                echo "setCorrections(3,0,2,0.106,0,0,180);"
            fi
            ;;
        20251020)
            if [ $time_int -lt 113040 ]; then
                echo "setCorrections(3,0,-2,-0.106,0,0,0);"
            else
                echo "setCorrections(3,0,2,0.106,0,0,180);"
            fi
            ;;
        20251021)
            if [ $time_int -lt 112631 ]; then
                echo "setCorrections(3,0,-2,-0.106,0,0,0);"
            else
                echo "setCorrections(3,0,2,0.106,0,0,180);"
            fi
            ;;
        20251022)
            if [ $time_int -lt 111929 ]; then
                echo "setCorrections(3,0,-2,-0.106,0,0,0);"
            else
                echo "setCorrections(3,0,2,0.106,0,0,180);"
            fi
            ;;
        20251023)
            if [ $time_int -lt 110000 ]; then
                echo "setCorrections(3,0,-2,-0.106,0,0,0);"
            else
                echo "setCorrections(3,0,2,0.106,0,0,180);"
            fi
            ;;
        20251024)
            if [ $time_int -lt 111407 ]; then
                echo "setCorrections(3,0,-2,-0.106,0,0,0);"
            else
                echo "setCorrections(3,0,2,0.106,0,0,180);"
            fi
            ;;
        20251028)
            if [ $time_int -lt 105835 ]; then
                echo "setCorrections(3,0,-2,-0.106,0,0,0);"
            else
                echo "setCorrections(3,0,2,0.106,0,0,180);"
            fi
            ;;
        20251029)
            if [ $time_int -lt 105435 ]; then
                echo "setCorrections(3,0,-2,-0.106,0,0,0);"
            else
                echo "setCorrections(3,0,2,0.106,0,0,180);"
            fi
            ;;
    esac

}

# Function to process each file
process_file() {
    local file="$1"
    local date_dir="$2"
    local base="${file%.array}"
    local corrections=$(get_corrections "$date_dir")

    # with corrections
    printf ".L panodisplay_REALDATA.C\nreadFile(\"%s\")\n%s\nparamCSV()\n" "$base" "$corrections" | root -l -b
    
    # without corrections
    # printf ".L panodisplay_REALDATA.C\nreadFile(\"%s\")\nparamCSV()\n" "$base" | root -l -b
}

# Export the function to make it available to GNU Parallel
export -f process_file
export -f get_corrections

# Loop over all date directories matching pattern
for date_dir in "$data_dir"/202510*/; do
    pcap_dir="${date_dir}pcap"
    echo "Processing directory: $pcap_dir"
    find "$pcap_dir" -type f -name "*.root.array" | parallel -j 8 process_file {} "$date_dir"
done
