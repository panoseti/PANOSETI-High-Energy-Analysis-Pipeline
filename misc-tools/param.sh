#!/bin/bash

# param.sh

# Assumes file structure like:
#   -data
#       --2026**
#           **.root.array

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

    # check flip
    local file="$2"
    local filename=$(basename "$file")
    local flip
    if [[ "$filename" == *_preflip* ]]; then
        flip="pre"
    elif [[ "$filename" == *_postflip* ]]; then
        flip="post"
    else
        echo ""
        return
    fi

    case "$date_str" in
        20260113)
            if [ "$flip" = "pre" ]; then
                echo "setCorrections(1,0,-0.85,2.63,0,0,0);"
                echo "setCorrections(2,0,-2.26,-11.49,0,0,0);"
                echo "setCorrections(3,0,-0.15,0.65,0,0,0);"
            else
                echo "setCorrections(1,0,0.62,-2.60,0,0,0);"
                echo "setCorrections(2,0,2.44,7.81,0,0,0);"
                echo "setCorrections(3,0,0.29,-1.15,0,0,0);"
            fi
            ;;
        20260114)
            if [ "$flip" = "pre" ]; then
                echo "setCorrections(1,0,-1.1,2.64,0,0,0);"
                echo "setCorrections(2,0,-2.27,-5.73,0,0,0);"
                echo "setCorrections(3,0,-0.19,0.57,0,0,0);"
            else
                echo "setCorrections(1,0,0.85,-2.54,0,0,0);"
                echo "setCorrections(2,0,2.44,3.37,0,0,0);"
                echo "setCorrections(3,0,0.07,-1.53,0,0,0);"
            fi
            ;;
        20260115)
            if [ "$flip" = "pre" ]; then
                echo "setCorrections(1,0,-0.98,2.84,0,0,0);"
                echo "setCorrections(2,0,-2.31,-6.51,0,0,0);"
                echo "setCorrections(3,0,-0.18,0.48,0,0,0);"
            else
                echo "setCorrections(1,0,0.90,-2.64,0,0,0);"
                echo "setCorrections(2,0,2.51,3.43,0,0,0);"
                echo "setCorrections(3,0,0.08,-1.31,0,0,0);"
            fi
            ;;
        20260116)
            if [ "$flip" = "post" ]; then
                echo "setCorrections(1,0,0.53,-2.31,0,0,0);"
                echo "setCorrections(2,0,2.23,1.16,0,0,0);"
                echo "setCorrections(3,0,0.09,-1.23,0,0,0);"
            fi
            ;;
        20260117)
            if [ "$flip" = "post" ]; then
                echo "setCorrections(1,0,0.82,-2.53,0,0,0);"
                echo "setCorrections(2,0,2.53,-1.03,0,0,0);"
                echo "setCorrections(3,0,0.44,1.08,0,0,0);"
            fi
            ;;
        20260118)
            if [ "$flip" = "post" ]; then
                echo "setCorrections(1,0,0.91,-2.48,0,0,0);"
                echo "setCorrections(2,0,2.53,-0.91,0,0,0);"
                echo "setCorrections(3,0,0.03,-1.33,0,0,0);"
            fi
            ;;
        20260119)
            if [ "$flip" = "pre" ]; then
                echo "setCorrections(1,0,-1.15,2.84,0,0,0);"
                echo "setCorrections(2,0,-2.27,0.77,0,0,0);"
                echo "setCorrections(3,0,-0.32,0.65,0,0,0);"
            else
                echo "setCorrections(1,0,0.68,-2.58,0,0,0);"
                echo "setCorrections(2,0,2.45,-1.17,0,0,0);"
                echo "setCorrections(3,0,0.03,-1.31,0,0,0);"
            fi
            ;;
        20260120)
            if [ "$flip" = "pre" ]; then
                echo "setCorrections(1,0,-1.5,2.84,0,0,0);"
                echo "setCorrections(2,0,-2.41,0.45,0,0,0);"
                echo "setCorrections(3,0,-0.29,0.66,0,0,0);"
            else
                echo "setCorrections(1,0,0.69,-2.57,0,0,0);"
                echo "setCorrections(2,0,2.48,-1.47,0,0,0);"
                echo "setCorrections(3,0,0.12,-1.23,0,0,0);"
            fi
            ;;
        20260121)
            if [ "$flip" = "pre" ]; then
                echo "setCorrections(1,0,-1.14,2.70,0,0,0);"
                echo "setCorrections(2,0,-2.42,0.34,0,0,0);"
                echo "setCorrections(3,0,-0.29,0.60,0,0,0);"
            else
                echo "setCorrections(1,0,0.88,-2.65,0,0,0);"
                echo "setCorrections(2,0,2.56,-1.65,0,0,0);"
                echo "setCorrections(3,0,0.28,-1.13,0,0,0);"
            fi
            ;;
    esac

}

# Function to process each file
process_file() {
    local file="$1"
    local date_dir="$2"
    local base="${file%.corr.array}"
    local corrections=$(get_corrections "$date_dir" "$file")

    # with corrections
    printf ".L panodisplay_REALDATA.C\nreadFile(\"%s\")\n%s\nparamCSV()\n" "$base" "$corrections" | root -l -b
    
    # without corrections
    # printf ".L panodisplay_REALDATA.C\nreadFile(\"%s\")\nparamCSV()\n" "$base" | root -l -b
}

# Export the function to make it available to GNU Parallel
export -f process_file
export -f get_corrections

# Loop over all date directories matching pattern
for date_dir in "$data_dir"/*/; do
    echo "Processing directory: $date_dir"
    find "$date_dir" -maxdepth 1 -type f -name "*flip*.root.corr.array" | parallel -j 8 process_file {} "$date_dir"
done
