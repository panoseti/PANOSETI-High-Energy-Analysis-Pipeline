# #!/bin/bash

# paramdata.sh
# apply pointing offset corrections to October 2024 Lick Data

# ---------------
# IMPORTANT: this script assumes your data is named like:
# Pause_onsky_Crabph0pe_ima0pe__20241030_092327_698.pcapng.root.array
# ---------------

dir=$1

# Check if directory argument is provided
if [ -z "$dir" ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

# Function to process each file
process_file() {
    local file="$1"
    #echo -e ".L panodisplay_REALDATA.C\nreadFile(\"${file::-6}\")\nparamCSV()" | root -l -b
    # 11/1
    # echo -e ".L panodisplay_REALDATA.C
    # readFile(\"${file::-6}\")
    # if(std::stoi(\"${file: -28:6}\")<110000){
    #     setCorrections(1,1730447174,0,0,1.013e-3,2.533e-4,0);
    #     setCorrections(2,1730447174,0,0,0,9.193e-5,0);
    #     setCorrections(3,1730447174,0,0,-2.058e-4,0,0);
    # }else{
    #     setCorrections(1,1730459429,0,0,-8.239e-4,2.060e-4,180);
    #     setCorrections(2,1730458800,0,0,0,0,180);
    #     setCorrections(3,1730458800,0,0,0,0,180);
    # }
    # paramCSV()" | root -l -b

    # 10/31
    # echo -e ".L panodisplay_REALDATA.C
    # readFile(\"${file::-6}\")
    # setCorrections(1,1730359112,0,1,5.359e-4,0,0);
    # setCorrections(2,1730359099,1,0,0,0,0);
    # setCorrections(3,1730359112,1,0,-2.975e-4,0,0);
    # paramCSV()" | root -l -b
    
    # 10/30
    echo -e ".L panodisplay_REALDATA.C
    readFile(\"${file::-6}\")
    if(std::stoi(\"${file: -28:6}\")<110000){
        setCorrections(1,1730276185,0,0,0,0,0);
        setCorrections(2,1730275403,0,0,-1.987e-4,0,90);
        setCorrections(3,1730282071,0,0,0,0,0);
    }else{
        setCorrections(1,1730286027,-1,0,0,0,180);
        if(std::stoi(\"${file: -28:6}\")<113600){
            setCorrections(2,1730286555,4,0,0,0,180);
        }else{
            setCorrections(2,1730288167,0,0,0,0,180);
        }
        setCorrections(3,1730287607,0,0,0,-2.139e-4,180);
    }
    paramCSV()" | root -l -b
}

# Export the function to make it available to GNU Parallel
export -f process_file

# Use GNU Parallel to process files in parallel
# find "$dir" -type f -name "*.root.array" | parallel -j 10 process_file {}
find "$dir" -type f -name "*Crab*.root.array" | parallel -j 15 process_file {}
# find "$dir" -type f -name "*Boomerang*.root.array" | parallel -j 15 process_file {}

# Note: Adjust the number after -j to control the parallelism (e.g., -j 4 for 4 concurrent jobs)