#!/bin/bash

# Function to display error message and usage instructions
display_error_and_exit() {
    echo -e "\033[31mERROR:\033[0m $1"
    echo -e "Format: cp2k-restart.sh restart_file input_file"
    exit 1
}

# Check if required arguments are provided
if [ $# -ne 2 ]; then
    display_error_and_exit "Missing file operand! Please identify the name of restart and target files."
fi

# Check if restart file exists
if [ ! -f "$1" ]; then
    display_error_and_exit "Restart file '$1' not found."
fi

# Check if target file exists
if [ ! -f "$2" ]; then
    display_error_and_exit "Target file '$2' not found."
fi

# Function to extract and replace sections
function extract_and_replace_section {
    local section="$1"
    local start_marker="$2"
    local end_marker="$3"
    local src_file="$4"
    local target_file="$5"

    local src_start_line=$(grep -n "${start_marker}\$" "$src_file" | grep -v "#" | head -1 | awk -F ':' '{print $1}')
    local src_end_line=$(grep -n "${end_marker}\$" "$src_file" | grep -v "#" | tail -1 | awk -F ':' '{print $1}')

    local target_start_line=$(grep -n "${start_marker}\$" "$target_file" | grep -v "#" | head -1 | awk -F ':' '{print $1}')
    local target_end_line=$(grep -n "${end_marker}\$" "$target_file" | grep -v "#" | tail -1 | awk -F ':' '{print $1}')

    # Extracting the section from the source file
    sed -n "${src_start_line},${src_end_line}p" "$src_file" > tmp_${section}

    # Deleting the existing section from the target file
    sed -i "${target_start_line},${target_end_line}d" "$target_file"

    # Inserting the new section from the temporary file into the target file
    new_line_number=$(expr "${target_start_line}" - 1)
    sed -i "${new_line_number}r tmp_${section}" "$target_file"

    # Removing the temporary file
    rm "tmp_${section}"
}

# Function to extract velocity from one file and replace it in another file
extract_replace_velocity() {
    
    # Checking if the target file contains a velocity section
    check_velocity_target=$(grep "&VELOCITY" "$2" | grep -v "#")
    if [ "$check_velocity_target" ]; then
        extract_and_replace_section "VELOCITY" "&VELOCITY" "&END VELOCITY" $1 $2
    else
        # If the target file doesn't contain a velocity section, append the velocity at the end of COORD block
        # Extracting the line numbers for the velocity section in the restart files
        local src_start_line=$(grep -n "&VELOCITY\$" "$1" | grep -v "#" | head -1 | awk -F ':' '{print $1}')
        local src_end_line=$(grep -n "&END VELOCITY\$" "$1" | grep -v "#" | tail -1 | awk -F ':' '{print $1}')

        # Extracting the section from the source file
        sed -n "${src_start_line},${src_end_line}p" $1 > tmp_vel

        # Inserting the new section from the temporary file into the target file
        new_velocity_line_number=$(grep -n "&END COORD\$" "$2" | grep -v "#" | tail -1 | awk -F ':' '{print $1}')
        sed -i "${new_velocity_line_number}r tmp_vel" "$2"

        # Removing the temporary file
        rm tmp_vel
    fi
}

check_RUNTYPE_restart=`grep "RUN_TYPE" $1 |  grep -v "#" | awk -F ' ' '{print $NF}'`
check_RUNTYPE_target=`grep "RUN_TYPE" $2 |  grep -v "#" | awk -F ' ' '{print $NF}'`
if [ "$check_RUNTYPE_restart" == "MD" ] && [ "$check_RUNTYPE_target" == "MD" ]; then
    extract_and_replace_section "CELL" "&CELL" "&END CELL" $1 $2
    extract_and_replace_section "COORD" "&COORD" "&END COORD" $1 $2
    extract_and_replace_section "MD" "&MD" "&END MD" $1 $2
    extract_replace_velocity $1 $2
elif [ "$check_RUNTYPE_restart" == "MD" ] && [ "$check_RUNTYPE_target" != "MD" ]; then
    extract_and_replace_section "CELL" "&CELL" "&END CELL" $1 $2
    extract_and_replace_section "COORD" "&COORD" "&END COORD" $1 $2
elif [ "$check_RUNTYPE_restart" != "MD" ] && [ "$check_RUNTYPE_target" == "MD" ]; then
	extract_and_replace_section "CELL" "&CELL" "&END CELL" $1 $2
    extract_and_replace_section "COORD" "&COORD" "&END COORD" $1 $2
else
    extract_and_replace_section "CELL" "&CELL" "&END CELL" $1 $2
    extract_and_replace_section "COORD" "&COORD" "&END COORD" $1 $2
    # This may cause issue when you have both GEO_OPT and MD sections in the input file
    src_step_start_val_line_number=$(grep -n "STEP_START_VAL" "$1" | awk -F ':' '{print $1}')
    if [ "$src_step_start_val_line_number" ]; then
        # Finding the line number for the STEP_START_VAL in the target file
        target_step_start_val_line_number=$(grep -n "STEP_START_VAL" "$2" | grep -v "#  " |  awk -F ':' '{print $1}')

        # Extracting the value of STEP_START_VAL from the restart file
        src_step_start_val=$(grep "STEP_START_VAL" "$1")

        if [ "$target_step_start_val_line_number" ]; then
            # Replacing the STEP_START_VAL in the target file if it already exists
            sed -i "${target_step_start_val_line_number}c $src_step_start_val" "$2"
        else
            # If STEP_START_VAL doesn't exist in the target file, adding it after RMS_FORCE
            target_rms_force_line_number=$(grep -n "RMS_FORCE" "$2" | grep -v "#  " | awk -F ':' '{print $1}')
            sed -i "${target_rms_force_line_number}a $src_step_start_val" "$2"
        fi
    fi
fi
echo "Done!"
