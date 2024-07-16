#!/bin/bash

# Function to display error message and usage instructions
display_error_and_exit() {
    echo -e "\033[31mERROR:\033[0m $1"
    echo -e "Format: ./diag2ot.sh target_file"
    exit 1
}

# Check if source file exists
if [ ! -f "$1" ]; then
    display_error_and_exit "Target file '$1' not found."
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

# Set a function
function Multiwfn_run (){
    local Multiwfn_function_choice section start_marker end_marker

    case "$function_choice" in
        1)
            Multiwfn_function_choice="4"
            section="SCF"
            start_marker="&SCF"
            end_marker="&END SCF"
            ;;
        2)
            Multiwfn_function_choice=$'\n'
            section="SCF"
            start_marker="&SCF"
            end_marker="&END SCF"
            ;;
        3)
            Multiwfn_function_choice="-1\n3"
            section="OPT"
            start_marker="&MOTION"
            end_marker="&END MOTION"
            sed -i 's/\(RUN_TYPE\).*/\1 GEO_OPT/' "$1"
            sed -i '/STRESS_TENSOR/d' "$1"
            ;;
        4)
            Multiwfn_function_choice="-1\n4"
            section="OPT"
            start_marker="&MOTION"
            end_marker="&END MOTION"
            sed -i 's/\(RUN_TYPE\).*/\1 CELL_OPT/' "$1"
            sed -i '/&END FORCE_EVAL/i STRESS_TENSOR ANALYTICAL #Compute full stress tensor analytically' "$1"
            ;;
        5)
            Multiwfn_function_choice="-1\n6\n10\n2"
            section="MD"
            start_marker="&MOTION"
            end_marker="&END MOTION"
            sed -i 's/\(RUN_TYPE\).*/\1 MD/' "$1"
            sed -i '/STRESS_TENSOR/d' "$1"
            ;;
        *)
            echo "Invalid choice"
            ;;
    esac


	# Run Multiwfn 
	echo -e "$1\ncp2k\ntemp.inp\n$Multiwfn_function_choice\n0\nq" | Multiwfn 1>>/dev/null

    extract_and_replace_section "$section" "$start_marker" "$end_marker" "temp.inp" $1

	rm temp.inp
}


# Select function type
echo "Please select a function type:"
echo "1. OT"
echo "2. DIAGONALIZATION"
echo "3. GEO_OPT"
echo "4. CELL_OPT"
echo "5. MD"

read function_choice

Multiwfn_run $1
