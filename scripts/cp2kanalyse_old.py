# Usage: python cp2kanalysis.py cp2k.out

"""
This script analyzes a CP2K output file and extracts information about the
geometry optimization, including the number of SCF cycles, energy change, and
convergence criteria.

Usage:
  cp2k-output-analyse.py OUTPUT_FILE.out
"""

import os
import sys
import datetime


def print_usage():
    print("""
***************************************
*** The format of the script: ***
cp2kanalyse.py OUTPUT_FILE.out
***************************************
    """)


def process_output_file(output_file):
    plot_file = output_file.split('.out')[0] + "__data.csv"
    starttime = ""
    TOTAL_TIME = 0
    line_number = 0
    MAX_D = "0.003"
    RMS_D = "0.0015"
    MAX_F = "0.00045"
    RMS_F = "0.0003"

    with open(output_file, 'r') as f, open(plot_file, 'w') as o:
        lines = f.readlines()
        num_lines = len(lines)
        MAX_SCF = 50
        OUTER_SCF_CHECK = "FALSE"
        SCF_OPTIMIZER = "DIAGONALIZATION"
        GEO_OPTIMIZER = "N/A"
        for line in lines:
            line_number += 1
            if "PROGRAM STARTED AT" in line:
                starttime = line.split('AT')[-1].strip()
            elif "PROGRAM ENDED AT" in line:
                endtime = line.split('AT')[-1].strip()
                TOTAL_TIME = (datetime.datetime.strptime(endtime, "%Y-%m-%d %H:%M:%S.%f") \
                        - datetime.datetime.strptime(starttime, "%Y-%m-%d %H:%M:%S.%f")).total_seconds()
            elif "Run type" in line:
                RUN_TYPE = line.split()[-1]
                if RUN_TYPE == "GEO_OPT":
                    line_added = 25
                elif RUN_TYPE == "CELL_OPT":
                    line_added = 32
                else:
                    print("\033[31mERROR:\033[0m This script can only be used for Geometry Optimization results!")
                    sys.exit()
            elif "eps_scf:" in line:
                EPS_SCF = line.split(':')[-1].strip()
            elif "Outer loop SCF in use" in line:
                OUTER_SCF_CHECK = "TRUE"
            elif "max_scf:" in line:
                MAX_SCF = line.split(':')[-1].strip()
            elif "STARTING GEOMETRY OPTIMIZATION" in line:
                line = line_number
                line += 1
                GEO_OPTIMIZER = lines[line_number].split('***')[1].strip()
            elif " OT " in line:
                SCF_OPTIMIZER = "OT"
            elif "outer SCF loop FAILED to converge" in line:
                SCF_NUMBER_OT = line.split()[-2]
                SCF_NUMBER = "!" + SCF_NUMBER_OT
            elif "SCF run NOT converged ***" in line:
                if OUTER_SCF_CHECK == "FALSE":
                    SCF_NUMBER = "!" + MAX_SCF
            elif "*** SCF run converged" in line:
                SCF_NUMBER = line.split()[-3]
            elif "outer SCF loop converged in" in line:
                SCF_NUMBER = line.split()[-2]
            elif "Informations at step" in line:
                CYCLE_NUMBER = line.split('=')[-1].split('-')[0].strip()
                top_line_number = line_number
                bottom_line_number = top_line_number + line_added
                for contents in lines[top_line_number:bottom_line_number]:
                    if "Decrease in energy " in contents:
                        ENERGY_CHANGE = contents.split('=')[-1].strip()
                    elif "Total Energy" in contents:
                        TOTAL_ENERGY = float(contents.split('=')[-1].strip())
                    elif "Real energy change" in contents:
                        ENERGY_CHANGE_VALUE = float(contents.split('=')[-1].strip())
                        ENERGY_CHANGE_VALUE = "{:.2e}".format(ENERGY_CHANGE_VALUE)
                    elif "Conv. limit for step size" in contents:
                        MAX_D = "0" + contents.split('=')[-1].strip().strip('0')
                        MAX_D = str(round(float(MAX_D),3))
                    elif "Conv. limit for RMS step" in contents:
                        RMS_D = "0" + contents.split('=')[-1].strip().strip('0')
                        RMS_D = str(round(float(RMS_D),4))
                    elif "Conv. limit for gradients" in contents:
                        MAX_F = "0" + contents.split('=')[-1].strip().strip('0')
                        MAX_F = str(round(float(MAX_F),5))
                    elif "Conv. limit for RMS grad" in contents:
                        RMS_F = "0" + contents.split('=')[-1].strip().strip('0')
                        RMS_F = str(round(float(RMS_F),4))
                    elif "Max. step size " in contents:
                        MAX_D_VALUE = float(contents.split('=')[-1].strip())
                    elif "RMS step size " in contents:
                        RMS_D_VALUE = float(contents.split('=')[-1].strip())
                    elif "Max. gradient " in contents:
                        MAX_F_VALUE = float(contents.split('=')[-1].strip())
                    elif "RMS gradient " in contents:
                        RMS_F_VALUE = float(contents.split('=')[-1].strip())
                    elif "Used time" in contents:
                        USEDTIME = contents.split('=')[-1].strip()
                        USEDTIME = str(round(float(USEDTIME)))
                        TOTAL_TIME += float(USEDTIME)
                try:
                    ENERGY_CHANGE = ENERGY_CHANGE
                    if ENERGY_CHANGE == "NO":
                        o.write("%1s %4s |%4s |%15.8f |%10s" % ("x", CYCLE_NUMBER, SCF_NUMBER, TOTAL_ENERGY,
                                                                  ENERGY_CHANGE_VALUE))
                    else:
                        o.write("%6s |%4s |%15.8f |%10s" % (CYCLE_NUMBER, SCF_NUMBER, TOTAL_ENERGY,
                                                            ENERGY_CHANGE_VALUE))
                except NameError:
                    o.write("%6s |%4s |%15.8f |%7s   " % (CYCLE_NUMBER, SCF_NUMBER, TOTAL_ENERGY, "N/A"))
                try:
                    o.write(" |%7.4f" % (MAX_D_VALUE))
                    if MAX_D_VALUE > float(MAX_D):
                        MAX_D_CONVERGENCE = "NO"
                    else:
                        MAX_D_CONVERGENCE = "YES"
                    o.write("%4s" % (MAX_D_CONVERGENCE))
                    o.write(" |%8.4f" % (RMS_D_VALUE))
                    if RMS_D_VALUE > float(RMS_D):
                        RMS_D_CONVERGENCE = "NO"
                    else:
                        RMS_D_CONVERGENCE = "YES"
                    o.write("%4s" % (RMS_D_CONVERGENCE))
                    o.write(" |%9.5f" % (MAX_F_VALUE))
                    if MAX_F_VALUE > float(MAX_F):
                        MAX_F_CONVERGENCE = "NO"
                    else:
                        MAX_F_CONVERGENCE = "YES"
                    o.write("%4s" % (MAX_F_CONVERGENCE))
                    o.write(" |%8.4f" % (RMS_F_VALUE))
                    if RMS_F_VALUE > float(RMS_F):
                        RMS_F_CONVERGENCE = "NO"
                    else:
                        RMS_F_CONVERGENCE = "YES"
                    o.write("%4s" % (RMS_F_CONVERGENCE))
                    o.write(" |%6s" % (USEDTIME))
                except NameError:
                    o.write(" |%8s    |%8s     |%8s      |%8s    " % ("N/A", "N/A", "N/A", "N/A"))
                    o.write(" |%6s" % (USEDTIME))
                o.write("\n")
            elif "OPTIMIZATION COMPLETED" in line:
                CYCLE_NUMBER = line.split('=')[-1].split('-')[0].strip()
                top_line_number = line_number
                bottom_line_number = num_lines
                for contents in lines[top_line_number:bottom_line_number]:
                    if "ENERGY" in contents:
                        TOTAL_ENERGY = float(contents.split(':')[-1].strip())

                o.write("%6s |%4s |%15.8f" % ("Final", SCF_NUMBER, TOTAL_ENERGY))
                o.write(" |%7s    |%8s    |%8s     |%8s      |%8s    " % ("N/A", "N/A", "N/A", "N/A", "N/A"))
                o.write(" |%6s" % ("N/A"))
                o.write("\n")
        TOTAL_TIME = str(datetime.timedelta(seconds=round(float(TOTAL_TIME))))
        o.write("# Done!")
 
    with open(plot_file, 'r+') as f:
        contents = f.read()
        f.seek(0, 0)
        f.write("# Job Starting Date: " + starttime \
                + "\n# Total used time: " + str(TOTAL_TIME) \
                + "\n# Directory: " + os.getcwd() \
                + "\n# RUN_TYPE: " + RUN_TYPE \
                + "\n# EPS_SCF: " + EPS_SCF \
                + "\n# MAX_SCF: " + MAX_SCF \
                + "\n# SCF_OPTIMIZER: " + SCF_OPTIMIZER \
                + "\n# OUTER_SCF: " + OUTER_SCF_CHECK \
                + "\n# GEO_OPTIMIZER: " + GEO_OPTIMIZER \
                + "\n# STEP | SCF |    E [a.u.]    |  Delta E  | M_D(" + MAX_D + ") | R_D(" + RMS_D + ") | M_F(" + MAX_F + ") | R_F(" + RMS_F + ") | TIME [s]" \
                + "\n" + contents)


def plot_cycle_vs_energy(plot_file):
    import math
    import re
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import MultipleLocator
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    with open(plot_file, 'r') as f:
        x = []
        y = []
        for lines in f:
            value = re.split(r'\s*\|\s*|\s+', lines.strip())
            if value[0].isdigit():
                x.append(int(value[0]))
                y.append(float(value[2]))
            elif value[0] == "x":
                x.append(int(value[1]))
                y.append(float(value[3]))
            elif value[0] == "Final":
                x.append(x[-1] + 1)
                y.append(float(value[2]))

        plt.figure(figsize=(10, 6))
        plt.grid(True)
        plt.plot(x, y, '-o', label='Energy')
        plt.xlabel("Ion steps")
        plt.ylabel("Energy (a.u.)")
        plt.xlim(-0.5, max(x)+0.5)
        plt.title('Energy of each ion steps')
        plt.legend()
        plt.savefig('./cp2k.png', bbox_inches='tight')


def main():
    if len(sys.argv) < 2:
        print("\033[31mERROR:\033[0m Missing file operand! Please identify the name of OUTPUT_FILE.out")
        print_usage()
        sys.exit(1)

    output_file = sys.argv[1]
    print("=======================================")
    print("In process...")
    print("...")

    try:
        process_output_file(output_file)
    except FileNotFoundError:
        print(f"Error: File '{output_file}' not found.")
        sys.exit(1)

    print("=======================================")
    print("Do you want to plot cycle vs. energy?\n(y/n)")
    plot_choice = input()
    if plot_choice.lower() == "y":
        plot_cycle_vs_energy(output_file.split('.out')[0] + "__data.csv")

    print("Done!")
    print("=======================================")


if __name__ == "__main__":
    main()


