#!/usr/bin/env python3
# summerizes a series of output file

import os, re, sys, math

# makes a list of all the parameters at the end of the file
def analyze_parameters(lines):
    
    # place parameters in dictionary
    pars = {}

    for line in lines:
        mobj = line.split(";")

        if len(mobj) > 1:
            pars[mobj[0]] = mobj[1].strip()

    return(pars)

# analyze a single file and get last line of output
# before the parameters start
def analyze_file(filename, print_header=False):

    # open the file and read in lines
    f = open(filename)
    fl = f.readlines()
    f.close

    # make list of line numbers and reverse it
    # so that we can go through all lines in reverse
    # looking for the last data line
    line_range = list(range(0,len(fl)))
    line_range.reverse()

    # aux variable to sot
    last_line_number = -1
    
    for line in line_range:

        # find a line that starts with a number (int or
        # floating point) and is finalized by a semicolon
        # after which another number appears -- definitely
        # output
        if re.match(r"\d+(\.?\d+)?;\d+.*",fl[line]) is not None:
            last_line_number = line
            break


    if last_line_number < 0:
        return()

    # get parameters present in file
    param_dict = analyze_parameters(fl[(last_line_number + 1):])
    
    # now output first and last lines

    # prepare the headers of this file
    if print_header:
        # first print the original headers that are already
        # in the output file
        header = fl[0].strip()

        # put parameter names as headers too
        param_header = ";".join(param_dict.keys())

        if header[-1] is not ";":
            header += ";"

        header += param_header

        # add additional values because we need to calculate
        # freq e2 and sbar from the s1 and s2 values
        header += ";freq_e2;sbar;"

        print(header)

    last_line = fl[last_line_number].strip()

    # add additional line for the x2 y2 values
    param_line = ";".join(param_dict.values()) 

    if last_line[-1] is not ";":
        last_line += ";"

    s1 = float(param_dict["s1"])
    s2 = float(param_dict["s2"])

    sbar = 0.5 * (math.log10(s1)  + math.log10(s2))

    freq_e2 = s1 / (s1 + s2)

    last_line += param_line 
    last_line += ";" + str(freq_e2) + ";" + str(sbar) + ";"

    print(last_line)

if len(sys.argv) < 2:
    print("Please provide a directory name in which the output files are located")
    exit(1)

first_file = True

# now loop through files
for root, dir, files in os.walk(sys.argv[1]):

    for file in files:
        # file is simulation file?
        if re.match("(sim|iter).*\d$",file) != None:
            analyze_file(
                    filename=os.path.join(root, file)
                    ,print_header=first_file)
            exit(1)
            if first_file:
                first_file = False
