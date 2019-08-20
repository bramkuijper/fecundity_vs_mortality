#!/usr/bin/env python3

import os,re,sys, traceback

# read files of Mathematica equation exports
# change variable names to corresponding things in C
# insert things in the corresponding C file

class Inserter:

    # the file used in C
    filename_cpp = "numsolve.cpp"

    paramstruct_name = "paramstruct"

    # specify the filenames containing the exports
    # from Mathematica
    filename_precur = "patchrecurs.txt" # patch frequency recursions
    filename_repval = "repvals.txt" # the recursions of the reproductive values
    filename_related = "relatedness.txt" # the relatedness coefficients
    filename_selgrads = "selgrad.txt" # the selection gradients
    filename_vars= "variables.txt" # a file containing all the variables used
    filename_params = "params.txt" # a file containing all the variables used

    # a list of regular expressions to transform Mathematica-expressions
    # to their C equivalents
    relist = ["r([am]{2,2})tplus1\((\d),(\d)\)","r([am]{2,2})\((\d),(\d)\)",  "s\((\d)\)", "f\((\d),(\d)\)", "ftplus1\((\d),(\d)\)", "v\((\d),(a|m),(\d)\)", "vtplus1\((\d),(a|m),(\d)\)", "Power", "sqrt", "ϵ"]
    sollist = ["r\\1tplus1\\2\\3", "r\\1\\2\\3", "s\\1", "f\\1\\2", "ftplus1\\1\\2", "v\\1\\2\\3","vtplus1\\1\\2\\3", "pow", "sqrt","epsilon"];

    # replaces all the Mathematica variables in the file filename
    # with their C equivalents according to the regexps above
    def transform(self, filename):

        # read in the file as a single string
        thefile = open(filename) 
        contents = thefile.read()
        thefile.close()

        # replace all occurrence of any of the patterns above
        for i in range(0, len(self.relist)):
            contents = re.sub(re.compile(self.relist[i]),self.sollist[i],contents)

        # return file contents
        return(contents)

    # function that inserts variable definitions at the required
    # places
    def insert_variables(self, strvars, strfile):

        # generate variable definition according to
        #double x = ((struct rparams *) params)->x;
        vardeffunc = re.sub("double\s(.*);","double \\1 = ((struct rparams *) params)->\\1;",strvars)

        strfile = re.sub("VARS","\n" + strvars,strfile)
        strfile = re.sub("VARFUNC","\n" + vardeffunc,strfile)
        
        return(strfile)

    # initialize the array of argc(x)
    # values at the start of the main() function
    # dependent on all the variables provided in vars(x)
    def initargument(self, strvars, strfile):
        
        strvars_split = strvars.split("\n")

        initstruct = ""

        for i in range(0, len(strvars_split)):
            initstruct += re.sub("double\s(.*);","\t\t" + self.paramstruct_name + ".\\1 = atof(argv[" + str(i + 2) + "]);\n",strvars_split[i].strip())


        strfile = re.sub("ARGINIT","\n" + initstruct, strfile)

        return(strfile)

    # generate the contents of the write params function
    # arguments: strvars - \n-separated string with all the variables
    # strfile: 
    def make_write_params_function(self, strparams, strfile):

        # split the string of params
        strparams_split = strparams.split("\n")

        filling = "DataFile << endl << endl "

        for param in strparams_split:

            if param == "":
                continue

            param_s = re.sub("double\s(.*);","\\1", param)

            # get all the parameters into a string with C++ code
            filling += " << \"" + param_s + ";\" << " + param_s + " << endl\n"

        filling += " << endl;"

        strfile = re.sub("WRITEPARAMS","\n" + filling, strfile)

        return(strfile)



    # generate the contents of the write data function
    # strvars: string with all the variables, still with types (i.e., double...)
    # strparams: a string with all the parameters, still with types (i.e., double...)
    # strfile: the file to write it to
    def make_write_data_function(self, strvars, strparams, strfile):
    
        # split all the variables according to line
        strvars_split = strvars.split("\n")

        # split all the parameters according to the line
        strparams_split = strparams.split("\n")

        # get those variables that are not parameters
        strvars_split = list(set(strvars_split) - set(strparams_split))


        # get all the variables in a string with C++ code
        function_filling = "DataFile << time << \";\" << "

        data_headers = "DataFile << \"time;"

        for var in strvars_split:
            if var == "":
                continue

            var_clean = re.sub("double\s(.*);","\\1", var) # strip the double before the variable declaration
            var_clean = var_clean.strip()


            function_filling += var_clean + " << \";\" << \n"
            data_headers += var_clean + ";"

        function_filling += " endl;"
        data_headers += "\" << endl;"

        strfile = re.sub("HEADERWRT","\n" + data_headers, strfile)

        strfile = re.sub("WRITEDATA","\n" + function_filling, strfile)

        return(strfile)


    # to calculate the eigenvalue, we need to specify the matrix
    # and then use gsl_eigen_nonsymmv
    # however we need to prepare this, which is done in this function
    def insert_matrix(self, dimension, strfile):

        strmat = "double data[] = {";

        for i in range(1,dimension+1):
            for j in range(1,dimension+1):
                strmat += "a" + str(i) + "_" + str(j) + ",\n";

        # cut off the trailing ",\n" that we needed in the for-loops
        strmat = strmat[0:-2] 

        # add a closing curly bracket
        strmat += "};\n\n"

        # replace it where necessary in the file
        strfile = re.sub("MATSPECIFY","\n" + strmat, strfile)

        return(strfile)

            
    def __init__(self):

        # transform all the files
        str_precur = self.transform(self.filename_precur)
        str_repval = self.transform(self.filename_repval)
        str_related = self.transform(self.filename_related)
        str_selgrads = self.transform(self.filename_selgrads)
        str_vars = self.transform(self.filename_vars)
        str_params = self.transform(self.filename_params)

        # open the c++ file
        cpp_f = open(self.filename_cpp)
        cpp_f_str = cpp_f.read()

        # put the contents of each of the Mathematica
        # file at their respective positions indicated in the c++ file
        cpp_f_str = re.sub("PATCHRECUR","\n" + str_precur,cpp_f_str)
#        cpp_f_str = re.sub("EVMAT","\n" + str_ev,cpp_f_str)
        cpp_f_str = re.sub("REPVALRECUR","\n" + str_repval,cpp_f_str)
        cpp_f_str = re.sub("RELRECUR","\n" + str_related,cpp_f_str)
        cpp_f_str = re.sub("SELGRADS","\n" + str_selgrads,cpp_f_str)

        # insert variable definitions
        cpp_f_str = self.insert_variables(str_vars, cpp_f_str)

        # insert matrix definition eigenvalue function
#        cpp_f_str = self.insert_matrix(4, cpp_f_str)

        cpp_f_str = self.initargument(str_vars, cpp_f_str)

        cpp_f_str = self.make_write_data_function(str_vars,str_params, cpp_f_str)
        cpp_f_str = self.make_write_params_function(str_params, cpp_f_str)
        print(cpp_f_str)

a = Inserter()
