# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 14:59:10 2022

@author: jkescher
"""
###############################################################################
# function created by Jan-Karl L. Escher for use in the course
# FENT2321/TEP4175-Design of a wind turbine, at NTNU.
# Modified by Loup Suja on 2022-08-31 to fir new airfoil polars format (for Ashes 3.18)
#
# DESCRIPTION: this function reads the polar files either from ASHES or airfoiltool,
# in order to create a database of Cl and Cd in function of Re-number and 
# angle of attack. This will allow to find the optimal alpha (this giving max Cl/Cd)
# for each Re-number.
#
# INPUT:
# The function takes two arguments: filename and source.
#
# filename: either a string or a list of strings. If the file originates from
# Ashes, the filename should containt the entire name of the file as Ashes files
# have data for all Reynolds numbers in the same file. In this case, it should
# be only one string.
# If it is for files originating from airfoiltools, there are two cases. If
# you want to import all files for the same airfoil from one directory, you
# enter the start of the name of all of the files (so without the number
# indicating the Reynolds number), and the function will import all text files
# in the current directory with the same start. If you want to specify a number
# of files to import, all of these have to be listed completely in a list.
#
# The files have to be located in the working directory, so in the folder from
# which you call the function.
#
# source: string that indicates whether the file originates form ashes or from
# airfoiltools. Should be either 'Airfoiltools' or 'Ashes'.
#
# OUTPUT:
#
# Clfun: function of type scipy.interpolate.interpnd.LinearNDInterpolator.
# it takes the Reynolds number and the angle of attack in the range of the
# interpolated data as arguments, and returns the corresponding Cl-value.
# Use in main module (script): give Re-number and alpha as input (e.g. Clfun(Re,alpha))
#
# Cdfun: same as Clfun, but with Cd
#
# AoAMaxfun: function which returns the optimal angle of attack at any
# given Re-number.
# Use in main module (script): give Re-number as input (e.g. AoAMaxfun(Re))
#
# Relims: returns the upper and lower limits of the Reynolds number.
# if you try to interpolate outside of this range, the functions will return
# errors.
#
# AoAlims: same as Relims, but with angle of attack

def importPolars(filename, source):
    # Import all of the relevant python modules
    import os
    import numpy as np
    import pandas as pd
    from scipy.interpolate import LinearNDInterpolator
    from scipy.interpolate import interp1d

    if source.upper() == 'ASHES':
        # Open the file and read all of its lines.
        with open(filename, "r") as bladefile:
            bladelines = bladefile.readlines()
        # Set initial values for the Re and the number of files to skip in the
        # header (Lines which are not commented out by an !, that is).
        ReStr = "Nan"
        # modification by Loup: bool to know whether the line being read is in the header or not
        headerBool = True
        # modification by Loup: bool to know whether we are reading a new Re number
        NewRe = False
        # create a temporary file, to be able to make a csv from our data later
        # The temporary file is supposed to contain five columns: Re, AoA, Cl,
        # Cd and Cm. Not that Cm is never used in this course.
        with open("temp.txt", "w") as tempTxt:
            # Iterate through all lines:
            for line in bladelines:
                # If the line is not commented out:
                # Modified by Loup: the comment symbol is no longer '!' but '#'
                if not (line[0] == "#"):
                    # if the line is "Re and corresponding polar", we're reading a new Re number and we're done with the header
                    if (line == "Re and corresponding polar\n"):
                        NewRe = True
                        headerBool = False
                    # the line following "Re and corresponding polar" is a Re number
                    # we add the Re number to ReStr (changing \n to \t)
                    elif(NewRe)==True:
                        ReStr = line[0:9]+"\t"
                        # we are not reading a Re anymore, we set the NewRe bool to false
                        NewRe = False
                    # if the line is not "Re and corresponding polar" and NewRe is false, then it contains aerodynamic data
                    elif(headerBool)==False:
                        if (line != "\n"):
                            tempTxt.write(ReStr+line)
        # now we interpret the txt as a comma separated values (csv) file, but
        # separated with tabs instead of commas.
        # We know what the columns are supposed to represent:
        colNames = ['Re', 'AoA', 'Cl', 'Cd', 'Cm']
        ClCdFrame = pd.read_csv("temp.txt", "\t",
                                names=colNames)

        ReList = np.asarray(ClCdFrame["Re"]).astype(float)
        AoAList = np.asarray(ClCdFrame["AoA"]).astype(float)
        ClList = np.asarray(ClCdFrame["Cl"]).astype(float)
        CdList = np.asarray(ClCdFrame["Cd"]).astype(float)
        AoAlims = [np.min(AoAList), np.max(AoAList)]
    # If the files are imported from airfoiltools, there is one file for each
    # Re-value. We will iterate through all of these, and extract data.
    elif source.upper() == 'AIRFOILTOOLS':
        AoA_max_list = []
        AoA_min_list = []
        # initialize all of the arrays
        ReList = np.array([])
        AoAList = np.array([])
        ClList = np.array([])
        CdList = np.array([])
        # two allowed cases: either all of the files to be imported are
        # Specified as lists, in that case the entire filename must be
        # specified, or only the filename without the reynolds number is given.
        # in that case, the program iterates through the current folder, to
        # find all .txt files starting with the same name.
        if type(filename) == list:
            file_iterator = filename
        elif type(filename) == str:
            file_iterator = []
            for file in os.listdir():
                if file.startswith(filename) and file.endswith('.txt'):
                    file_iterator.append(file)
        # Then we iterate through all the files:
        for file in file_iterator:
            lineCount = 0
            AoA_min = 360
            AoA_max = -360
            with open(file, 'r') as bladeFile:
                lines = bladeFile.readlines()
                for line in lines:
                    # we know on which line the Reynolds number is given:
                    if lineCount == 8:
                        Re = float(''.join(line.split()[5:8]))
                    # and on which line the actual blade data starts:
                    elif lineCount > 11:
                        # all of this is then saved in the lists.
                        ReList = np.append(ReList, Re)
                        AoA, Cl, Cd = line.split()[0:3]
                        AoAList = np.append(AoAList, float(AoA))
                        ClList = np.append(ClList, float(Cl))
                        CdList = np.append(CdList, float(Cd))
                        if float(AoA) < AoA_min:
                            AoA_min = float(AoA)
                        if float(AoA) > AoA_max:
                            AoA_max = float(AoA)

                    lineCount += 1
            AoA_min_list.append(AoA_min)
            AoA_max_list.append(AoA_max)
        AoAlims = [max(AoA_min_list), min(AoA_max_list)]

    else:
        raise NameError(f'{source} is not a valid value for source. Use \'ashes\' or \'airfoiltools\' instead')
    # Create a list of points with the Parameters, Re and AoA
    points = [(x, y) for (x, y) in zip(ReList, AoAList)]

    # Return functions for the lift and drag at different Re and different AoA
    Clfun = LinearNDInterpolator(points, ClList)
    Cdfun = LinearNDInterpolator(points, CdList)

    # In order to get the optimal angle of attack, we'll have to
    # do some more.
    Aoa = []
    Re = []
    RateOld = -1
    first = True
    # Entering a for loop, iterating through all the points.
    for (i, j) in points:
        # For the first iteration, we'll have to remember the first Re
        if first is True:
            ReOld = i
            first = False
            # If the next iteration has the same Re as the last iteration, we
            # check the Lift-drag ratio and keep this AoA if it's better than
            # the last best rate
        if i == ReOld:
            rate = Clfun(i, j)/Cdfun(i, j)
            if rate > RateOld:
                Aoa_current = j
                RateOld = rate
                # Otherwise save the Reynolds number together
                # with the best AoA and update ReOld.
        else:
            Re.append(ReOld)
            ReOld = i
            Aoa.append(Aoa_current)
            RateOld = Clfun(i, j)/Cdfun(i, j)
            # Save all of those values at the end of the loop.
    Re.append(ReOld)
    Aoa.append(Aoa_current)
    AoaMaxFun = interp1d(Re, Aoa)

    # And the limits for the reynolds number (upper and lower)
    Relims = [np.min(ReList), np.max(ReList)]

    return [Clfun, Cdfun, AoaMaxFun, Relims, AoAlims]
