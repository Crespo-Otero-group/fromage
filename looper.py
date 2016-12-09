# this scripts loops cryspy to feed it new excited geometries

import sys
import numpy
import subprocess
import os
from readfile import *
from editinputs import *

# convergence criterion for average difference in Mulliken charge
convCrit = 0.001

# index of the loop
loopN = 0

# name of project, should match the one in cryspy
name = "naphthalene222"
# same story with Gaussian dir
here = os.path.dirname(os.path.realpath(__file__))
gaussianDir = "GAUSSIAN"
gaussianPath = os.path.join(here, gaussianDir)

# will contain a list of lists of Mulliken charges
innerCharges = []

# carry on looping
stopBool = False
while stopBool == False:

    # run cryspy
    subprocess.call("python cryspy.py", shell=True)

    # prepare a new cryspy input
    subprocess.call("mv " + name + ".vasp " + name +
                    ".vasp." + str(loopN), shell=True)
    subprocess.call("mv " + name + ".new.vasp "+name + ".vasp", shell=True)

    # add the new Mulliken charges
    os.chdir(gaussianPath)
    newCharges = readGMull(name)["charges"]
    os.chdir(here)

    innerCharges.append(newCharges)

    # if it's not the first loop
    if loopN > 1:

        # compare the last two sets of charges
        ultiChar = innerCharges[-1]
        penultiChar = innerCharges[-2]

        # get average difference
        diffList = [abs(a - b) for a, b in zip(ultiChar, penultiChar)]
        avgDiff = sum(diffList) / len(diffList)

        # if the average is less than the tolerance, stop the loop
        if avgDiff < convCrit:
            stopBool = True
            print("Thank you for using Fruit Loops")
    loopN += 1
    print("Loop completed: "+str(loopN))
