__author__ = 'Alessandro Orchini'

#! /usr/bin/python
"""Linear ThermoAcoustic Solver by Alessandro Orchini -
Started on 22-03-2015
Last Updated on 22-03-2015
This programs runs a GUI that enables the user to evaluate a Flame Transfer Function that then can be coupled with an acoustic network, or to directly build an eigenvalue problem that calculates the thermoacoutic modes with the largest eigenvalues.
"""

def main():
    varList = {'beta':1,'convSpeed':1,'Mark':0.02,'axi':1,'acModes':4,'Nr':800,'Tf':600,'xf':0.15}
    print varList

if __name__ == "__main__":
    main()
