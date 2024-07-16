#!/usr/bin/env python
"""
A class to define the property of an atom.
"""
from Utils import *
from Element import *
from Bond import *

def CheckAtomLine(line):
    """
    Check if the line about atom is correct.
    """
    array = line.split()
    if len(array) < 4:
        print(("This line is too short : %s" % (line)))
        return False
    else:
        if (not array[0].isdigit()) and (not CheckElement(array[0])):
            print(("Failed to get atomic number or element symbol from line: %s" % (line)))
            return False
        for i in range(1, 4):
            try:
                float(array[i])
            except ValueError as err:
                print(("Error: %s. Incorrect coordinates " % (err)))
                print(("Failed to read atom in line: %s" % (line)))
                return False
            except:
                print(("Unknown Exception while reading atoms in line: %s " % (line)))
                return False
    return True

class Atom(Element):
    """
    Define a single atom with cartesian coordinates.
    """
    def __init__(self):
        """
        Define an empty atom.
        """
        Element.__init__(self)
        self.Coord = np.array([0.00]*3)
        self.Force = np.array([0.00]*3)
        # this will be useful in the molecule class.
        self.Index = 0
        self.CoordNum = 0
        # If the atom is adsorb site, then the AdsorbType if "A",
        # If the atom should be kept frozen, then the type should be "F".
        # If the atom belongs to the slab, then the type shoube be "S"
        # Else it is empty.
        # Only the first letter is meaningful.
        # Always use uppercase.
        self.Type = ""

    def ReadAtom(self, line):
        """
        Read the information of one atom from a line and store it as an Atom instance.
        Return False if something goes wrong; else return True.
        """
        # debug 
        # print "Reading line: %s" %(line)
        if not CheckAtomLine(line):
            print("Not a coorect atom line.")
            return False
        array = line.split()
        # Atom symbol or atomic number should be used.
        # Convert atomic number to atomic symbol.
        # Actually,is_integer function is better than isdigit, but is_integer is a new function in python 2.6
        # if (array[0].is_integer() == True):
        if array[0].isdigit():
            self.Symbol = ElementList[int(array[0]) - 1].Symbol
        else:
            self.Symbol = array[0]
        # Read the coordinate part of atom
        self.Coord = np.array([float(array[1]), float(array[2]), float(array[3])])
        self.AtomicMass = PeriodicTable[self.Symbol].AtomicMass
        if len(array) > 4:
            self.Type = array[4].upper()
        else:
            self.Type = ""
            # self.PrintAtom()
        return True

    def PrintAtom(self):
        """
        Friendly print the information of an atom,
        including atomic symbol, cartesian coordinates and atom type
        """
        print(("%s %.6f %.6f %.6f  %s" % (self.Symbol, self.Coord[0], self.Coord[1], self.Coord[2], self.Type)))

    def GetDist(self, atom2):
        """
        Get the distance between two atoms using their cartesian coordinates.
        """
        Dist = Length(self.Coord - atom2.Coord)
        return Dist

    def GetBond(self, atom2):
        """
        Calculate the bond type and bond distance of this atom and the other atom, and return a bond instance.
        """
        MyBond = Bond()
        atomnum1 = min(PeriodicTable[self.Symbol].AtomicNum, PeriodicTable[atom2.Symbol].AtomicNum)
        atomnum2 = max(PeriodicTable[self.Symbol].AtomicNum, PeriodicTable[atom2.Symbol].AtomicNum)
        BondType = ElementList[atomnum1 - 1].Symbol + "-" + ElementList[atomnum2 - 1].Symbol
        BondLength = self.GetDist(atom2)
        # Note: the vector is from atom to atom2.
        BondVector = Normalize(atom2.Coord - self.Coord)
        MyBond = Bond(BondType, BondLength, BondVector)
        return MyBond
