#!/usr/bin/env python
"""
Define the property of an element.
Read the PeriodicTable from external file.
The concept of Covalent radiusis from http://en.wikipedia.org/wiki/Covalent_radius
Data are from http://www.ccdc.cam.ac.uk/products/csd/radii/table.php4
These date will offer the first guess for the threshold of atomic bonds.
When the distance between two atoms is larger than the sum of their Covalent Radius,
there is a chemical bond between these two atoms.
"""
import os
DataDir = os.path.abspath(os.path.dirname(__file__))
PeriodicTableFileName = os.path.join(DataDir, "PeriodicTableOfElements.txt")

class Element(object):
    """
    Define an element
    """
    def __init__(self, AtomicNum=0, Symbol="", Name="", AtomicMass=0.0, CovalentRadius=0.0):
        """
        Define an empty element.
        """
        self.AtomicNum = AtomicNum
        self.Symbol = Symbol
        self.Name = Name
        self.AtomicMass = AtomicMass
        self.CovalentRadius = CovalentRadius

    def ReadElement(self, line):
        """
        Read element from line. Example:
        AtomicNumber, AtomicSymbol, ElementName, AtomicMass, CovalentRadius.
        1 H Hydrogen 1.0079 0.23
        """
        array = line.split()
        self.AtomicNum = int(array[0])
        self.Symbol = str(array[1])
        self.Name = str(array[2])
        self.AtomicMass = float(array[3])
        self.CovalentRadius = float(array[4])

    def PrintElement(self):
        """
        Friendly print the information about element
        """
        line = "%d\t%s\t%s\t%.4f\t%.2f" % (self.AtomicNum, self.Symbol, self.Name, self.AtomicMass, self.CovalentRadius)
        print(("%s" % (line)))

def ReadPeriodicTable(FileName):
    """
    Read all the information about elements from an external file.
    """
    # debug
    # print "Reading file %s" %(FileName)
    inp = open(FileName, "r")
    lines = inp.readlines()
    tmplist = []
    #Pop the first line.This is a descriptor
    #AtomicNumber Symbol Name AtomicMass CovalentRadius
    lines.pop(0)
    for line in lines:
        line = line.strip()
        array = line.split()
        if len(array) == 5:
            atom = Element()
            atom.ReadElement(line)
            tmplist.append(atom)
    tmpElementList = sorted(tmplist, key=lambda Element: Element.AtomicNum)
    inp.close()
    return tmpElementList

def PrintPeriodicTable(tmpElementList):
    """
    A function to print the Periodic Table.
    Just for debug.
    """
    print("AtomicNumber Symbol Name AtomicMass CovalentRadius")
    for element in tmpElementList:
        element.PrintElement()

def CreatePeriodicTable(tmpElementList):
    """
    Create the Periodic Table dict from the element list.
    """
    tmpPeriodicTable = {}
    for atom in tmpElementList:
        tmpPeriodicTable[atom.Symbol] = atom
    return tmpPeriodicTable

#Two Global Variables.
ElementList = ReadPeriodicTable(PeriodicTableFileName)
PeriodicTable = CreatePeriodicTable(ElementList)

def CheckElement(name):
    """
    Check if the name is the name of an element.
    """
    for atom in ElementList:
        if name == atom.Symbol:
            return True
    return False
