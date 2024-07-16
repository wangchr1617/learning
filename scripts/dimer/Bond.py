#!/usr/bin/env python
"""
A class to describe chemical bond.
Not very useful.
"""
from Utils import *
from Element import *

def InitBondRange():
    """
    Two atoms are considered to be bonded when their distance is smaller than the bondrange.
    A tolerence of about 0.4 angstrom should be allowed.
    If this tolerence is not enough, then the bondrange should be defined elsewhere.
    Make a table of BondRange using the CovalentRadius map.
    Just include the first 100 elements(From H to Fm). I don't care about other heavier elements.
    """
    BondRange = {}
    for i in range(0, 101):
        for j in range(i, 101):
            atom1 = ElementList[i].Symbol
            atom2 = ElementList[j].Symbol
            BondType1 = atom1 + "-" + atom2
            BondType2 = atom2 + "-" + atom1
            BondLength = ElementList[i].CovalentRadius + ElementList[j].CovalentRadius
            BondRange[BondType1] = BondLength
            BondRange[BondType2] = BondLength
    return BondRange

class Bond(object):
    """
    A class to define distance in molecule.
    """
    def __init__(self, BondType="", Dist=0, Vec=np.array([0.0]*3)):
        self.BondType = BondType
        self.BondLength = Dist
        self.BondVector = Vec

    def PrintBond(self):
        """
        Print the information about one atom-pair.
        """
        print(("Bond %s: %f" % (self.BondType, self.BondLength)))

    def HaveBond(self, BondRange, Tol):
        """
        Check if the bond is really bonded.
        """
        if self.BondLength < BondRange[self.BondType] + Tol:
            return True
        else:
            return False
