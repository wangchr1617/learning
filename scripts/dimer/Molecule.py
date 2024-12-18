#!/usr/bin/env python
"""
A class to describe a molecule built by atoms.
A molecule is kept fixed during the basin hopping process.
No bond breaking is allowed.
Class layers: Element --> Atom --> Molecule --> MoleculeSet.
"""
import random
from Atom import *

class Molecule(object):
    """
    Define a molecule made by atoms.
    """
    def __init__(self):
        """
        Create an empty molecule.
        """
        self.NumofAtoms = 0
        self.Atoms = []
        # The adsorb type is determined by the property of atoms.
        # If there are any frozen atoms in the molecule, it is also Frozen.
        self.Type = ""
        self.Center = np.array([0.0]*3)
        self.CenterOfMass = np.array([0.0]*3)
        self.MolName = ""
        self.Force = np.array([0.0]*3)
        self.TotalForce = np.array([0.0]*3)

    def PrintMol(self):
        """
        Simply print the molecule and all its coordinates in XYZ format.
        This function is for debugging.
        """
        print(("%s" % (self.MolName)))
        for myatom in self.Atoms:
            line = "%s %.6f %.6f %.6f %s " % (myatom.Symbol, myatom.Coord[0], myatom.Coord[1], myatom.Coord[2], myatom.Type)
            print(line)

    def SetMolType(self):
        """
        Set the molecule Type property according to the AdsorbType of atoms.
        The whole molcule should be marked as fixed if one atom in this molecule is marked as F.
        """
        for myatom in self.Atoms:
            if myatom.Type == "F":
                self.Type = "F"
                return
        self.Type = ""

    def MoveMol(self, shiftvector):
        """
        Move the coordinates of molecule by a shift vector.
        """
        for i in range(0, self.NumofAtoms):
            self.Atoms[i].Coord += shiftvector

    def GetCenter(self):
        """
        Calculate the shape center of a molecule.
        """
        center = np.array([0.00]*3)
        if self.NumofAtoms == 0:
            print("Warning: No atoms in this molecule!\n")
        else:
            for i in range(0, self.NumofAtoms):
                center += self.Atoms[i].Coord
            center /= self.NumofAtoms
        return center

    def SetCenter(self):
        """
        Set the center of the molecule.
        """
        self.Center = self.GetCenter()

    def GetCenterOfMass(self):
        """
        Get the center of Mass of this Molecule
        """
        center = np.array([0.0]*3)
        totalmass = 0
        if self.NumofAtoms == 0:
            print("Warning: No atoms in this molecule!\n")
        else:
            for i in range(0, self.NumofAtoms):
                center += self.Atoms[i].Coord*self.Atoms[i].AtomicMass
                totalmass += self.Atoms[i].AtomicMass
        center /= totalmass
        return center

    def SetCenterOfMass(self):
        """
        Set the CenterOfMass of the molecule.
        """
        self.CenterOfMass = self.GetCenterOfMass()

    def CenterMol(self):
        """
        Move the shape center of a molecule to (0,0,0).
        """
        shiftvector = np.array([0.0]*3)
        shiftvector -= self.GetCenter()
        self.MoveMol(shiftvector)

    def GetShiftVector(self, BHStep, ShiftAngleRange, MoleculeSetCenter, IsPlanarSearch):
        """
        Get a shift vector with step and the allowed angle range.
        This function is never used. So I'd better try to use it, or just ignore it.
        Do it later.
        """
        ShiftValue = np.array([0.0]*3)
        self.SetCenter()
        Vec2Center = MoleculeSetCenter - self.Center
        #debug
        #print "Vec2Center is %s" % (Vec2Center)
        RandSize = random.uniform(0, BHStep)
        RandAngle = random.uniform(ShiftAngleRange[0] - 90, ShiftAngleRange[1] - 90)
        if IsPlanarSearch:
            ShiftValue = Normalize(RandomPerpendicularVectorXY(Vec2Center) - Normalize(Vec2Center)*math.tan(math.radians(RandAngle)))*RandSize
        else:
            ShiftValue = Normalize(RandomPerpendicularVector(Vec2Center) - Normalize(Vec2Center)*math.tan(math.radians(RandAngle)))*RandSize
        return ShiftValue

    def HaveBond(self, atomnum1, atomnum2, BondRange, tol):
        """
        Check if atomnum1 and atomnum2 in this molecule are bonded.
        When the distance of atomnum1 and atomnum2 is smaller than the bondrange + tol, they have bond.
        """
        atomnum1 -= 1
        atomnum2 -= 1
        mybond = self.Atoms[atomnum1].GetBond(self.Atoms[atomnum2])
        if mybond.HaveBond(BondRange, tol):
            print(("Distance of Atom %d and %d is %f, shorter than %f. They're bonded." % (atomnum1, atomnum2, mybond.BondLength, BondRange[mybond.BondType] + tol)))
            return True
        else:
            return False

    def GetNearestAtomPair(self, mol2):
        """
        Get the nearest atom-pair from this molecule and the other molecule.
        Stupid algorithm, but hopefully it works.
        """
        atomnum1 = 0
        atomnum2 = 0
        mybond = self.Atoms[atomnum1].GetBond(mol2.Atoms[atomnum2])
        dist = mybond.BondLength
        for i in range(0, self.NumofAtoms):
            for j in range(0, mol2.NumofAtoms):
                newbond = self.Atoms[i].GetBond(mol2.Atoms[j])
                if newbond.BondLength < dist:
                    dist = newbond.BondLength
                    atomnum1 = i
                    atomnum2 = j
        #debug
        #print "Distance in atom %d %s in this molecule and atom %d %s in mol2  is %.6f." % (atomnum1, self.Atoms[atomnum1].Symbol, atomnum2, mol2.Atoms[atomnum2].Symbol, dist)
        pair = [atomnum1, atomnum2]
        return pair

    def CountAtoms(self):
        """
        A function to count the number of atoms of each kind.
        """
        mydict = {}
        for myatom in self.Atoms:
            symbol = myatom.Symbol
            if not symbol in mydict:
                mydict[symbol] = 1
            else:
                mydict[symbol] += 1
        return mydict

    def SetMolName(self):
        """
        This function will count the number of atoms of each kind of elements.
        Then a name is generated from the number of atoms.
        eg: CH3OH should be named as H4CO.
        """
        mydict = self.CountAtoms()
        molname = ""
        for element in list(mydict.keys()):
            if 1 == mydict[element]:
                molname += element
            else:
                molname += element + str(mydict[element])
        self.MolName = molname

    def SameMolName(self, mol):
        """
        a simple function to compare two the name of two molecules.
        If the name of the molecules are the same, they are considered to be the same.
        """
        self.SetMolName()
        mol.SetMolName()
        if self.MolName == mol.MolName:
            return True
        else:
            return False

    def SplitMol(self, atomnum1, atomnum2, BondRange, tol):
        """
        Split a chemical bond.
        Algorithm Description:
        For atom n and atom m,
        Find all atoms connectted to atom n and atoms m,
        then move the molecule along the bond vector in opposite direction.
        """
        if not self.HaveBond(atomnum1, atomnum2, BondRange, tol):
            print(("Atom N.O. %d and %d are not bonded." % (atomnum1, atomnum2)))
            return
        #First, split the molecule into two parts:
        #Part 1 connectted with atomnum1 and part2 connected with atomnum2
        full = set(range(0, self.NumofAtoms))
        set1 = set([atomnum1])
        set2 = set([atomnum2])
        part1 = set()
        part1 = part1.union(set1)
        part2 = set()
        part2 = part2.union(set2)
        rest = full - part1 - part2
        neighbors = set()
        #Get all atoms connected with atomnum1
        while True:
            for i in set1:
                for j in rest:
                    if self.HaveBond(i, j, BondRange, tol):
                        neighbors.add(j)
            if 0 == len(neighbors):
                break
            else:
                part1 = part1.union(neighbors)
                set1 = neighbors
                neighbors = set()
                rest = full - part1 - part2
        #Get all atoms connected with atomnum2
        while True:
            for i in set2:
                for j in rest:
                    if self.HaveBond(i, j, BondRange, tol):
                        neighbors.add(j)
            if 0 == len(neighbors):
                break
            else:
                part2 = part2.union(neighbors)
                set2 = neighbors
                neighbors = set()
                rest = full - part1 - part2
        if 0 == len(rest):
            print(("Molecule can be divided to two parts by splitting bond between atom %d and %d " % (atomnum1, atomnum2)))
        else:
            print("Warning: Molecule is more than two parts")

    def SetTotalForce(self):
        """
        Set the total force of the whole molecule
        """
        for i in range(0, self.NumofAtoms):
            self.TotalForce += self.Atoms[i].Force
        return True

    def RotateMol(self, axis, angle):
        """
        Rotate a molecule along given axis by certain angle
        """
        #First, put the molecule center to (0,0,0)
        center = self.GetCenter()
        self.CenterMol()
        #Then rotate the molecule
        R = GetRotationMatrix(axis, angle)
        for i in range(0, self.NumofAtoms):
            self.Atoms[i].Coord = RotateByMatrix(self.Atoms[i].Coord, R)
        #Finally, move the molecule back to the center
        self.MoveMol(center)
