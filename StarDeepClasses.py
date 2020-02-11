# -*- coding: utf-8 -*-

"""
Created on Tue Dec 24 20:50:16 2019

@author: Red (Marco Rosso)
"""
from dataclasses import dataclass

from typing import List

import numpy as np
import pandas as pd


import scipy.integrate as integrate
from time import time
from StarDeepDataImporter import getStarYield
from random import random
from random import choice

# import numpy as np
TIMETICK = 1  # the tick per update call
yearPERTICK = 100000000  # the year per tick in the calculation, 1m is the default
cUnit = 1000000000  # the scale of 1 solar mass in the composition, default is 1b
cDP = 0  # the number of decimal places in the composition, default is 0

# The base element dictionary all the compositions will use
allowed_elements = pd.read_csv("elements_and_isotopes.csv")
_base_elements = list(allowed_elements["name"])
_metal_elements = list(allowed_elemens.query("is_metal")["name"])


class Composition(dict):
    """
    Composition, a class composed by a dictionary with a key Element(), 
    and a int value indicating his mass  (10**-6 solar masses).
    is iterable, and has functions to update the current mass, by adding or settings single or multiple element.
    Initializided with a the current periodic table + isotopes, that retrieve from the file "Elements&Isotopes.csv",
    to be done apropialy, the file must be a .csv, the separator must be a tabulation (the tab key), and the text must                      
    NOT be ""ended.
    in case, a new element can to be added.
    """

    def __init__(self, new_mass: dict = None):
        super.__init__(self)

        if new_mass:
            for element, mass in new_mass.items():
                self[element] = mass

    def __set_item__(self, key, value):
        if key not in _base_elements:
            raise KeyError(
                f"Element {key} not valid. See .allowed_elements for list of allowed elements."
            )
        if value < 0:
            raise ValueError("Mass can not be smaller than zero")

        super().__setitem__(key, value)

    def __add__(self, other):
        joint_elements = set(self.keys()).union(other.keys())
        new_composition = {}

        for element in joint_elements:
            new_composition[element] = self.get(element, 0) + other.get(element, 0)

        return Composition(new_composition)

    def convert(self, element1, element2, mass):
        """Convert a certain amount of a element 1 to element 2"""

        self[element1] += -mass
        self[element2] += mass

    def extract(self, elem, q):
        """
        extract:
        elem, the key element to extract from the , can be an array
        q, the quantity to extract,
        remove a quantity from the current composition and return it as
        a new composition
        """

        a = Composition()
        if type(elem) == list and type(q) == list and len(elem) == len(q):
            nC = Composition()
            for e, qa in zip(elem, q):
                e = self.getElKey(e)
                self.elements.addSingleMass(elem, -qa)
                nC.addSingleMass(elem, q)
            return nC
        elif type(elem) == list:
            raise KeyError("element list and quantity list must be of the same size")
        elem = self.getElKey(elem)
        for k in self.elements.keys():
            if elem == k:
                self.addSingleMass(elem, -q)
                return a.addSingleMass(elem, q)

    def genericEx(self, mass):
        """
        Generic extract an  amount of mass picking a percentile of all elements
        """
        exCom = Composition()
        a = self.totalMass()
        assert mass <= a, "tring to extract more mass that what there is"
        try:
            for k, v in self.elements.items():
                nv = round((v / a) * mass, cDP)
                self.addSingleMass(k, -nv)
                exCom.addSingleMass(k, nv)
        except ZeroDivisionError:
            return exCom
        return exCom

    def metal_mass(self):
        """Return the total metal mass in the composition"""

        metal_elements = [
            element for element in self.keys() if element in _metal_elements
        ]

        metal_mass = 0
        for metal_el in metal_elements:
            metal_mass += self[metal_el]

        return metal_mass

    def total_mass(self):
        """Return the total mass in the composition"""

        mass = 0
        for el, el_mass in self.items():
            mass += el_mass

        return mass


@dataclass
class Star:
    """
    Abstract Class Stars, will be used in all other star classes
    mass, a float defining the mass of the stars contained
    number, an int defining the number of similar star
    composition, a compositionClass defining the elements contained in the class
    sisters, a list of other star if it is a binary or more system
    lifespanDef, a string defining the equation for the calculate the star lifespan
        chose from "Default", "simple", "SimpleClassic"
    """

    sisters: List[Star] = []
    mass: float
    number: float
    lifespanDef: str = "Default"
    position: float
    ID: str

    def __init__(
        self, mass, number, composition, Position, sisters=[], lifespanDef="Default"
    ):
        self.sisters = sisters
        self.mass = mass
        self.number = number
        self.composition = composition
        self.lifespanDef = lifespanDef
        self.lifespan = self.iniLifeSpan()
        self.pos = Position
        self.ID = self.genID()
        pass

    def iniLifeSpan(self):
        """
        get the stellar lifespam (stellar time), based on the mass
        """
        m = self.mass
        if self.lifespanDef == "Default":
            if m <= 10:
                a = 1.2 * (10 ** 10) * (m ** -2.78)
                return int(a // yearPERTICK)
            else:
                a = 1.1 * (10 ** 8) * (m ** -0.75)
                return int(a // yearPERTICK)
        elif self.lifespanDef == "Simple":
            a = 11.7 * ((10 ^ 9) / (m ^ 2))
            return int(a // yearPERTICK)
        elif self.lifespanDef == "SimpleClassic":
            return 1

    def genID(self):
        if self.__class__.__name__ == "LowStars":
            o = "1"
        elif self.__class__.__name__ == "HighStars":
            o = "2"
        elif self.__class__.__name__ == "WDStars":
            o = "3"
        elif self.__class__.__name__ == "NeutronStars":
            o = "4"
        else:
            o = "5"
        try:
            if self.sisters[0].__class__.__name__ == "LowStars":
                t = "1"
            elif self.sisters[0].__class__.__name__ == "HighStars":
                t = "2"
            elif self.sisters[0].__class__.__name__ == "WDStars":
                t = "3"
            elif self.sisters[0].__class__.__name__ == "NeutronStars":
                t = "4"
            else:
                t = "5"
        except:
            t = "0"

        try:
            if self.sisters[1].__class__.__name__ == "LowStars":
                tr = "1"
            elif self.sisters[1].__class__.__name__ == "HighStars":
                tr = "2"
            elif self.sisters[1].__class__.__name__ == "WDStars":
                tr = "3"
            elif self.sisters[1].__class__.__name__ == "NeutronStars":
                tr = "4"
            else:
                tr = "5"
        except:
            tr = "0"

        pos = "{:06d}".format(self.pos)
        return "".join(o + t + tr + pos)

    def update(self):
        """
        Update the star lifespan
        return None if the star or her's sisters didn't die
        return a tuple of (yield, renmant) if the star died
        return a tuple of (yield, None) if just a sister died
        """
        empty = Composition()
        res = Composition()
        res2 = Composition()
        for s in self.sisters:
            yie = s.update()
            if yie != None:
                res += yie[0]
                if yie[1] != None:
                    self.sisters[self.sisters.index(s)] = yie[1]
        self.lifespan -= TIMETICK
        if self.lifespan <= 0:
            res2, ren = self.endStars()
            if res != empty:
                res2 += res
                return res2, ren
            return res2, ren
        if res != empty:
            return res, None
        return None, None

    def getNum(self):
        """
        return the number of star in the class
        """
        return self.number

    def getPos(self):
        """
        Return the index position in the star group array
        """
        return int(self.ID[3:])

    def setPos(self, pos):
        self.pos = pos
        self.ID = self.genID()

    def endStars(self):
        """
        Stars is an abstract class, use a child
        class for have any effect
        """
        pass


class LowStars(Stars):
    """
    Stars of mass between 0.08 and 8 solar masses,
    leave behid a whitedwarf as renmant
    """

    def endStars(self):
        """
        Return a tuple of the product of the star yield and the renmant
        star
        """
        sYield = getStarYield(self.composition.metalMassFract(), self.getMass())
        renmant = sYield.pop("rem")
        self.composition.genericEx(renmant * self.getNum() * cUnit)
        for k, v in sYield.items():
            if v < 0:
                try:
                    self.composition.extract(k, -v)
                    continue
                except ArithmeticError:
                    self.composition.setSingleMass(k, 0)
                    continue
            self.composition.convert("h", k, v * cUnit * self.getNum())

        product = self.composition
        self.composition = Composition()
        renmantComp = Composition()
        # it should never spawn a bHoles or a neutron star, but i leave this here for safety
        if renmant > 2.4:
            return (
                product,
                BHoles(
                    renmant,
                    self.getNum(),
                    renmantComp.addSingleMass("bm", renmant * self.getNum() * cUnit),
                    self.pos,
                    self.sisters,
                    self.lifespanDef,
                ),
            )
        elif renmant > 1.44:
            return (
                product,
                NeutronStars(
                    renmant,
                    self.getNum(),
                    renmantComp.addSingleMass("gm", renmant * self.getNum() * cUnit),
                    self.pos,
                    self.sisters,
                    self.lifespanDef,
                ),
            )
        return (
            product,
            WDStars(
                renmant,
                self.getNum(),
                renmantComp.addSingleMass("wm", renmant * self.getNum() * cUnit),
                self.pos,
                self.sisters,
                self.lifespanDef,
            ),
        )


class HighStars(Stars):
    """
    Stars of mass over the 8 solar masses,
    leave behid a neutron stars and black hole as renmant,
    the mainb difference is that they cause supernoave event
    """

    def endStars(self):
        """
        Return a tuple of the product of the star yield and the renmant
        star
        """
        sYield = getStarYield(self.composition.metalMassFract(), self.getMass())
        renmant = sYield.pop("rem")
        self.composition.genericEx(renmant * self.getNum() * cUnit)
        for k, v in sYield.items():
            if v < 0:
                try:
                    self.composition.extract(k, -v)
                    continue
                except ArithmeticError:
                    self.composition.setSingleMass(k, 0)
                    continue
            self.composition.convert("h", k, v * cUnit * self.getNum())
        product = self.composition
        self.composition = Composition()
        renmantComp = Composition()
        if renmant > 2.4:
            return (
                product,
                BHoles(
                    renmant,
                    self.getNum(),
                    renmantComp.addSingleMass("bm", renmant * self.getNum() * cUnit),
                    self.pos,
                    self.sisters,
                    self.lifespanDef,
                ),
            )
        elif renmant > 1.44:
            return (
                product,
                NeutronStars(
                    renmant,
                    self.getNum(),
                    renmantComp.addSingleMass("gm", renmant * self.getNum() * cUnit),
                    self.pos,
                    self.sisters,
                    self.lifespanDef,
                ),
            )
        return (
            product,
            WDStars(
                renmant,
                self.getNum(),
                renmantComp.addSingleMass("wm", renmant * self.getNum() * cUnit),
                self.pos,
                self.sisters,
                self.lifespanDef,
            ),
        )


class Renmants(Stars):
    """
    Abstract class for inerte renmant Stars
    in the update def there will
    be no reducing lifespan
    """

    def update(self):
        """
        return None if the star or her's sisters didn't die
        return a tuple of (yield, None) if just a sister died
        """
        empty = Composition()
        res = Composition()
        for s in self.sisters:
            yie = s.update()
            if yie != None:
                res += yie[0]
                if yie[1] != None:
                    self.sisters[self.sisters.index(s)] = yie[1]
        if res != empty:
            return res, None
        return None, None


class WDStars(Renmants):
    """
    White Dwarf Stars
    """


class NeutronStars(Renmants):
    """
    Neutron Stars
    """

    def SpinVeryFast(self):
        """
        It will spin Very Fast
        """
        pass


class BHoles(Renmants):
    """
    Black Holes
    """


class SimpleModelUnit(object):
    """
    Single Simple model:
    composition: a Composition class defining the initial "free" mass of the sistem
    name: a string defining the name (used to keep track of more than 1 model during massive computation)
    surface: an int defining the parsec squared surfaces
    SFR: stellar formation rate, defined in solar masses per mil. years per squared parsec
    age: the initial timeTICK age of the simulation
    pB: binary probability, a float between 0 and 1
    edges: link to the other simple model
    stars: a dictionary containig the initial star status
    IMF: a string defining the chosen IMF definition (chose between: "Salpeter", "Scalo", "Kroupa")  
    SMLmethod, define the solar masses range to integrate for the IMF
    sT, define the calculation for the star lifespam (chose between: "Default", "Dimple", "SimpleClassic")
    """

    def __init__(
        self,
        composition,
        name="SimpleModel",
        surface=1,
        SFR=1,
        bP=0.85,
        age=0,
        edges=None,
        stars=None,
        IMF="Salpeter",
        SMLmethod="default",
        sT="default",
    ):

        self.IMF = IMF.lower().capitalize()
        imfModels = ["Salpeter", "Scalo", "Kroupa"]
        if self.IMF not in imfModels:
            raise ValueError("imf models acepted: Salpeter, Scalo, Kroupa")

        self.name = name
        self.composition = composition
        self.surface = surface
        self.edges = []
        self.SMLmethod = SMLmethod
        self.age = age
        self.sT = sT.lower().capitalize()

        if edges != None:
            self.edges = edges

        self.stars = {}
        if stars != None:
            self.stars = stars
        self.totalMassStarBorn = 0
        if stars != None:
            for s in stars:
                self.totalMassStarBorn += s.getMass()

        self.SFR = SFR
        """gsm, generated star masses, normalize the SFR in the current situation.
        moltiplied by the current surface density of gas of the simple model"""
        self.gsm = round(
            (SFR * (surface * (yearPERTICK / 1000000000)))
            * self.getSurfaceMassDensity(),
            3,
        )

        """The percentage of stars born in a binary system"""
        self.bP = bP

        """bS, borned star (mass), collect the result of the stellar formation and keep storred
        the values not reaching the unit"""
        self.bS = [0] * len(self.getSML(self.SMLmethod))
        """bBS, borned Binary Star, collect the result of the stellar formation for the binary system"""
        self.bBS = [0] * len(self.getSML(self.SMLmethod))
        pass

    def totalMass(self):
        """
        Return the total mass of the model in solar masses
        """
        mass = self.composition.totalMass() / cUnit
        for g in self.stars:
            for s in g:
                mass += s.getMass()
        return mass

    def gasMass(self):
        """
        Return the gas mass of the model (the "free" mass)
        """
        return self.composition.totalMass()

    def fractMass(self):
        """
        fraction of mass in gas
        """
        return self.composition.totalMass() / self.totalMass()

    def globMetal(self):
        """
        global metallicity
        """
        met = 0
        met += self.composition.metalMass()
        for g in self.stars:
            for s in g:
                met += s.composition.metalMass
        return met / self.gasMass()

    def getSurface(self):
        """
        Get the surface (pc)**2
        """
        return self.surface

    def getTotalStarMassBorn(self):
        """
        return the current total mass of all born star
        """
        return self.totalMassStarBorn

    def starCount(self):
        """
        return a dictionari withe the number of stars for each category,
        by mass if still alive,
        by the kind of renmant if dead
        """

        def insertCount(s):
            if s.id()[0] == "3":
                try:
                    res["WhiteDwarfs"] += s.getNum()
                except:
                    res.update({"WhiteDwarfs": s.getNum()})
            elif s.id()[0] == "4":
                try:
                    res["NeutronStars"] += s.getNum()
                except:
                    res.update({"NeutronStars": s.getNum()})
            elif s.id()[0] == "5":
                try:
                    res["BlackHoles"] += s.getNum()
                except:
                    res.update({"BlackHoles": s.getNum()})
            else:
                try:
                    res[s.getMass()] += s.getNum()
                except:
                    res.update({s.getMass(): s.getNum()})

        res = {}
        tC = 0
        for k, g in self.stars.items():
            for s in g:
                for sis in s.getSisters():
                    insertCount(sis)
                    tC += sis.getNum()
                insertCount(s)
                tC += s.getNum()
        res.update({"tC": tC})
        return res

    def printStars(self, full=False):
        """
        Print the current star population
        """
        res = self.starCount()
        if full == True:
            for k, v in res.items():
                if type(k) != str:
                    print(str(k) + "M: " + "{:,}".format(v))
                elif k == "WhiteDwarfs":
                    print("WhiteDwarfs: " + "{:,}".format(v))
                elif k == "NeutronStars":
                    print("NeutronStars: " + "{:,}".format(v))
                elif k == "BlackHoles":
                    print("BlackHoles: " + "{:,}".format(v))
                elif k == "tC":
                    print("TotalCount: " + "{:,}".format(v))
        else:
            for k, v in res.items():
                if type(k) != str:
                    print(str(k) + "M: " + "{:,}".format(v))
                elif k == "tC":
                    print("TotalCount: " + "{:,}".format(v))

    def IMFdef(self, SM):
        """
        return the Initial Mass Function for each solar mass range,
        depending by the definition
        """
        if self.IMF == "Salpeter":
            return 1.35
        if self.IMF == "Scalo":
            if SM <= 2:
                return 1.35
            if SM > 2:
                return 1.7
        if self.IMF == "Kroupa":
            if SM <= 0.5:
                return 0.3
            if SM > 0.5 and SM <= 1:
                return 1.2
            if SM > 1:
                return 1.7

    def getSurfaceMassDensity(self):
        """
        surface mass density of gas (solar masses)
        """
        return self.gasMass() / cUnit / self.surface

    def getSML(self, method="default"):
        """
        return the solar masses list in consideration, the format is: ([1]:
        """
        if method == "default":
            return [
                (0.08, 0.1, 0.12),
                (0.12, 0.14, 0.16),
                (0.16, 0.18, 0.20),
                (0.20, 0.24, 0.28),
                (0.28, 0.32, 0.36),
                (0.36, 0.40, 0.44),
                (0.44, 0.52, 0.60),
                (0.60, 0.68, 0.76),
                (0.76, 0.84, 0.92),
                (0.92, 1.08, 1.24),
                (1.24, 1.40, 1.56,),
                (1.56, 1.72, 1.88),
                (1.88, 1, 94, 2.00),
                (2.00, 2.25, 2.50),
                (2.50, 2.75, 3.00),
                (3.00, 3.25, 3.50),
                (3.50, 3.75, 4.00),
                (4.00, 4.50, 5.00),
                (5.00, 5.50, 6.00),
                (6.00, 6.50, 7.00),
                (6.00, 7.00, 8.00),
                (8.00, 9.00, 10.00),
                (10.00, 10.00, 12.00),
                (12.00, 14.00, 16.00),
                (16.00, 18.00, 20.00),
                (20.00, 22.00, 24.00),
                (24.00, 28.00, 32.00),
                (32.00, 36.00, 40.00),
                (40.00, 50.00, 60.00),
                (60.00, 70.00, 80.00),
                (80.00, 90.00, 100.00),
                (100.00, 110.00, 120.00),
            ]  # debug test
        pass

    def getInitialMass(self):
        """
        get initial masses:
        it will return a list with the number of star in each mass category
        """
        ms = self.getSML()
        pr = []
        for m in ms:
            # pr.append(1/(m**(1+self.IMFdef(m))))
            pr.append(
                float(
                    integrate.quadrature(
                        lambda x: 1 / (x ** (1 + self.IMFdef(x))), m[0], m[2]
                    )[0]
                )
            )
        a = sum(pr)
        for p, b in zip(pr, self.bS):
            self.bS[self.bS.index(b)] += (p / a) * self.gsm
        return self.bS

    def stellarGen(self):
        """
        Generate the single Stars. Should always be called before the binary generation
        """
        for s, m in zip(self.getInitialMass(), self.getSML(self.SMLmethod)):
            if s / m[1] >= 1:
                c = int(s // m[1])
                si = self.bS.index(s)
                self.bS[si] = s % m[1]
                self.bBS[si] = int(c * self.bP)
                c = c - int(c * self.bP)
                try:
                    a = self.composition.genericEx(round(c * m[1] * cUnit, cDP))
                except AssertionError:
                    c *= self.composition.totalMass() / (c * m[1] * cUnit)
                    c = int(c // 1)
                    a = self.composition.genericEx(round(c * m[1] * cUnit, cDP))
                if m[1] < 8:
                    try:
                        p = len(self.stars[m[1]]) - 1
                        self.stars[m[1]].append(
                            LowStars(m[1], c, a, p, lifespanDef=self.sT)
                        )
                    except:
                        self.stars.update(
                            {m[1]: [LowStars(m[1], c, a, 0, lifespanDef=self.sT)]}
                        )
                else:
                    try:
                        p = len(self.stars[m[1]]) - 1
                        self.stars[m[1]].append(
                            HighStars(m[1], c, a, p, lifespanDef=self.sT)
                        )
                    except:
                        self.stars.update(
                            {m[1]: [HighStars(m[1], c, a, 0, lifespanDef=self.sT)]}
                        )
                self.totalMassStarBorn += m[1] * c
        self.gsm = round(
            (self.SFR * (self.surface * (yearPERTICK / 1000000000)))
            * self.getSurfaceMassDensity(),
            3,
        )
        pass

    def bStellarGen(self, bBS, SML):
        """
        Generate the binary stars
        """
        c = 0
        while sum(bBS) > 1:
            a = choice(range(len(bBS)))
            if bBS[-1] > 0:
                x = SML[-1][1]
                y = SML[a][1]
                c = int(bBS[-1] * random()) if int(bBS[-1] * random()) > 1 else 1
                if bBS[a] < c and bBS[a] != bBS[-1]:
                    c = bBS[a]
                elif bBS[a] == bBS[-1]:
                    c = int(bBS[a] // 2)
                if c == 0:
                    continue
                bBS[a] -= c
                bBS[-1] -= c
                xComp = self.composition.genericEx(round(c * x * cUnit, cDP))
                yComp = self.composition.genericEx(round(c * y * cUnit, cDP))

                # Generate the minor stars
                if y < 8:
                    m = LowStars(y, c, yComp, 0, lifespanDef=self.sT)
                else:
                    m = HighStars(y, c, yComp, 0, lifespanDef=self.sT)

                # Generate the Major Stars
                p = len(self.stars[x]) - 1
                if x < 8:
                    M = LowStars(x, c, xComp, p, sisters=[m], lifespanDef=self.sT)
                else:
                    M = HighStars(x, c, xComp, p, sisters=[m], lifespanDef=self.sT)

                # Store the result in the stars list
                try:
                    self.stars[M.getMass()].append(M)
                except:
                    self.stars.update({M.getMass(): [M]})
            else:
                self.bStellarGen(bBS[:-1], SML[:-1])
                break
        pass

    def starsYield(self, g, v):
        """
        Calculate the stellar yield of a list of stars
        v, a list of stars
        g, the indices in which the stars are located
        """
        for s in v:
            a, r = s.update()
            if a != None:
                self.composition += a
            if r != None:
                self.stars[g].pop(v.index(s))
                try:
                    r.setPos(len(self.stars["Renmant"]) - 1)
                    self.stars["Renmant"].append(r)
                except:
                    r.setPos(0)
                    self.stars.update({"Renmant": [r]})
            self.gsm = round(
                (self.SFR * (self.surface * (yearPERTICK / 1000000000)))
                * self.getSurfaceMassDensity(),
                3,
            )

    def update(self):
        """
        Update the current simple unit composition
        firstly, elaborate the yield of the current star population, both
        composition and stellar renmant creation
        secondly, calculate the stars born in this unitTime
        do not return anithing
        """
        self.age += 1

        for g, v in list(self.stars.items()):
            self.starsYield(g, v)

        # Generate the Stars
        self.stellarGen()
        self.bStellarGen(self.bBS, self.getSML(self.SMLmethod))
        pass

    def Run(self, times):
        for t in range(times):
            start = time()
            self.update()
            print(
                str(round(t / times * 100, 2))
                + "% - "
                + str(round(time() - start, 3))
                + "s"
            )
        self.printStars(full=True)
        print()
        pass


a = Composition()
a.addSingleMass("H", 100000000000 * cUnit)
fi = SimpleModelUnit(a, surface=75398, SFR=0.1)
fi.Run(5)