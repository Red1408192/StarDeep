# -*- coding: utf-8 -*-

"""
Created on Tue Dec 24 20:50:16 2019

@author: Red (Marco Rosso)
"""
import scipy.integrate as integrate
from StarDeepDataImporter import elIso
from StarDeepDataImporter import getStarYield
from random import random
from random import choice
# import numpy as np
TIMETICK = 1 # the tick per update call
yearPERTICK = 1000000 #the year per tick in the calculation, 1m is the default
cUnit = 1000000000 #the scale of 1 solar mass in the composition, default is 1b
cDP = 0 #the number of decimal places in the composition, default is 0



class Element():
    '''
    Element: A class that store all elements and Isotopes data and names,
    Usualy you shall see just the english name, symbol and massID from the next class,
    but you can retrieve many other attribuites
    '''
    def __init__(self, n, sy= None, aN = None, gN = '', fN = '', sN = '', iN = '',
                 aM = '', me = True, headers = None):
        self.rAV = 999
        if sy == None:
            sy = n[::2]
        if aN == None:
            aN = self.rAV
            self.rAV -= 1
        
        self.keys = [n.lower().capitalize(),
                     sy.lower().capitalize(),
                     str(aN).lower().capitalize(),
                     gN.lower().capitalize(),
                     fN.lower().capitalize(),
                     sN.lower().capitalize(),
                     iN.lower().capitalize()]
        self.data = [aM, True if me == "True" else False]
        if headers == None or len(headers) < len(self.keys)+len(self.data):
            headers = ["English Name", "Symbol", "Atomic Number"] + ["" for c in range(len(self.keys)+len(self.data))]   
        self.header = headers
        pass
    
    def __contains__(self, k):
        k = k.lower().capitalize()
        if k in self.keys:
            return True
        return False
    
    def __eq__(self, other):
        if self.__class__ == other.__class__:
            if self.keys == other.keys:
                return True
            return False
        if type(other) == str:
            if other.lower().capitalize() in self.keys:
                return True
            return False
        raise NotImplementedError

    def __ne__(self, other):
        if self.__class__ == other.__class__:
            if self.keys != other.keys:
                return True
            return False
        if type(other) == str:
            if not other.lower().capitalize() in self.keys:
                return True
            return False
        raise NotImplementedError
    
    def __hash__(self): #doesn't work as intended
        return hash((self.keys[0], self.keys[1], self.keys[2], self.keys[3], self.keys[4], self.keys[6]))
    
    def __str__(self):
        return " ".join([self.keys[0], self.keys[1], self.keys[2]])
    
    def getAM(self):
        '''
        Get Atomic Mass
        '''
        return self.data[0]

    def isMetal(self):
        '''
        Return True if element is metallic (astrophisics)
        '''
        return self.data[1]
    def ShowKeys(self):
        '''
        Show the keys on which the element is defined
        '''
        return self.keys
    
    def ShowAll(self):
        '''
        Show all the element attribuites
        '''
        return "".join(h + " : " + a + "\n" for a, h in zip(self.keys+self.data, self.header))


class Composition():
    '''
    Composition, a class composed by a dictionary with a key Element(), 
    and a int value indicating his mass  (10**-6 solar masses).
    is iterable, and has functions to update the current mass, by adding or settings single or multiple element.
    Initializided with a the current periodic table + isotopes, that retrieve from the file "Elements&Isotopes.csv",
    to be done apropialy, the file must be a .csv, the separator must be a tabulation (the tab key), and the text must                      
    NOT be ""ended.
    in case, a new element can to be added.
    '''
    def __init__(self, elements = None): #never set as a default atribuite a list or a dictionary, they will have the same location between all the class instances
        self.index = 0
        self.elements = {}
        for i in range(len(elIso)): #Collect the data from the elIso global variable
            b = Element(elIso[i][0], elIso[i][1], elIso[i][2],
                        elIso[i][3], elIso[i][4], elIso[i][5],
                        elIso[i][6], elIso[i][7], elIso[i][8])
            self.elements.update({b : 0})
        
        if elements != None:
            self.aM(elements)
        
        self.rAV = 999
        pass
        
        
    # def nI(self, name):
    #     '''
    #     normalize the input
    #     '''
    #     if type(name) == int:
    #         name = str(name)
    #     elif type(name) == tuple:
    #         return name
    #     for k in self.elements.keys():
    #         if name.lower().capitalize() in k:
    #             return k
    #     name = (name, name[::2], self.rAV)
    #     self.rAV -= 1
    #     return name

    def aMS(self, elem, mass):
        """
        addMassSingle:
        elem, an atribuite defining the element to add (or subtract) mass on, could be:
            the name (English, German, French, Spanish, Italian work fine),
            the simbol,
            the Atomic value
            or another Element Class
            IF AN ISOTOPE, add underscore+the isotope massID (E.g. Helium_3, He_3, 2_3; check the csv file)
        mass, a real number defining the mass to add(or subtract), it will never leave the mass go below zero
        
        if the elements is not present, and defined with an element class, it will add as a new one in the dictionary
        """
        # elem = self.nI(elem)
        elem = self.gEl(elem)
        try:
            if elem in self.elements.keys() and self.elements[elem] + mass >= 0:
                self.elements.update({elem : round(self.elements[elem] + mass, cDP)})
            elif elem in self.elements.keys():
                raise ArithmeticError("element ", elem, " has gone below zero mass, has been setted to zero")
            else:
                raise KeyError("Periodic Table aMS() failed to find the element ", elem)
        except KeyError:
            if mass > 0 and type(elem)==Element:
                self.elements.update({elem:mass})
                print("it has been added as new Element")
            elif type(elem)==Element:
                self.elements.update({elem:0})
                print("it has been added as new Element with a mass of 0 (because negative)")
            pass
        except ArithmeticError:
            self.elements.update({elem : 0})
            pass
        pass
    
    def aM(self, elems, masses = None):
        """
        AddMasses:
        cases
        1- elems is a LIST: list for the masses to add or subtract,
        in this case masses need to be a list of equal lenght
            in the exact order of the previus list
        2- elems is a DICTIONARY: compile it with the keys defining the elements to
        elaborate, and the values as the masses to add or subtracts (leave the masses field empty)
        3- elems is a Composition() class: it will merge the two composition
        """
        if type(elems)==list and masses!=None and len(masses)==len(elems):
            for e in elems:
                e = self.gEl(e)
                try:
                    self.aMS(e, masses[elems.index(e)])
                except KeyError:
                    print("Periodic Table AMS() failed to find the element ", e)
                    # if masses[elems.index(e)] > 0:
                    #     self.elements.update({e:masses[elems.index(e)]})
                    #     self.metalEl.append(True)
                    #     print("it has been added as a metallic new key ")
                    pass
        elif type(elems)==list and len(masses)!=len(elems):
            print("periodic Table aM() error: the two list doesn't match, operation aborted")

        elif type(elems)==dict:
            for k,v in elems.items():
                k = self.gEl(k)
                try:
                    self.aMS(k, v)
                except KeyError:
                    pass

        elif type(elems)==Composition:
            for e,v in elems.elements.items():
                try:
                    self.aMS(e, v)
                except KeyError:
                    pass
        pass
    
    def sMS(self, elem, mass):
        """
        SetMassSingle:
        elem, an atribuite defining the element to set mass on, could be:
            the name (English, German, French, Spanish, Italian work fine),
            the simbol,
            the Atomic value
            or another Element Class
            IF AN ISOTOPE, add underscore+the isotope massID (E.g. Helium_3, He_3, 2_3; check the csv file)
        mass, a real number defining the new mass, it will never leave the mass go below zero
        
        if the elements is not present, and defined with an element class, it will add as a new one in the dictionary
        """
        elem = self.gEl(elem)
        try:
            if elem in self.elements.keys() and mass >= 0:
                self.elements.update({elem : round(mass, cDP)})
            elif elem in self.elements.keys():
                raise ArithmeticError("element ", elem, " has gone below zero mass, has been setted to zero") 
                #print("element ", elem, " has gone below zero mass, has been setted to zero")
            else:
                raise KeyError
        except KeyError:
            print("Periodic Table aMS() failed to find the element ", elem)
            if mass > 0 and type(elem)==Element:
                self.elements.update({elem:mass})
                print("it has been added as new Element")
            elif type(elem)==Element:
                self.elements.update({elem:0})
                print("it has been added as new Element with a mass of 0 (because negative)")
            pass
        except ArithmeticError:
            self.elements.update({elem : 0})
            pass
    pass

    def sM(self, elems, masses = None):
        """
        SetMass:
        cases
        1- elems is a LIST: add another list for the masses to set in the exact order of the previus list
        2- elems is a DICTIONARY: compile it with the keys defining the elements to elaborate,
        and the values as the masses set (leave the masses field empty)
        """
        if type(elems)==list and masses!=None and len(masses)==len(elems):
            for e in elems:
                e = self.gEl(e)
                try:
                    self.sMS(e, masses[elems.index(e)])
                except KeyError:
                    print("Periodic Table AMS() failed to find the element ", e)
                    # if masses[elems.index(e)] > 0:
                    #     self.elements.update({e:masses[elems.index(e)]})
                    #     print("it has been added as new key")
        elif type(elems)==list and len(masses)!=len(elems):
            print("periodic Table aM() error: the two list don't have the same lenght, operation aborted")
        
        elif type(elems)==dict:
            for k,v in elems.items():
                k = self.gEl(k)
                try:
                    self.sMS(k, v)
                except KeyError:
                    pass

        elif type(elems)==Composition:
            for e,v in elems.elements.items():
                try:
                    self.sMS(e, v)
                except KeyError:
                    pass
        pass
    
    def gM(self, elem):
        """
        Get the mass of a single element
        """
        elem = self.gEl(elem)
        for e in self.elements.keys():
            if elem == e:
                return self.elements[e]
        raise KeyError
    
    def addElement(self, n, sy= None, aN = None, gN = '', fN = '', sN = '',
                   iN = '', aM = '', me = True, headers = None):
        '''
        A new element can be created, it will be instantieted only in this class,
        it requires only the english name, all other attrbuites can be ignored
        not recomended, YOU SHOULD ALWAIS MODIFY THE CSV FILE
        '''
        if sy == None:
            sy = n[::2]
        if aN == None:
            aN = self.rAV
            self.rAV -= 1
        if headers == None or len(headers) < len(self.keys)+len(self.data):
            headers = ["English Name", "Symbol", "Atomic Number"] + ["" for c in range(len(self.keys)+len(self.data))]   
        header = headers
        self.elements.update({Element(n, sy, aN, gN, fN, sN,
                                      iN, aM, me, header):0})
        pass
        
    def mergeE(self, key1, key2):
        '''
        mergeElements:
        define an element 1 (key 1), to collect the mass of a element 2 (key2), the last will be removed
        '''
        key1 = self.gEl(key1)
        key2 = self.gEl(key2)
        self.aMS(key1, self.gM(key2))
        self.elements.pop(key2)
        pass
        
    def splitE(self, key1, key2, weight=1):
        """
        SplitElement:
        the key 1 will be splittet by a weight between 0 and 1 in a key2 depositorary, 
        """
        assert weight <= 1, "weight must fall between 0 and 1"
        key1 = self.gEl(key1)
        key2 = self.gEl(key2)
        self.aMS(key2, self.gM(key1)*(weight))
        self.aMS(key1, -self.gM(key1)*(weight))
        pass
    
    def con(self, key1, key2, mass):
        """
        convert a certain amount of a key to another
        """
        key1 = self.gEl(key1)
        key2 = self.gEl(key2)
        self.aMS(key2, mass)
        self.aMS(key1, -mass)
        pass
        
    
    def gEl(self, key):
        '''
        Normalize the input
        '''
        if type(key)==str:
            key=key.lower().capitalize()
            for el in self.elements.keys():
                if key in el:
                    return el
        return key
    
    def elInComp(self, elem):
        """
        Return true if the element is in the composition
        """
        elem = self.gEl(elem)
        for k in self.elements.keys():
            return True
        return False
    
    def readElementKey(self, elem, full=False):
        """
        Read the element key, set sull to True for show all data about the key
        """
        elem = self.gEl(elem)
        for k in self.elements.keys():
            if elem == k:
                if full:
                    return k.ShowAll()
                return k
        raise KeyError("element not present")

    def ex(self, elem, q):
        """
        extract:
        elem, the key element to extract from the , can be an array
        q, the quantity to extract,
        remove a quantity from the current composition and return it as
        a new composition
        """
        a = Composition()
        if type(elem) == list and type(q) == list and len(elem)==len(q):
            nC = Composition()
            for e, qa in zip(elem, q):
                e = self.gEl(e)
                self.elements.aMS(elem, -qa)
                nC.aMS(elem, q)
            return nC
        elif type(elem) == list:
            raise KeyError("element list and quantity list must be of the same size")
        elem = self.gEl(elem)
        for k in self.elements.keys():
            if elem == k:
                self.aMS(elem, -q)
                return a.aMS(elem, q)
    
    def gEx(self, mass):
        """
        Generic extract an  amount of mass picking a percentile of all elements
        """
        exCom = Composition()
        a = self.totalMass()
        assert mass <= a, "tring to extract more mass that what there is"
        try:
            for k, v in self.elements.items():
                nv = round((v/a)*mass, cDP)
                self.aMS(k, -nv)
                exCom.aMS(k, nv)
        except ZeroDivisionError:
            return exCom
        return exCom

    def metalMass(self):
        """
        Return the total metal mass in the composition
        """
        mM = 0
        for v, b in zip(self.elements.values(), self.elements.keys()):
            if b.isMetal():
                mM += v
        return mM
    
    def nonMetalMass(self):
        """
        Return the total non metal mass in the composition
        """
        nMM = 0
        for v, b in zip(self.elements.values(), self.elements.keys()):
            if not b.isMetal():
                nMM += v
        return nMM
    
    def mMF(self):
        """
        Return the metal mass fraction of the composition
        """
        mM = self.metalMass()
        tM = self.totalMass()
        return round(mM/tM, 6)

    def totalMass(self):
        """
        Return the total mass in the composition
        """
        tM = 0
        for v in self.elements.values():
            tM += v
        return tM
    
    def compositionList(self):
        """
        Return the full list of the element in the composition
        """
        return [str(k) for k in self.elements.keys()]
    
    def __iter__(self):
        return self
    
    def __next__(self):
        """
        every iterable is a list of two elements, the SIMBOL and the MASS
        """
        self.index += 1
        if self.index < len(self.elements)+1:
            return [list(self.elements.items())[self.index-1][0], list(self.elements.items())[self.index-1][1]]
        raise StopIteration
        
    def __add__(self, other):
        for k, v in other.elements.items():
            self.aMS(k, v)
        return Composition(self.elements)
            
    def __sub__(self, other):
        for k, v in other.elements.items():
            self.aMS(k, -v)
        return Composition(self.elements)
        
    def __str__(self):
        return "".join(str(k) + ": " + str(v) + "\n" for k, v in self.elements.items() if v != 0)


class Stars():
    """
    Abstract Class Stars, will be used in all other star classes,
    create a list where similar star will be stored, when a star has the same lifespam
    and the same mass, instead of add another element, increase the number.
    ade, a list where the stargroup go to die.
    """
    def __init__(self, initial=None):
              
        if initial == None:
            self.SG=[]
        else:
            self.SG = initial
            
        self.defLife = 10
        pass
    
    def fSis(self, mass, lifespan=None):
        """
        Find Sister, return a int defining the index where a similar stargroup is in self.Stars
        """
        lifespan = self.defLife if lifespan == None else lifespan
        for s in self.SG:
            if s[0] == mass and s[2] == lifespan:
                return self.SG.index(s)
        raise KeyError
        
    def addStars(self, Comp, mass, number, lifespan=None):
        """
        addStars: it will try to find a sister with the same int mass, the same
        int lifespam,and if there is one, it will increase the number in that group
        Otherwise, create a new group.
        
        By default, addStars should receive mass from the update call by an equal ammount of mass removed from the
        system
        """
        lifespan = self.defLife if lifespan == None else lifespan
        # try:
        #     self.SG[self.fSis(mass, lifespan)][1] += number
        #     self.SG[self.fSis(mass, lifespan)][3] += Comp
        # except KeyError:
        #     self.SG.append([mass, number, lifespan, Comp])
        # pass
        self.SG.append([mass, number, lifespan, Comp])
    
    def getTotComp(self):
        """
        Return the total composition of the the same star class
        """
        res = Composition()
        for s in self.SG:
            res += s[3]
        return self.composition
    
    def getMass(self):
        """
        Get the total mass of the star class, in solar masses
        """
        mass = 0
        for s in self.SG:
            mass += s[0]*s[1]
        return mass
    
    def getPop(self):
        """
        Return the total number of star of the same star class
        """
        number = 0
        for s in self.SG:
            number += s[1]
        return number
    
    def getAverageLS(self):
        """
        Get the average lifespan of all star in the same star class
        """
        lS = 0
        for s in self.SG:
            lS += [2]*[1]
        return lS/self.getPop()

    def endStars(self, starGroup):
        # self.ade.append(self.SG.pop(self.SG.index(starGroup)))
        pass

    def update(self):
        """
        Update every star group in the class, reducing by a Time Tick the life of each of them,
        return a tuple of the product of the died star as a composition and the mass of the renmant
        """
        rn = []
        upP = Composition()
        if len(self.SG) == 0:
            return upP, rn
        for s in self.SG.copy():
            s[2] -= TIMETICK
            if s[2] <= 0:
                p, r = self.endStars(s)
                rn.append(r)
                upP += p
        return upP, rn
    
    def __str__(self):
        # res = []
        # for s in self.SG:
        #         res.append("Mass: " + str(s[0]) + ", Pop: " + str(s[1]) + ", Life: " + str(s[2]) + "\n")
        # a = "".join(res)
        res = {}
        for s in self.SG:
            if s[0] not in res:
                res.update({s[0]:s[1]})
                continue
            res[s[0]] += s[1]
        nRes = [str(k) + "M = " + str(v) + "\n" for k,v in res.items()]
        nRes = "".join(nRes)
        return nRes
    
    def printStars(self):
        res = {}
        for s in self.SG:
            if s[0] not in res:
                res.update({s[0]:s[1]})
                continue
            res[s[0]] += s[1]
        nRes = str("M= " + k + " pop= " + v + "\n" for k,v in res)
        return nRes

class NormalStars(Stars):
    """
    Star under the 8 solar masses threshold
    """
    def __init__(self):
        super().__init__()
        self.defLife = None
        pass

    def endStars(self, starGroup):
        sYield = getStarYield(starGroup[3].mMF(), starGroup[0])
        renmant = sYield.pop("rem")
        starGroup[3].gEx(renmant*cUnit*starGroup[1])
        for k, v in sYield.items():
            if v < 0:
                try:
                    starGroup[3].ex(k, -v)
                    continue
                except ArithmeticError:
                    starGroup[3].sMS(k, 0)
                    continue
            starGroup[3].con("h", k, v*cUnit*starGroup[1])
        product = starGroup.pop(3)#starGroup[3].gEx(starGroup[1]*starGroup[0])
        # self.ade.append(self.SG.pop(self.SG.index(starGroup)))
        self.SG.pop(self.SG.index(starGroup))
        return product, (renmant, starGroup[1])


class HugeStars(Stars):
    """
    Star over the 8 solar masses threshold
    """
    def __init__(self):
        super().__init__()
        self.defLife = 1
        pass
    
    def endStars(self, starGroup):
        """generate the supernovae"""
        sYield = getStarYield(starGroup[3].mMF(), starGroup[0])
        renmant = sYield.pop("rem")
        starGroup[3].gEx(renmant*cUnit)
        for k, v in sYield.items():
            if v < 0:
                try:
                    starGroup[3].ex(k, -v)
                    continue
                except ArithmeticError:
                    starGroup[3].sMS(k, 0)
                    continue
            starGroup[3].con("h", k, v*cUnit*starGroup[1])
        product = starGroup.pop(3)#starGroup[3].gEx(starGroup[1]*starGroup[0])
        # self.ade.append(self.SG.pop(self.SG.index(starGroup)))
        self.SG.pop(self.SG.index(starGroup))
        return product, (renmant, starGroup[1])


class NeutronStars(Stars):
    def __init__(self):
        super().__init__()
        self.defLife = None
        pass
    def spinVeryFast(self):
        pass
    def update(self):
        """
        Update every star group in the class, reducing by a Time Tick the life of each of them,
        return the product of the died star as a composition
        """
        rn = []
        upP = Composition()
        return upP, rn

class WhiteDwarf(Stars):
    def __init__(self):
        super().__init__()
        self.defLife = None
        pass
    def update(self):
        """
        Update every star group in the class, reducing by a Time Tick the life of each of them,
        return the product of the died star as a composition
        """
        rn = []
        upP = Composition()
        return upP, rn

class BinarySystemSC():
    """
    Class for bynary system that are Snee 1a candidate
    Initial: a list of Star Groups to initialize the class with 
    """
    def __init__(self, initial=None):
        if initial == None:
            self.SG=[]
        else:
            self.SG = initial
            
        self.defLife = 10
        pass
    
    def fSis(self, mass1, mass2, lifespan=None):
        """
        Find Sister, return a int defining the index where a similar stargroup is in self.Stars
        """
        lifespan = self.defLife if lifespan == None else lifespan
        for s in self.SG:
            if s[0] == mass1 and s[1] == mass2 and s[3] == lifespan:
                return self.SG.index(s)
        raise KeyError
        
    def addStars(self, Comp1, Comp2, mass1, mass2, number, lifespan1=None, lifespan2=None):
        """
        addStars: it will try to find a sister with the same int mass, the same
        int lifespam,and if there is one, it will increase the number in that group
        Otherwise, create a new group.
        
        By default, addStars should receive mass from the update call by an equal ammount of mass removed from the
        system
        """
        lifespan1, lifespan2 = (self.defLife, self.defLife) if lifespan1 == None or lifespan2 == None else lifespan1, lifespan2
        # try:
        #     self.SG[self.fSis(mass, lifespan)][1] += number
        #     self.SG[self.fSis(mass, lifespan)][3] += Comp
        # except KeyError:
        #     self.SG.append([mass, number, lifespan, Comp])
        # pass
        self.SG.append([mass1, mass2, number, lifespan1, lifespan2, Comp1, Comp2])
        pass

    def getTotComp(self):
        """
        Return the total composition of the the same star class
        """
        res = Composition()
        for s in self.SG:
            res += s[5]
            res += s[6]
        return res
    
    def getMass(self):
        """
        Get the total mass of the star class, in solar masses
        """
        mass = 0
        for s in self.SG:
            mass += (s[0]+s[1])*s[2]
        return mass
    
    def getPop(self):
        """
        Return the total number of star of the same star class
        """
        number = 0
        for s in self.SG:
            number += s[2]
        return number

    def endStars(self, starGroup):
        # self.ade.append(self.SG.pop(self.SG.index(starGroup)))
        pass

    def update(self):
        """
        Update every star group in the class, reducing by a Time Tick the life of each of them,
        return a tuple of the product of the died star as a composition and the mass of the renmant
        """
        rn = []
        upP = Composition()
        if len(self.SG) == 0:
            return upP, rn
        for s in self.SG.copy():
            s[2] -= TIMETICK
#            if s[2] <= 0:
#                p, r = self.endStars(s)
#                rn.append(r)
#                upP += p
        return upP, rn



class SimpleModelUnit(object):
    """
    Single Simple model:
    composition: a Composition class defining the initial "free" mass of the sistem
    name: a string defining the name (used to keep track of more than 1 model during massive computation)
    surface: an int defining the parsec squared surfaces
    SFR: stellar formation rate, defined in solar masses per mil. years per squared parsec
    age: the initial timeTICK age of the simulation
    edges: link to the other simple model
    stars: a list containig the initial star status
    IMF: a string defining the chosen IMF definition (chose between: "Salpeter", "Scalo", "Kroupa")  
    SMLmethod, define the solar masses range to integrate for the IMF
    sT, define the calculation for the star lifespam (chose between: "Default", "Dimple", "SimpleClassic")
    """
    def __init__(self, composition, name="SimpleModel", surface=1, SFR=1, bP = 0.85, age=0,
                 edges = None, stars = None, IMF = "Salpeter", SMLmethod = "default", sT = "default"):
        
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
        
        self.stars = [NormalStars(), HugeStars(), BinarySystemSC(), NeutronStars(), WhiteDwarf()]
        if stars != None:
            self.stars = stars
        self.totalMassStarBorn = 0
        if stars != None:
            for s in stars:
                self.totalMassStarBorn += s.getMass()
        
        self.SFR = SFR
        """gsm, generated star masses, normalize the SFR in the current situation.
        moltiplied by the current surface density of gas of the simple model"""
        self.gsm = round((SFR*(surface*(yearPERTICK/1000000000)))*self.getSMDG(), 3)
        
        """The percentage of stars born in a binary system"""
        self.bP = bP
        
        """bS, borned star (mass), collect the result of the stellar formation and keep storred
        the values not reaching the unit"""
        self.bS = [0]*len(self.getSML(self.SMLmethod))
        """bBS, borned Binary Star, collect the result of the stellar formation for the binary system"""
        self.bBS = [0]*len(self.getSML(self.SMLmethod))
        pass
    
    def totalMass(self):
        """
        Return the total mass of the model in solar masses
        """
        mass = self.composition.totalMass()/cUnit
        for s in self.stars:
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
        return self.composition.totalMass()/self.totalMass()
    
    def globMetal(self):
        """
        global metallicity
        """
        met = 0
        met += self.composition.metalMass()
        for s in self.stars:
            met += s.composition.metalMass
        return met/self.gasMass()
    
    def getSurface(self):
        """
        Get the surface (pc)**2
        """
        return self.surface
    
    def getTMSB(self):
        """
        return the current total mass of all born star
        """
        return self.totalMassStarBorn
    
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
    
    def getSMDG(self):
        """
        surface mass density of gas (solar masses)
        """
        return self.gasMass()/cUnit/self.surface
    
    def getSML(self, method="default"):
        """
        return the solar masses list in consideration, the format is: ([1]:
        """
        if method == "default":
            return [(0.08, 0.1, 0.12), (0.12, 0.14, 0.16), (0.16, 0.18, 0.20),
                    (0.20, 0.24, 0.28), (0.28, 0.32, 0.36), (0.36, 0.40, 0.44),
                    (0.44, 0.52, 0.60), (0.60, 0.68, 0.76), (0.76, 0.84, 0.92),
                    (0.92, 1.08, 1.24), (1.24, 1.40, 1.56, ), (1.56, 1.72, 1.88),
                    (1.88, 1,94 , 2.00), (2.00, 2.25 , 2.50), (2.50, 2.75 , 3.00),
                    (3.00, 3.25 , 3.50), (3.50, 3.75 , 4.00), (4.00, 4.50, 5.00),
                    (5.00, 5.50, 6.00), (6.00, 6.50, 7.00), (6.00, 7.00, 8.00),
                    (8.00, 9.00, 10.00), (10.00, 10.00, 12.00), (12.00, 14.00, 16.00),
                    (16.00, 18.00, 20.00), (20.00, 22.00, 24.00), (24.00, 28.00, 32.00),
                    (32.00, 36.00, 40.00), (40.00, 50.00, 60.00), (60.00, 70.00, 80.00),
                    (80.00, 90.00, 100.00), (100.00, 110.00, 120.00)] #debug test
        pass
    
    def getIM(self):
        """
        get initial masses:
        it will return a list with the number of star in each mass category
        """
        ms = self.getSML()
        pr = []
        for m in ms:
            # pr.append(1/(m**(1+self.IMFdef(m))))
            pr.append(float(integrate.quadrature(lambda x : 1/(x**(1+self.IMFdef(x))),m[0], m[2])[0]))
        a = sum(pr)
        for p, b in zip(pr, self.bS):
           self.bS[self.bS.index(b)] += (p/a)*self.gsm
        return self.bS
    
    def getST(self, m):
        """
        get the stellar lifespam (stellar time), based on the mass
        """
        # if m <= 6.6:
        #     a = (10**((1.338-(1.790-(0.2232*(7.764-np.log10(m))))**
        # (1/2))/0.116)-9)*(1000000000/yearPERTICK)
        #     return int(a // yearPERTICK)
        if self.sT == "Default":
            if m <= 10:
                a = 1.2*(10**10)*(m**-2.78)
                return int(a // yearPERTICK)
            else:
                a = 1.1*(10**8)*(m**-0.75)
                return int(a // yearPERTICK)
        elif self.sT == "Simple":
            a = 11.7*((10^9)/(m^2))
            return int(a // yearPERTICK)
        elif self.sT == "SimpleClassic":
            return 1
        
    def stellarGen(self, s, m, index):
        c = int(s//m[1])
        si = self.bS.index(s)
        self.bS[si] = s%m[1]
        self.bBS[si] = int(c*self.bP)
        c = c - int(c*self.bP)
        try:
            a = self.composition.gEx(round(c*m[1]*cUnit, cDP))
        except AssertionError:
            c *= self.composition.totalMass()/(c*m[1]*cUnit)
            c = int(c//1)
            a = self.composition.gEx(round(c*m[1]*cUnit, cDP))
        self.stars[index].addStars(a, m[1], c, self.getST(m[1]))
        self.totalMassStarBorn += (m[1]*c)//m[1]
        self.gsm = round((self.SFR*(self.surface*(yearPERTICK/1000000000)))*self.getSMDG(), 3)
        pass

    def bStellarGen(self, bBS, SML):
        c = 0
        while sum(bBS) > 1:
            a = choice(range(len(bBS)))
            if len(bBS) == 1:
                a = choice(range(len(bBS)))
                pass
            if bBS[-1] > 0:
                x = SML[-1][1]
                y = SML[a][1]
                c = 1 + int(bBS[-1]*random()) if 1 + int(bBS[-1]*random()) < bBS[-1] else bBS[-1]
                if bBS[a] < c and bBS[a] != bBS[-1]:
                    c = bBS[a]
                elif bBS[a] == bBS[-1]:
                    c = int(bBS[a]//2)
                bBS[a] -= c
                bBS[-1] -= c
                xComp = self.composition.gEx(round(c*x*cUnit, cDP))
                yComp = self.composition.gEx(round(c*y*cUnit, cDP))
                self.stars[2].addStars(xComp, yComp, x, y, c, self.getST(x), self.getST(y))
            else:
                self.bStellarGen(bBS[:-1], SML[:-1])
                break
        pass

    def update(self):
        """
        Update the current simple unit composition
        firstly, elaborate the yield of the current star population, both
        composition and stellar renmant creation
        secondly, calculate the stars born in this unitTime
        do not return anithing
        """
        self.age += 1
        
        for s in self.stars:
            a, r = s.update()
            self.composition += a
            self.gsm = round((self.SFR*(self.surface*(yearPERTICK/1000000000)))*self.getSMDG(), 3)
            if s.__class__.__name__ == "NormalStars":
                for rn in r:
                    c = Composition()
                    self.stars[4].addStars(c.aMS("h", rn[0]*rn[1]*cUnit), rn[0], rn[1])
            if s.__class__.__name__ == "HugeStars":
                for rn in r:
                    c = Composition()
                    self.stars[3].addStars(c.aMS("h", rn[0]*rn[1]*cUnit), rn[0], rn[1])
        for s, m in zip(self.getIM(), self.getSML(self.SMLmethod)):
            if m[1] <= 8 and s/m[1] >= 1:
                self.stellarGen(s,m,0)
            elif s >= 1:
                self.stellarGen(s,m,1)
        self.bStellarGen(self.bBS, self.getSML(self.SMLmethod))
        pass
    
    def Run(self, times):
        for t in range(times):
            self.update()
            # if t % 10 == 0:
            #     print(str(t) + "\n"+ str(self.composition))
            #     print(str(self.totalMass()) + " "
            #           + str(self.gsm) + " " + str(self.composition.mMF()))
            #     print()
            if t % 100 == 0:
                print("pop at " + str(t) + "M years:")
                print(self.stars[0])
                print(self.stars[1])
                print("Total Pop = " + str(self.stars[0].getPop() + self.stars[0].getPop()) + "\n")
        pass

a = Composition()
a.aMS("H", 100000000000*cUnit)
fi = SimpleModelUnit(a, surface = 75398, SFR=0.1)
fi.Run(13000)
