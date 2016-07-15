# generic objects and utilities for dealing with pdb files
import re
from numpy import *

# regular expression for an atom in a pdb file
atomR = re.compile('^(ATOM|HETATM) *([0-9]+) *(\S+) +([A-Z]+) +([a-zA-Z]?) *\
([-0-9]+) *(-?[0-9.]+) *(-?[0-9.]+) *(-?[0-9.]+)(.*$)')

# regular expression for a biological symmetry transform line
transformR = re.compile('^REMARK +350 +BIOMT[1-3] *(\d+) +(-?[0-9.]+) +(-?[0-9.]+) +(-?[0-9.]+) +(-?[0-9.]+) *$')


class Atom:
    """This object represents information concerning a single atom"""
    def __init__(self,pdbmatch = None,coords = None, name = 'X',num=0):
        """Fill in information using a match object from a pdb file"""
        if pdbmatch != None:
            # base the atom information on a match object from a pdb line
            self.type = pdbmatch.group(1)
            self.num = int(pdbmatch.group(2))
            self.name = pdbmatch.group(3)
            self.res = pdbmatch.group(4).strip()
            self.chain = pdbmatch.group(5)
            self.resnum = int(pdbmatch.group(6))
            self.x = float(pdbmatch.group(7))
            self.y = float(pdbmatch.group(8))
            self.z = float(pdbmatch.group(9))
            self.coords = array([self.x,self.y,self.z])
            self.occupancy = None
            self.bfactor = None
            self.tail = pdbmatch.group(10)

        elif coords != None:
            # give the atom the desired coordinates and type
            self.num = num
            self.coords = coords
            [self.x,self.y,self.z] = coords
            self.name = name
            self.res = 'SSN'
            self.resnum = 0            
            self.chain = 'X'
            self.occupancy = 1
            self.bfactor = 1
            self.tail = '            C'
            self.type = 'HETATM'

        self.coords = array(self.coords)
        self.chainobj = None
        self.resobj = None
        self.conect = []
        
    def __repr__(self):
        return "<Atom %d, %s; res %d; chain %s>" %(self.num, self.name, self.resnum, self.chain)

    def pdbline(self):
        """Get the pdb line corresponding to this atom"""
        
        if (self.occupancy==None or self.bfactor == None):            
            line = '%s%5d  %s%3s%2s%4d%12.3f%8.3f%8.3f%s\n' \
               %(self.type.ljust(6),self.num,self.name.ljust(4),self.res,self.chain,self.resnum,\
                 self.coords[0],self.coords[1],self.coords[2],self.tail)
        else:
            line = '%s%5d  %s%3s%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%s\n' \
               %(self.type.ljust(6),self.num,self.name.ljust(4),self.res,self.chain,self.resnum,\
                 self.coords[0],self.coords[1],self.coords[2],self.occupancy,self.bfactor,self.tail)

        return line

    def conline(self):
        """Get a CONECT line for this atom"""
     
        str = 'CONECT%5d' %(self.num)
        for a in self.conect:
            str += '%5d' %a.num
        str += '\n'
        return str
    
class Residue:
    """Residue or nucleotide object"""
    def __init__ (self,atoms=None):
        if atoms == None:
            self.atoms = []
            self.num = None
            self.chain = None
            self.name = None
        else:
            self.atoms = atoms
            self.chain = atoms[0].chain
            self.num = atoms[0].resnum
            self.name = atoms[0].res

    def __repr__(self):
        return "<Residue %d: %s>" %(self.num, self.name)

    def getAtomByName(self,name,regexp=0):
        """Return the first atom with the given name found in the residue. If regexp is true, then get the atom whose name matches the regular expression"""
        if regexp:
            rN = re.compile(name)
            for a in self.atoms:
                if rN.search(a.name) != None:
                    return a
        else:
            for a in self.atoms:
                if a.name == name:
                    return a
            
        return None
        
class Chain:
    """Chain object"""
    def __init__(self,residues = None):
        if residues == None:
            self.residues=[]
            self.name = None
        else:
            self.residues = residues
            self.name = self.residues[0].chain
    def __repr__(self):
        return "<Chain %s>" %self.name

    def fromAtomList(self,atmlist):
        # set up a chain from a list of atom objects
        # split them up into residues
        # raises error if they don't all have the same chain name

        self.residues = []
        self.name = atmlist[0].chain

        start = 0
        for a in atmlist:           
            if a.chain != self.name:
                raise ValueError("Atmlist does not all have the same chain name %s %s" %(self.name, a.chain))

            if not start or a.resnum!= res.num:
                if start:
                    self.residues.append(res)
                res = Residue()
                res.name = a.res; res.num = a.resnum
                res.chain = a.chain;
                res.atoms = [a]
                start = 1
            else:
                res.atoms.append(a)

        self.residues.append(res)

        return 0
class BasePair:
    """DNA base pair object. Contains 2 residues and a list of atoms"""
    def __init__(self,residues = None):
        self.atoms = []
        if residues == None:
            self.residues = []
        else:
            self.residues = residues
            for r in residues:
                self.atoms.extend(r.atoms)

        if len(self.residues) != 0 and len(self.residues) != 2:
            print "WARNING: creating base-pair with neither 0 nor 2 residues"

    def __repr__(self):
            return "<BasePair: %s %s>" %(self.residues[0].name,self.residues[1].name)

def goodPair(r1,r2):
    """Check if two DNA/RNA residues make a correct basepair"""

    return (r1.name[-1]  == 'A' and r2.name[-1] == 'T') \
        or (r1.name[-1]  == 'T' and r2.name[-1] == 'A') \
        or  (r1.name[-1]  == 'G' and r2.name[-1] == 'C') \
        or  (r1.name[-1]  == 'C' and r2.name[-1] == 'G') \
        or (r1.name=='THY' and r2.name=='ADN') \
        or (r1.name=='ADN' and r2.name=='THY') \
        or (r1.name=='CYT' and r2.name=='GUA') \
        or (r1.name=='GUA' and r2.name=='CYT')

class Structure:
    """Class containing an molecule or multi-molecule structure"""
    def __init__(self,infile = None):
        self.initvars()

        if infile != None:
            self.structFromFile(infile)

    def renumRes(self):
        """Renumber the residues"""
        for c in range(len(self.residues)):
            self.residues[c].num = c+1
            for a in self.residues[c].atoms:
                a.resnum = c+1

    def resetFromResidues(self):
        # reset the chains and atoms in a structure when the residue list is set

        # rebuild the list of chains
        # dictionary mapping chain objects to their names
        self.chains = []; 
        chainnames = {}

        for r in self.residues:
            if r.chain in chainnames.keys():
                chainnames[r.chain].residues.append(r)
                chainnames[r.chain].atoms.extend(r.atoms)
            else:
                chainnames[r.chain] = Chain()
                self.chains.append(chainnames[r.chain])
                chainnames[r.chain].name = r.chain
                chainnames[r.chain].residues = [r]
                chainnames[r.chain].atoms = r.atoms[:]
            for a in r.atoms:
                a.chain = r.chain

        # rebuild atom list 
        self.atoms = []
        [self.atoms.extend(r.atoms) for r in self.residues]

    def resetFromChains(self):
        # reset the residues and atoms based on the chains

        self.residues = []; self.atoms = []
        for c in self.chains:
            self.residues.extend(c.residues)
            for r in c.residues:
                r.chain = c.name
                for a in r.atoms:
                    a.chain = c.name
            [self.atoms.extend(r.atoms) for r in c.residues]

        return 0

    def PCA(self):
        """Perform a principal component analysis on the atom coordinates. Returns eigenvals and eigenvecs. Sorted from largest to smallest eigenvalue"""
        
        M = array([a.coords for a in self.atoms])
        # covariance matrix
        covmat = cov(M,rowvar=0)
        
        # find the eigenvalues and eigenvecs
        (eval,evec) = linalg.eig(covmat)

        # sort eval and evec
        ind = argsort(-eval)
        eval = eval[ind]
        evec = evec[:,ind]

        # make sure it forms a right-handed coordinate system
        check = dot(cross(evec[:,0],evec[:,1]),evec[:,2])
        
        if (check < 0):
            ind = [1,0,2]
            evec = evec[:,ind]

        return (eval,evec)

    def rotateM(self,M):
        """Rotate entire structure by the given rotation matrix"""

        rotmat = matrix(M).T
        for a in self.atoms:            
            a.coords = array(a.coords*rotmat)[0]
            
    def atomByBum(self,num):
        # get the atom of the given number
        for a in self.atoms:
            if a.num == num:
                return a
        return None
        
    def chainByName(self,name):
        # get the first chain with the given name

        for c in self.chains:
            if c.name==name:
                return c

    def initvars(self):
        """Initialize various variables"""

        # extra lines at start/end of structure
        self.startlines = []
        self.endlines = []
        self.chains = []
        self.residues = []
        self.atoms = []
       
        # transformations
        self.transforms = []

    def outputPDB(self,outfile,append=0,ident=1):        
        # output a pdb file for this structure
        # if append then append to the pdb file
        # ident is the model identifier to use

        if (append):
            OF = open(outfile,'a')
        else:
            OF = open(outfile,'w')

        OF.write('MODEL     %4d\n' %ident)
        for l in self.startlines:
            OF.write(l)

        if len(self.chains) == 0:
            [OF.write(a.pdbline()) for a in self.atoms]
        else:
            # write in all the chain atoms
            for c in self.chains:
                for r in c.residues:
                    for a in r.atoms:
                        OF.write(a.pdbline())
                        #if hasattr(c,'ter'):
                        #    OF.write(c.ter)
                            
            # write in all the extra atoms (not part of chains)
            [OF.write(a.pdbline()) for a in self.atoms if a.chainobj == None]

        # write out any conect records
        [OF.write(a.conline()) for a in self.atoms if len(a.conect) > 0]     
                
        for l in self.endlines:
            OF.write(l)
        
        OF.write('ENDMDL\n')
        OF.close()
    
    def structFromFile(self,infile):
        # given an input file, get structure based on the chains and residues found in that file
        # WARNING: for now preserves only direct connectivity information
        # no hydrogen bonds or salt bridges
            
        IF = open(infile)
        lines = IF.readlines()
        IF.close()

        self.transforms = []

        started = 0
        for line in lines:
            # check for symmetry transform lines
            m = transformR.match(line)
            if m != None:
                tn = int(m.group(1))
                matrow = [float(m.group(i)) for i in range(2,5)]
                shift = float(m.group(5))
                if tn > len(self.transforms):
                    self.transforms.append([[matrow],[shift]])
                else:
                    self.transforms[tn-1][0].append(matrow)
                    self.transforms[tn-1][1].append(shift)
                self.startlines.append(line)
                continue

            m = atomR.match(line)
            #m2 = conectR.match(line)
            if m == None:
                if not started:
                    self.startlines.append(line)
                elif line[:3] == 'TER':
                    newchain.ter = line
                elif line[:6] =='CONECT':
                    # get the atom to which this record pertains
                    atm = self.atomByNum(int(line[6:11]))
                    if atm == None:
                        print "Cannot find atom number : ", int(line[6:11])
                    else:
                        # get the atoms connected to it (max of 4)
                        for a in range(4):
                            str = line[11+5*a:16+5*a]
                            try:
                                atm.conect.append(self.atomByNum(int(str)))
                            except ValueError:
                                break                    
                else:
                    self.endlines.append(line)

                continue
            newatom = Atom(m)
            if not started:
                newres = Residue([newatom])
                newchain = Chain([newres])
                self.residues.append(newres)
                self.chains.append(newchain)
                started = 1
            else:
                if newatom.resnum != self.atoms[-1].resnum:
                    newres = Residue([newatom])
                    self.residues.append(newres)
                    if newatom.chain != self.atoms[-1].chain:
                        newchain = Chain([newres])
                        if newatom.chain.strip() != '':
                            self.chains.append(newchain)
                    else:
                        newchain.residues.append(newres)
                else:
                    newres.atoms.append(newatom)

            self.atoms.append(newatom)
            if newatom.chain.strip() != '':
                newatom.chainobj = newchain
            newatom.resobj = newres
            
            self.endlines = []
                
            for t in self.transforms:
                t[0] = array(t[0]); t[1] = array(t[1])

    def renumAtoms(self,start=1):
        """Renumber all atoms from 1 or another starting point"""
        for i in range(len(self.atoms)):
            self.atoms[i].num = start + i

    def atomByNum(self,num):
        """Get the atom with the corresponding number"""
        for a in self.atoms:
            if a.num == num:
                return a
        return None
        
    def getCOM(self,countHet=0):
        """Get the center of mass of all the atoms in the structure. If countHet is true, include heteroatoms as well"""
        
        if countHet:
            alist = self.atoms
        else:
            alist = [a for a in self.atoms if a.type == 'ATOM']
            
        com = [0,0,0]
        for i in range(3):
            com[i] = sum([a.coords[i] for a in alist])/float(len(alist))
            
        return com

    def recenter(self):
        """Move entire structure to shift COM to 0"""
        com = self.getCOM(countHet=1) # center of mass of whole thing
        # shift everything by center of mass
        for a in self.atoms:
            a.coords = a.coords - com

def makeBareStruct():
    """Make a barebones structure with nothing in it, but containing appropriate start/end stuff so that it can be output to pdb"""

    struct = Structure()
    #struct.startlines = ['HET    SSN  A   1    1402     Pseudo atom representation of DNA\n', \
 #                            'HETNAM     SSN Body and ribbon spatial coordinates\n',\
  #                           'FORMUL  1   SSN    C20 N20 P21\n']
    struct.startlines = ['HEADER    structure\n']
    struct.endlines = ['END\n']

    return struct
