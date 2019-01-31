"""Simulate stocatic positioning of loop extrussion factors.

This module implements a Gillespie algorithm to simulate positioning of
loop extrusion factors on a chromatin chain.

Example:
    beadsPerLoop = 500
    processivity = 500.0
    nloops = int(np.floor(L/beadsPerLoop))
    chain = loopSim.LoopExtrusionChain(L, nloops, CTCF_file)

    forward_rate = 1.0
    reverse_rate = 0.0
    falloff_rate = (forward_rate-reverse_rate)/processivity
    delta_t = 20.0/falloff_rate
    chain.run_for_time(delta_t, forward_rate, reverse_rate,
                           falloff_rate)
"""
import random
import numpy as np

class CohesinMolecule:
    """Structure to keep track of a single Cohesin molecule

    Attributes:
        position (:obj:`list` of int): [left, right] position of feet
        for_ind (:obj:`list` of int): [left, right] index in forward_moves
        rev_ind (:obj:`list` of int): [left, right] index in reverse_moves
    """
    def __init__(self,left,right):
        self.position = [left,right]
        self.for_ind = [None,None]
        self.rev_ind = [None,None]

    def __str__(self):
        return str(self.for_ind[0]) + "<-" + str(self.position[0]) + "->" \
                + str(self.rev_ind[0]) + "\t || \t" +  str(self.rev_ind[1]) \
                + "<-" + str(self.position[1]) + "->" + str(self.for_ind[1])

    def is_forward(self,location,foot):
        """Returns true if foot steping to location is a forward move
        """
        forward = (location-self.position[foot])*(foot*2-1)
        if forward == 1:
            return True
        else:
            return False

    def next_position(self,foot,reverse=False):
        if reverse:
            return self.position[foot] - [-1,1][foot]
        else:
            return self.position[foot] + [-1,1][foot]



class LoopExtrusionChain:
    """Data Structure for Gillespie algorithm for cohesins walking along a
    strand of DNA.


    Attributes:
        cohesins (:obj:`list` of :obj:`CohesinMolecule`): List of
            CohesinMolecule objects.
        fowrad_moves (:obj:`list` of :obj:`tuple` of :obj:`int`): Each element
            is of the form (cohesin index, foot)
        reverse_moves (:obj:`list` of :obj:`tuple` of :obj:`int`): Each element
            is of the form (cohesin index, foot)
        beads (:obj:`list`): None if no foot is on that bead
        blocks (:obj:`list`): -1 for block left feet, +1 for block right feet, 2 for both
    """
    def __init__(self, nbeads, ncohesin=None, CTCF_file=None):
        """ Create a chain with nbeads beads

        Args:
            nbeads (int): Number of beads in chain
            ncohesin (int, optional): Number of cohesin to insert
            CTCF_file (str, optional): File name of CTCF locations, see
                load_CTCF_file
        """
        self.beads = [None] * nbeads
        self.blocks = [None] * nbeads
        if type(CTCF_file) != type(None):
            self.load_CTCF_file(CTCF_file)
        self.forward_moves = []
        self.reverse_moves = []
        self.cohesins = []
        if (ncohesin==None):
            return
        for ii in range(ncohesin):
            self.add_cohesin()

    def load_CTCF_file(self,file_name):
        """Load CTCF locations from text file.

        """
        nbeads = len(self.beads)
        if not hasattr(self, 'blocks'):
            self.blocks = [None] * nbeads

        with open(file_name) as f:
            for line in f:
                cols = line.split()
                try:
                    index = int(cols[0])
                except ValueError:
                    raise ValueError(
                        "Expected int in first column of CTCF file, got "
                         + cols[0])
                direction = cols[1]
                if direction in [">", "0", "L"]:
                    # Forward CTCF motif, prevents leftward motion
                    # 5'- CCACNAGGTGGCAG - 3'
                    blocks[index] = 0
                elif direction in ["<", "1", "R"]:
                    # Reverse CTCF motif, prevents rightward motion
                    # 5' - CTGCCACCTNGTGG - 3'
                    blocks[index] = 1
                elif direction in ["x", "2", "B"]:
                    # Both CTCF motif dicrections on same bead
                    blocks[index] = 2
                else:
                    raise ValueError(direction+" is not a valid direction")

    def check_consistancy(self):
        """Does internal check to see if structure is self consistant.
        """
        ii=0
        for for_move in self.forward_moves:
            if self.cohesins[for_move[0]].for_ind[for_move[1]] != ii:
                print("for_move",for_move)
                print("ii",ii)
                print("should be ",self.cohesins[for_move[0]].for_ind[for_move[1]])
                raise ValueError("Internal Error")
            ii=ii+1
        ii=0
        for rev_move in self.reverse_moves:
            if self.cohesins[rev_move[0]].rev_ind[rev_move[1]] != ii:
                print("rev_move",rev_move)
                print("ii",ii)
                print("should be ",self.cohesins[rev_move[0]].rev_ind[rev_move[1]])
                raise ValueError("Internal Error")
            ii=ii+1
        ii=0
        for my_cohesin in self.cohesins:
            for foot in [0,1]:
                if my_cohesin.for_ind[foot] == None:
                    next_bead = my_cohesin.next_position(foot)
                    if self.is_available(next_bead,foot):
                        raise ValueError("Move shoud be OK")
                else:
                    if self.forward_moves[my_cohesin.for_ind[foot]] != (ii,foot):
                        print(self)
                        self.display()
                        print("cohesin",ii,"=",str(my_cohesin))
                        print("LHS ",self.forward_moves[my_cohesin.for_ind[foot]])
                        print("RHS ", (ii,foot))
                        raise ValueError("forward moves inconsistant with cohesins")
                if my_cohesin.rev_ind[foot] == None:
                    next_bead = my_cohesin.next_position(foot,reverse=True)
                    if self.is_available(next_bead,foot):
                        print(self)
                        self.display()
                        print("ii",ii,"foot",foot)
                        print("next",my_cohesin.next_position(foot,reverse=True))
                        print("LHS",self.beads[my_cohesin.next_position(foot,reverse=True)])
                        raise ValueError("Move shoud be OK")
                else:
                    if self.reverse_moves[my_cohesin.rev_ind[foot]] != (ii,foot):
                        raise ValueError("reverse moves inconsistant with cohesins")
            ii=ii+1
        ii=0
        for bead in self.beads:
            if bead != None:
                if self.cohesins[bead[0]].position[bead[1]] != ii:
                    raise ValueError("beads not consistant with cohesins")
            ii=ii+1




    def __str__(self):
        """Produces a string represntation of the cohesin loops

        e.g.  () (4)  ((4)56) 16(3)
        Each ( and ) correispond to left and right feet.
        A space prepresents a empty bead
        A number represents that number of empty beads
        """
        indexes = []
        cohesin = []
        for ii in range(0,len(self.beads)):
            if self.beads[ii] != None:
                indexes.append(ii)
                if self.beads[ii][1]==0:
                    cohesin.append("(")
                elif self.beads[ii][1]==1:
                    cohesin.append(")")
                else:
                    cohesin.append("?")
            if self.blocks[ii]==0:
                indexes.append(ii)
                cohesin.append(">")
            elif self.blocks[ii]==1:
                indexes.append(ii)
                cohesin.append("<")
            elif self.blocks[ii]==2:
                indexes.append(ii)
                cohesin.append("x")
        iostr = ""

        for ii in range(len(cohesin)-1):
            iostr = iostr + cohesin[ii]
            if indexes[ii+1]==indexes[ii]+1:
                continue
            if indexes[ii+1]==indexes[ii]+2:
                iostr=iostr+" "
                continue
            if indexes[ii+1]==indexes[ii]+3:
                iostr=iostr+"  "
                continue
            iostr=iostr+str(indexes[ii+1]-indexes[ii]-1)
        iostr = iostr + cohesin[-1]
        return iostr

    def display(self):
        """Prints a multi-line discription of the cohesins on the chain.
        """
        print("cohesins")
        ii =0
        for me in self.cohesins:
            print(me.for_ind[0], "<-", me.position[0], "->", me.rev_ind[0],
                  "\t |"+str(ii)+"| \t", me.rev_ind[1],"<-" ,me.position[1],
                  "->", me.for_ind[1])
            ii=ii+1
        print("forward moves")
        for me in self.forward_moves:
            print("Cohesin Index",me[0]," "+['left','right'][me[1]])
        print("reverse moves")
        for me in self.reverse_moves:
            print("Cohesin Index",me[0]," "+['left','right'][me[1]])

    def total_spand_beads(self):
        total = 0
        for my_cohesin in self.cohesins:
            total = total + my_cohesin.position[1] - my_cohesin.position[0] + 1
        return total

    def print_for_Monte_Carlo(self,filename='bindpairs',offset=1):
        my_file = open(filename,mode='w')
        for ii in range(len(self.beads)):
            if self.beads[ii] == None:
                print("-1",file=my_file)
            else:
                (coh_ind,foot) = self.beads[ii]
                print(str(self.cohesins[coh_ind].position[1-foot]+offset),
                      file=my_file)
        my_file.close()

    def is_available(self,bead,LR=None):
        """True if the bead is available for a cohesin to step onto.

        Args:
            bead (:obj:`int`): The index of the bead to be tested.
            LR (bool, optional): 0 for Left cohesin foot, 1 for right cohesin
                foot.
        Returns:
            bool: True if bead is available for cohesin to step there.
        """
        if (bead<0):
            return False
        if (bead>=len(self.beads)):
            return False
        if (self.beads[bead] != None):
            return False
        if (LR == self.blocks[bead] or self.blocks[bead]==2):
            return False
        return True

    def cohesin_located_at(self,bead):
        """Returns the (index, foot) of cohesin located at bead
        """
        if (bead<0):
            return None
        if (bead>=len(self.beads)):
            return None
        if (self.beads[bead] != None):
            return self.beads[bead]
        return None


    def delete_forward_move(self,to_delete):
        """Remove to_delete from forward_moves

        This is accomplied efficently by moving the end element of the list to
        the  position of the element being deleted.
        Also updates the forward move index in the cohesin list for the move
        which was previously at the end of the list.

        Args:
            to_delete (int): Index of forward move to be deleted
        """
        # copy move from end of list to position being deleted
        self.forward_moves[to_delete] = self.forward_moves[-1]
        # The cohesin which used to have its forward move at the end of
        # the forward_moves list now must point to to_delete
        (coh_ind,foot) = self.forward_moves[-1]
        self.cohesins[coh_ind].for_ind[foot] = to_delete
        # remove the duplicate from the end of the list
        del self.forward_moves[-1]

    def delete_reverse_move(self,to_delete):
        """Remove to_delete from reverse_moves

        This is accomplied efficently by moving the end element of the list to
        the  position of the element being deleted.
        Also updates the reverse move index in the cohesin list for the move
        which was previously at the end of the list.

        Args:
            to_delete (int): Index of reverse move to be deleted
        """
        # copy move from end of list to position being deleted
        self.reverse_moves[to_delete] = self.reverse_moves[-1]
        # The cohesin which used to have its reverse move at the end of
        # the reverse_moves list now must point to to_move
        (coh_ind,foot) = self.reverse_moves[-1]
        self.cohesins[coh_ind].rev_ind[foot] = to_delete
        # remove the duplicate from the end of the list
        del self.reverse_moves[-1]



    def prevent_motion(self,from_bead,to_bead):
        """Prevents feet from moving from from_bead to to_bead.
        Edits: self.cohesins, self.forward_moves, self.reverse_moves
        """
        to_stop = self.cohesin_located_at(from_bead)
        if (to_stop != None):
            # The cohesin you ran into can move
            # which cohesin did you run into
            (coh_ind_other,foot_other)=to_stop
            # Prevent a forward or backward move
            if self.cohesins[coh_ind_other].is_forward(to_bead,foot_other):
                self.delete_forward_move(self.cohesins[coh_ind_other].for_ind[foot_other])
                self.cohesins[coh_ind_other].for_ind[foot_other]=None
            else:
                self.delete_reverse_move(self.cohesins[coh_ind_other].rev_ind[foot_other])
                self.cohesins[coh_ind_other].rev_ind[foot_other]=None

    def allow_motion(self,from_bead,to_bead):
        """Allos feet to move from from_bead to to_bead.
        Edits: self.cohesins, self.forward_moves, self.reverse_moves
        """
        # Motion of a cohesin may allow another to move to follow it
        freed = self.cohesin_located_at(from_bead)
        if (freed != None):
            (coh_ind_other,foot_other)=freed
            if self.cohesins[coh_ind_other].is_forward(to_bead,foot_other):
                if self.cohesins[coh_ind_other].for_ind[foot_other] != None:
                    raise ValueError("Should be None!")
                self.cohesins[coh_ind_other].for_ind[foot_other]=len(self.forward_moves)
                self.forward_moves.append(freed)
            else:
                if self.cohesins[coh_ind_other].rev_ind[foot_other] != None:
                    raise ValueError("Should be None!")
                self.cohesins[coh_ind_other].rev_ind[foot_other]=len(self.reverse_moves)
                self.reverse_moves.append(freed)


    def add_cohesin(self, with_index = None):
        """Insert a cohesin molecule to the chain randomly.

        Avoids placing cohesin on top of one another.
        Adds and removes possible moves as needed.

        Args:
            with_index (int, optional): The index in the cohesin list to insert
            the cohesin at.  If None it will be appended to the end.
        """
        assert (len(self.beads)>1), "must have non-trival chain for a bead"
        while (True):
            left = random.randint(0,len(self.beads)-2)
            right = left+1
            if not self.is_available(left,LR=0):
                continue
            if not self.is_available(right,LR=1):
                continue

            # my_cohesin[left,right]
            my_cohesin = CohesinMolecule(left,right)
            if with_index == None:
                cohesin_ind = len(self.cohesins)
            else:
                cohesin_ind = with_index
            # OK to insert
            self.beads[left] = (cohesin_ind,0)
            self.beads[right] = (cohesin_ind,1)

            # Can the input feet move forward?  (They can move backward.)
            if self.is_available(left-1,LR=0):
                my_cohesin.for_ind[0] = len(self.forward_moves)
                self.forward_moves.append((cohesin_ind,0))

            if self.is_available(right+1,LR=1):
                my_cohesin.for_ind[1] = len(self.forward_moves)
                self.forward_moves.append((cohesin_ind,1))

            #add to cohesin list
            if with_index == None:
                self.cohesins.append(my_cohesin)
            else:
                self.cohesins[with_index] = my_cohesin

            # Existing cohesin my nolonger be able to move
            self.prevent_motion(left-1,left)
            self.prevent_motion(right+1,right)

            return

    def get_total_rate(self,rate_forward,rate_reverse,rate_falloff):
        """Get total rate for Gillespie algorithm.

        Args:
            rate_forward (float): Rate at which cohesin feet step forward.
            rate_reverse (float): Rate at which cohesin feet step backward.
            rate_falloff (float): Rate at which cohesin fall off.
        Returns:
            float: Total rate
        """
        return (len(self.forward_moves)*rate_forward
               + len(self.reverse_moves)*rate_reverse
               + len(self.cohesins)*rate_falloff)


    def step(self,to_move,reverse=False):
        """ Move cohesin foot forward.

        Args:
            to_move (int): Index of move in self.forward_moves
            reverse (bool): Wheter to take a reverse step
        """
        # ----- Prepare for Move -----
        # what is the cohesin index and foot
        if reverse:
            (coh_ind,my_foot)=self.reverse_moves[to_move]
        else:
            (coh_ind,my_foot)=self.forward_moves[to_move]
        this_cohesin = self.cohesins[coh_ind] # not a copy

        # Make sure foot is a valid value
        assert(my_foot in [0,1]), "Not a valid value for foot"
        # Get current position of foot
        current_bead = this_cohesin.position[my_foot]
        # Get next position of foot (left foot moves -1, right foot moves +1)
        if reverse:
            # Make sure cohesin data is consistant with reverse move data
            assert(this_cohesin.rev_ind[my_foot] == to_move), "inconsisent"
            next_bead = current_bead + [1,-1][my_foot]
            freed_bead = current_bead - [1,-1][my_foot]
            nnbead = current_bead + [2,-2][my_foot]
        else:
            # Make sure cohesin data is consistant with forward move data
            assert(this_cohesin.for_ind[my_foot] == to_move), "inconsisent"
            next_bead = current_bead + [-1,1][my_foot]
            freed_bead = current_bead - [-1,1][my_foot]
            nnbead = current_bead + [-2,2][my_foot]
        # Make sure that the step is OK
        assert(self.is_available(next_bead,LR=my_foot)), "Not available"+\
            str(next_bead)+" "+str(my_foot)

        # ---- Now Make move ----
        # Update foot position in cohesin
        this_cohesin.position[my_foot]=next_bead
        # Update bead list
        self.beads[next_bead]= self.beads[current_bead]
        self.beads[current_bead]=None

        # update forward(reverse)_moves list and the pointers to it in the
        # cohesin
        if reverse:
            if this_cohesin.for_ind[my_foot] == None:
                # You can do a reverse move next time
                this_cohesin.for_ind[my_foot] = len(self.forward_moves)
                self.forward_moves.append((coh_ind,my_foot))
            if not self.is_available(nnbead,LR=my_foot):
                self.delete_reverse_move(to_move)
                # this cohesin nolonger can move backward
                this_cohesin.rev_ind[my_foot] = None

        else:
            if this_cohesin.rev_ind[my_foot] == None:
                # You can do a reverse move next time
                this_cohesin.rev_ind[my_foot] = len(self.reverse_moves)
                self.reverse_moves.append((coh_ind,my_foot))
            if not self.is_available(nnbead,LR=my_foot):
                self.delete_forward_move(to_move)
                # this cohesin nolonger can move forward
                this_cohesin.for_ind[my_foot] = None

        self.prevent_motion(nnbead,next_bead)
        self.allow_motion(freed_bead,current_bead)


    def falloff(self,to_move):
        """Takes a cohesin off and re-inserts it randomly
        """
        doomed_cohesin = self.cohesins[to_move]
        [left,right] = doomed_cohesin.position
        self.beads[left]=None
        self.beads[right]=None
        for foot in [0,1]:
            if doomed_cohesin.for_ind[foot] != None:
                self.delete_forward_move( doomed_cohesin.for_ind[foot] )
            if doomed_cohesin.rev_ind[foot] != None:
                self.delete_reverse_move( doomed_cohesin.rev_ind[foot] )

        #other feed my nolonger be block by removed cohesin
        self.allow_motion(left+1,left)
        self.allow_motion(left-1,left)
        self.allow_motion(right+1,right)
        self.allow_motion(right-1,right)

        #add the cohesin back on in a random position
        self.add_cohesin(with_index = to_move)


    def makeMove(self,rate_forward,rate_reverse,rate_falloff,rtotal=None):
        """Make Gillespie move

        Randomly chooses which move to make in proportion to it's rate and
        makes that move.
        Args:
            rate_forward (float): Rate at which cohesin feet step forward.
            rate_reverse (float): Rate at which cohesin feet step backward.
            rate_falloff (float): Rate at which cohesin fall off.
            rtotal (float, optional): Total rate
        """
        if not rtotal:
            rtotal = self.get_total_rate(rate_forward,rate_reverse,rate_falloff)
        rx=random.random()*rtotal
        total_forward = len(self.forward_moves)*rate_forward
        total_reverse = len(self.reverse_moves)*rate_reverse
        if rx < total_forward and len(self.forward_moves)>0:
            # move forward
            to_move = random.randint(0,len(self.forward_moves)-1)
            #print("forward move",to_move)
            self.step(to_move)
        elif rx < total_forward + total_reverse and len(self.reverse_moves)>0:
            # move backward
            to_move = random.randint(0,len(self.reverse_moves)-1)
            #print("reverse_move",to_move)
            self.step(to_move,reverse=True)
        elif len(self.cohesins):
            to_move = random.randint(0,len(self.cohesins)-1)
            #print("falloff",to_move)
            self.falloff(to_move)
        else:
            print("Warning, move not taken")


    def run_for_time(self,time,rate_forward,rate_reverse,rate_falloff):
        """Run Gillespi algorithm for given time period

        Args:
            time (float): Time period over which to run algorithm
            rate_forward (float): Rate at which cohesin feet step forward.
            rate_reverse (float): Rate at which cohesin feet step backward.
            rate_falloff (float): Rate at which cohesin fall off.
        """
        t=0.0
        while True:
            total_rate = self.get_total_rate(rate_forward, rate_reverse,
                                             rate_falloff)
            dt = np.random.exponential(1/total_rate)
            if t+dt>time:
                return
            t=t+dt
            self.makeMove(rate_forward, rate_reverse, rate_falloff,
                          rtotal=total_rate)

