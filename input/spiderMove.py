import numpy as np

def blank_slide_section(center_bead):
    #print('blank section ',center_bead,' created')
    return {'L_slide_bead':center_bead,
            'R_slide_bead':center_bead,
            'L_knee':None,
            'L_toe':None,
            'R_knee':None,
            'R_toe':None}

def combine(subsumed_section_id,subsuming_section_id,sections,membership,
            bindpairs,leglength,moveable_loop):
    #print("subsumed_section_id",subsumed_section_id,"subsuming_section_id",subsuming_section_id)
    subsumed = sections[subsumed_section_id]
    subsuming = sections[subsuming_section_id]

    #print("combining section",print_section(subsumed,True),' into',
    #      print_section(subsuming, True))
    # make entire memebership region subsumed_section_id
    if (subsumed['L_slide_bead'] > subsuming['L_slide_bead']):
        left = subsuming['L_slide_bead']
        right = subsumed['R_slide_bead']
        subsuming['R_slide_bead'] = right
        membership[left:right+1] = subsuming_section_id
        if subsumed['R_knee'] != None: # attached leg from subsumed onto subsuming
            subsuming['R_knee'] = subsumed['R_knee']
            subsuming['R_toe'] = subsumed['R_toe']
        elif subsuming['R_knee'] != None: # grow new leg to replace old one
            subsuming['R_knee'] = None
            subsuming['R_toe'] = None
            growLeg(subsuming_id,bindpairs,membership,sections,
                    leglength=leglength,moveable_loop=moveable_loop,
                    direction=1)
    else:
        right = subsuming['R_slide_bead']
        left = subsumed['L_slide_bead']
        subsuming['L_slide_bead'] = left
        membership[left:right+1] = subsuming_section_id
        if subsumed['L_knee'] != None: # attached leg from subsumed onto subsuming
            subsuming['L_knee'] = subsumed['L_knee']
            subsuming['L_toe'] = subsumed['L_toe']
        elif subsuming['L_knee'] != None: # grow new leg to replace old one
            subsuming['L_knee'] = None
            subsuming['L_toe'] = None
            growLeg(subsuming_id,bindpairs,membership,sections,
                    leglength=leglength,moveable_loop=moveable_loop,
                    direction=-1)

    #print("combined_section",print_section(subsuming,True))
    del sections[subsumed_section_id]
    #print("combining done",membership,"    combinedSection:")
    #print_section(subsuming)

def checkForCollision(bead,membership,section_id,slide_section,sections,
                      polymerLength,bindpairs,leglength,moveable_loop):
    #polymerLength=len(bindpairs)
    if bead<0: # encounter end of polymer
        membership[0:slide_section['L_slide_bead'] ]=section_id
        slide_section['L_slide_bead']=0
        return 'PolymerEnd'
    if bead>=polymerLength: # encounter end of polymer
        membership[slide_section['R_slide_bead']:]=section_id
        slide_section['R_slide_bead']=polymerLength-1
        return 'PolymerEnd'
    if membership[bead]>=0: # combine with another section
        combine(section_id,membership[bead],sections,membership,
                bindpairs,leglength,moveable_loop)
        return 'Combined'
    return False

def print_section(section,return_string=False):
    iostr = (str(section['L_toe'])+'-'+str(section['L_knee'])+'-'
            +str(section['L_slide_bead'])+'='+str(str(section['R_slide_bead']))
            +'-'+str(section['R_knee'])+'-'+str(section['R_toe']))
    if return_string:
        return iostr
    else:
        print(iostr)


def is_independent_polymer(end1,end2,bindpairs):
    low=min(end1,end2)
    high=max(end1,end2)
    if np.min(bindpairs[low:high+1])<low:
        return False
    if np.max(bindpairs[low:high+1])>high:
        return False
    return True


def growLeg(section_id,bindpairs,membership,sections,leglength=10,moveable_loop=50,direction=1):
    slide_section=sections[section_id]
    if sections[section_id][{1:'R_toe',-1:'L_toe'}[direction]] != None:
        return 'LegExists'
    if moveable_loop<4*leglength:
        print('moveable_loop',moveable_loop)
        print('leglength',leglength)
        raise ValueError('moveable_loop must be >= 4*leglength')
    if direction not in [-1,1]:
        raise ValueError('direction must be in [-1,1]')

    LorR = {1:'R_slide_bead',-1:'L_slide_bead'}[direction]

    #print('Growing Leg of section',section_id,' from',slide_section[LorR],{1:'up',-1:'down'}[direction])
    knee=slide_section[LorR]
    for ii in range(0,leglength):
        knee=knee+direction # grow leg by one
        collision = checkForCollision(knee,membership,section_id,slide_section,
                                      sections,len(bindpairs),bindpairs,
                                      leglength,moveable_loop)
        if collision:
            return collision
        otherend = bindpairs[knee]
        if otherend == -1: # normal chain
            continue
        if (abs(otherend-knee) <= moveable_loop
            # loop is short enough
            and is_independent_polymer(knee,otherend,bindpairs)
            # loop isn't interanlly connected elsewhere
            and (otherend-knee)*direction >=0): # loop is in same direction as leg
            knee = otherend
            collision = checkForCollision(knee,membership,section_id,
                                          slide_section,sections,
                                          len(bindpairs),bindpairs,leglength,
                                          moveable_loop)
            if collision:
                return collision
            continue

        # didn't have normal step or subsumed loop
        # extend slide region
        # create new section with two more legs
        membership[min(slide_section[LorR],knee):max(slide_section[LorR],knee)+1] = section_id
        membership[otherend]=otherend
        slide_section[LorR]=knee
        sections[otherend]=blank_slide_section(otherend) # name the section by its starting bead

        whatHappend = growLeg(otherend,bindpairs,membership,sections,
                              leglength=leglength,moveable_loop=moveable_loop,direction=1)
        if whatHappend != 'Combined':
            growLeg(otherend,bindpairs,membership,sections,
                    leglength=leglength,moveable_loop=moveable_loop,direction=-1)

        # the leg we've built so far must be slid instead, continue building from there
        #slide_section[LorR]=knee
        #print('Retry Growing Leg of section',section_id,' from',slide_section[LorR],{1:'up',-1:'down'}[direction])
        return growLeg(section_id,bindpairs,membership,sections,leglength=leglength,moveable_loop=moveable_loop,direction=direction)


    toe=knee
    for ii in range(0,leglength):
        toe=toe+direction # grow leg by one
        collision = checkForCollision(toe,membership,section_id,slide_section,
                                      sections,len(bindpairs),bindpairs,
                                      leglength,moveable_loop)
        if collision:
            return collision
        #print('toe',toe,'membership',membership[toe],'other end',bindpairs[toe])
        otherend = bindpairs[toe]
        if otherend == -1: # normal chain
            continue

        # loop subsumed into leg
        if (abs(otherend-toe) <= moveable_loop
            # loop is short enough
            and is_independent_polymer(toe,otherend,bindpairs)
            # loop isn't interanlly connected elsewhere
            and (otherend-toe)*direction >=0): # loop is in same direction as leg
            toe = otherend
            collision = checkForCollision(toe,membership,section_id,
                                          slide_section,sections,len(bindpairs),bindpairs,
                                          leglength,moveable_loop)
            if collision:
                return collision
            continue

        # didn't have normal step or subsumed loop
        # create new section with two more legs
        membership[min(slide_section[LorR],toe):max(slide_section[LorR],toe)+1] = section_id
        membership[otherend]=otherend
        slide_section[LorR]=toe
        sections[otherend]=blank_slide_section(otherend) # name the section by its starting bead
        whatHappend = growLeg(otherend,bindpairs,membership,sections,
                              leglength=leglength,moveable_loop=moveable_loop,direction=1)
        if whatHappend != 'Combined':
            growLeg(otherend,bindpairs,membership,sections,
                    leglength=leglength,moveable_loop=moveable_loop,direction=-1)

        # the leg we've built so far must be slid instead, continue building from there
        #print('Retry Growing Leg of section',section_id,' from',slide_section[LorR],{1:'up',-1:'down'}[direction])
        return growLeg(section_id,bindpairs,membership,sections,leglength=leglength,moveable_loop=moveable_loop,direction=direction)

    LorRknee = {1:'R_knee',-1:'L_knee'}[direction]
    LorRtoe = {1:'R_toe',-1:'L_toe'}[direction]
    slide_section[LorRknee]=knee
    slide_section[LorRtoe]=toe


    return 'normal'


def makeSections(bindpairs,end1,leglength=5,moveable_loop=25):
    end2=bindpairs[end1]
    if end2==-1:
        raise ValueError('end1 must be bound to something')
    sections={}
    membership=np.zeros(len(bindpairs)).astype(int)-1


    sections[end1] = blank_slide_section(end1)
    membership[end1]=end1
    sections[end2] = blank_slide_section(end2)
    membership[end2]=end2
    keys = list(sections.keys())

    for section_id in keys:
        whatHappend = growLeg(section_id,bindpairs,membership,sections,
                              leglength=leglength,moveable_loop=moveable_loop,direction=1)
        if whatHappend != 'Combined':
            growLeg(section_id,bindpairs,membership,sections,
                    leglength=leglength,moveable_loop=moveable_loop,direction=-1)

    ID_list = np.sort(list(sections.keys()))
    for ii in range(len(ID_list)-1):
        rightID= ID_list[ii+1]
        leftID=ID_list[ii]
        seperation = sections[rightID]['L_slide_bead'] - sections[leftID]['R_slide_bead']
        overlapped_legs = sections[rightID]['L_toe'] <= sections[leftID]['R_toe']
        if seperation < moveable_loop or overlapped_legs:
            combine(leftID,rightID,sections,membership,
                    bindpairs, leglength, moveable_loop)

    for key in sections:
        section=sections[key]
        if section['L_toe']==None:
            raise ValueError('missing left leg')
        if section['R_toe']==None:
            raise ValueError('missing right leg')
        if section['L_toe'] >= section['L_knee']:
            raise ValueError('left knee and toe out of order')
        if section['L_knee'] >= section['L_slide_bead']:
            raise ValueError('left hip and knee out of order')
        if section['R_toe'] <= section['R_knee']:
            raise ValueError('right knee and toe out of order')
        if section['R_knee'] <= section['R_slide_bead']:
            raise ValueError('right hip and knee out of order')
    return sections


def print_sections(sections):
    for ID in sections:
        print_section(sections[ID])


import sys
from contextlib import contextmanager
@contextmanager
def file_or_stdout(file_name):
    if file_name is None:
        yield sys.stdout
    else:
        with open(file_name, 'w') as out_file:
            yield out_file

def print_spiders_for_fortran(spiders,filename=None,offset=0):
    #file = open(filename,"w")
    with file_or_stdout(filename) as file:
        print(len(spiders),file=file)
        for spider in spiders:
            print(len(spider),file=file)
            legs=[]
            moved_regions=[]
            for key in spider:
                section=spider[key]
                moved_retion = [section['L_slide_bead'],section['R_slide_bead']]
                print(section['L_slide_bead']+offset,' ',
                      section['R_slide_bead']+offset,file=file)
                if section['L_toe'] != None:
                    legs.append([section['L_slide_bead'],
                                 section['L_knee'],
                                 section['L_toe']])
                    moved_retion[0]=section['L_toe']
                if section['R_toe'] != None:
                    legs.append([section['R_slide_bead'],
                                 section['R_knee'],
                                 section['R_toe']])
                    moved_retion[1]=section['R_toe']
                moved_regions.append(moved_retion)
            print(len(legs),file=file)
            for leg in legs:
                for val in leg:
                    print(val+offset,end=' ',file=file)
                print('',file=file)
            for moved_retion in moved_regions:
                print(moved_retion[0]+offset,' ',
                      moved_retion[1]+offset,file=file)


def makeSpiders(bindpairs,leglength=5,moveable_loop=25):
    min_slide_loop=moveable_loop
    spiders = []
    for ii in range(0,len(bindpairs)):
        if bindpairs[ii] == -1:
            continue

        if bindpairs[ii] < ii+min_slide_loop:
            continue

        sections = makeSections(bindpairs,ii,leglength=leglength,moveable_loop=moveable_loop)
        spiders.append(sections)
        #print_sections(sections)
        #print('-------')
    return spiders

