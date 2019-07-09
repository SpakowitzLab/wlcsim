from polymerMole.r2pdb import r2pdb
from polymerMole.makepmlFile import makepmlFile
import os
direct = "../"
baseNames = {direct: len(open(direct+'data/energies').readlines())-2 }
kwargs ={}
for baseName in baseNames.keys():
    savept_max = baseNames[baseName]
    for savept in [savept_max]: #[2,3,4,7,8,9,12]: #[savept_max]: #[0,5,9]: #range(0,111,20):

        for rep in [0]: #range(1,16): #range(1,8): #range(2,10): #[2,4,6,7,8,9,10,11,12,14]:
            if False:
                suffix = "v"+str(rep)
            else:
                suffix = ""

            xyzFileName=baseName+"data/r"+str(savept)+suffix
            if baseName in ["../../multiTenkb_woCTCF/Tenkb_woCTCF_2_cupn0p0001_chin240/"]:
                continue
            kwargs['skip'] = 1

            kwargs['closePymol']=True
            image = "Cube_full_CTCFs"
            image = "chrom22"
            image = "EpiColor22"
            image = "Cube_full"
            image = "Luke"
            image = "bind2boundary"
            image = "single_polymer"
            image = "Cube_recolor"
            image = "Cube"

            if (image=="single_polymer"):
                kwargs['skip']=1
                kwargs['Ncolors'] = 10
                kwargs['color_type'] = "sequential"
                kwargs['colorOption'] = "Aseries"
                kwargs['maxpoints'] = 99999
                kwargs['circles'] = None
                kwargs['xlimits'] = None
                kwargs['ylimits'] = None
                kwargs['zlimits'] = None
                kwargs['scalebar'] = None
                kwargs['methFileName'] = None
                kwargs['color_cohisn']=False
                kwargs['bindFileName']=None
                kwargs['color_palette']=None
                kwargs['ball_radius'] = 2.0
                kwargs['stick_radius']=1.0
                kwargs['view']="Luke"
                kwargs['polymerLengthFile'] = None
                kwargs['recenter'] = True

            if (image == "bind2boundary"):
                kwargs['skip']=1
                kwargs['Ncolors'] = None
                kwargs['color_type'] = "meth"
                kwargs['colorOption'] = "H3K9me3"
                kwargs['circles'] = None
                kwargs['scalebar']=100/28.7
                kwargs['methFileName'] = baseName+"input/meth"
                kwargs['color_cohisn']=False
                kwargs['bindFileName']=None
                kwargs['color_palette']=None
                kwargs['ball_radius'] = 0.198
                kwargs['stick_radius']=0.05
                kwargs['polymerLengthFile'] = None
                kwargs['highlight_file'] = baseName+"data/bind_points"
                kwargs['ylimits'] = [26,38]
                zoom = False
                if zoom:
                    kwargs['xlimits'] = [0,15]
                    kwargs['zlimits'] = [16,31]
                    kwargs['view']="zoom_boundary"
                    kwargs['cube'] = [[0,26,16],[15,38,31]]
                    kwargs['boundary_thickness'] = 0.05
                else:
                    kwargs['view']="cube"
                    kwargs['cube'] = {'1':[0,64],'2':[[0,26,16],[15,38,31]]}

            if (image=="Luke"):
                kwargs['skip']=1
                kwargs['Ncolors'] = None
                kwargs['color_type'] = "firstFraction"
                kwargs['colorOption'] = "highlight_homopoly"
                kwargs['maxpoints'] = 5*3000
                kwargs['fractionType1'] = 0.0
                kwargs['circles'] = None
                kwargs['xlimits'] = None
                kwargs['ylimits'] = None
                kwargs['zlimits'] = None
                kwargs['scalebar'] = None
                kwargs['methFileName'] = None
                kwargs['color_cohisn']=False
                kwargs['bindFileName']=None
                kwargs['color_palette']=None
                enlarge = 2
                kwargs['ball_radius'] = 0.048*enlarge
                kwargs['stick_radius']=0.048*enlarge
                kwargs['cube']=[[0,0,0],[80,80,3]]
                kwargs['view']="Luke"
                kwargs['polymerLengthFile'] = 5
                kwargs['period'] = (80.0, 80.0, None)

            if (image=="CTCFs_in_sphere"):
                kwargs['skip']=5
                kwargs['Ncolors'] = None
                kwargs['color_type'] = "meth"
                kwargs['colorOption'] = "H3K9me3"
                kwargs['circles'] = [(31,[32.0,32.0,32.0])]
                kwargs['xlimits'] = None
                #kwargs['scalebar'] = None
                kwargs['methFileName'] = baseName+"input/meth"
                kwargs['color_cohisn']=True
                kwargs['bindFileName']=baseName+"input/bindpairs"
                kwargs['color_palette']=None
                kwargs['ball_radius'] = 0.198
                kwargs['stick_radius']=0.05
                #kwargs['cube']=[1,63]
                kwargs['view']="cube"
                kwargs['polymerLengthFile'] = None
                kwargs['ylimits'] = None

            if (image=="chrom22"): #color polymers
                kwargs['polymerLengthFile'] = baseName+"input/polyLengths"
                kwargs['Ncolors'] = sum(1 for line in
                                        open(kwargs['polymerLengthFile']))
                kwargs['color_type'] = "polymer"
                kwargs['colorOption']="Aseries"
                kwargs['circles'] = [(15,[16.0,16.0,16.0])]
                kwargs['xlimits']=[16, 20]
                #kwargs['scalebar']=500/100
                kwargs['methFileName'] = None
                kwargs['bindFileName'] = None
                kwargs['color_cohisn'] = False
                kwargs['color_palette']="hls"
                kwargs['ball_radius']=0.07
                kwargs['stick_radius']=0.07
                kwargs['view']="chrom22"
                kwargs['highlightPolymers'] = [1, 15, 18]

            if (image == "EpiColor22"): #
                kwargs['polymerLengthFile'] = baseName+"input/polyLengths"
                kwargs['Ncolors'] = 50
                kwargs['color_type'] = "meth10"
                kwargs['colorOption']="Aseries"
                kwargs['circles'] = [(15,[16.0,16.0,16.0])]
                kwargs['xlimits']=[16, 20]
                #kwargs['scalebar']=500/100
                kwargs['methFileName'] = baseName+"input/ab"
                kwargs['bindFileName'] = baseName+"input/bindpairs"
                kwargs['color_cohisn'] = False
                kwargs['color_palette']="coolwarm"
                kwargs['ball_radius']=0.07
                kwargs['stick_radius']=0.07
                kwargs['view']="chrom22"

            if (image=="New"):
                kwargs['skip']=1
                kwargs['Ncolors'] = None
                kwargs['color_type'] = "meth"
                kwargs['colorOption']= "H3K9me3"
                kwargs['circles'] = None
                kwargs['xlimits'] = [0,15]
                kwargs['scalebar'] = None
                kwargs['methFileName'] = baseName+"input/meth"
                kwargs['color_cohisn']=False
                kwargs['bindFileName']=None
                kwargs['color_palette']=None
                kwargs['ball_radius'] = 0.198
                kwargs['stick_radius']=0.05
                kwargs['cube']=[1,63]
                kwargs['view']="cube"
                kwargs['polymerLengthFile'] = None
                kwargs['ylimits'] = None


            if (image=="Cube"): # Cube
                kwargs['skip']=1
                kwargs['Ncolors'] = None
                kwargs['color_type'] = "meth"
                kwargs['colorOption'] = "H3K9me3"
                kwargs['circles'] = None
                kwargs['xlimits'] = None
                kwargs['scalebar']=100/28.7
                kwargs['methFileName'] = baseName+"input/meth"
                kwargs['color_cohisn']=False
                kwargs['bindFileName']=None
                kwargs['color_palette']=None
                kwargs['ball_radius'] = 0.198
                kwargs['stick_radius']=0.05
                kwargs['cube']=[0,64]
                kwargs['view']="cube"
                kwargs['polymerLengthFile'] = None
                kwargs['ylimits'] = [26,38] # what I use most of the time
                #kwargs['ylimits'] = [29,35] # for when I run out of beads
                #kwargs['zlimits'] = [0,32]
                #kwargs['filter_meth'] = 'PNAS_window'
                #kwargs['mirror'] = [None, None, 32]

            if (image=="Cube_recolor"): # Cube
                kwargs['skip']=1
                kwargs['Ncolors'] = None
                kwargs['color_type'] = "meth"
                kwargs['colorOption'] = "H3K9me3"
                kwargs['circles'] = None
                kwargs['xlimits'] = None
                #kwargs['scalebar']=100/28.7
                kwargs['methFileName'] = baseName+"input/meth"
                kwargs['color_cohisn']=False
                kwargs['bindFileName']=None
                kwargs['color_palette']=None
                kwargs['ball_radius'] = 0.198
                kwargs['stick_radius']=0.05
                kwargs['cube']=[0,64]
                kwargs['view']="cube"
                kwargs['polymerLengthFile'] = None
                c=10
                kwargs['ylimits'] = [24+c,38+c] # what I use most of the time
                kwargs['filter_meth'] = 'PNAS_window'

            if (image=="Cube_full"): # Cube
                kwargs['skip']=5
                kwargs['Ncolors'] = None
                kwargs['color_type'] = "meth"
                kwargs['colorOption'] = "H3K9me3"
                kwargs['circles'] = None
                kwargs['xlimits'] = None
                kwargs['scalebar'] = None
                kwargs['methFileName'] = baseName+"input/meth"
                kwargs['color_cohisn']=False
                kwargs['bindFileName']=None
                kwargs['color_palette']=None
                kwargs['ball_radius'] = 0.198
                kwargs['stick_radius']=0.05
                #cube=[1,63]
                kwargs['cube']=[0,64]
                kwargs['view']="cube"
                kwargs['polymerLengthFile'] = None
                kwargs['ylimits'] = None

            if (image=="Cube_full_CTCFs"): # Cube
                kwargs['skip']=5
                kwargs['Ncolors'] = None
                kwargs['color_type'] = "meth"
                kwargs['colorOption'] = "H3K9me3"
                kwargs['circles'] = None
                kwargs['xlimits'] = None
                kwargs['scalebar'] = None
                kwargs['methFileName'] = baseName+"input/meth"
                kwargs['color_cohisn']=True
                kwargs['bindFileName']=baseName+"input/bindpairs"
                kwargs['color_palette']=None
                kwargs['ball_radius'] = 0.198
                kwargs['stick_radius']=0.05
                kwargs['cube']=[1,63]
                kwargs['view']="cube"
                kwargs['polymerLengthFile'] = None
                kwargs['ylimits'] = None

            if (image=="PolyColor"): #color polymers
                kwargs['polymerLengthFile'] = baseName+"input/polyLengths"
                kwargs['Ncolors'] = sum(1 for line in
                                        open(kwargs['polymerLengthFile']))
                kwargs['color_type'] = "polymer"
                kwargs['colorOption']="Aseries"
                kwargs['circles'] = [(19.5,[20.5,20.5,20.5])]
                kwargs['xlimits']=[19.5,25.0]
                kwargs['scalebar']=500/100
                kwargs['methFileName'] = None
                kwargs['bindFileName'] = None
                kwargs['color_cohisn'] = False
                kwargs['color_palette']="hls"
                kwargs['ball_radius']=0.07
                kwargs['stick_radius']=0.07
                kwargs['view']="chrom10"

            if (image == "EpiColor"): #
                kwargs['Ncolors'] = 13
                kwargs['color_type'] = "meth10"
                kwargs['colorOption']="Aseries"
                kwargs['circles'] = [(19.5,[20.5,20.5,20.5])]
                kwargs['xlimits']=[19.5,25.0]
                kwargs['scalebar']=500/100
                kwargs['methFileName'] = baseName+"input/ab"
                kwargs['bindFileName'] = baseName+"input/bindpairs"
                kwargs['color_cohisn'] = True
                kwargs['color_palette']="coolwarm"
                kwargs['ball_radius']=0.07
                kwargs['stick_radius']=0.07
                kwargs['view']="chrom10"
                kwargs['polymerLengthFile'] = None

            if (image == "tri_color"):
                kwargs['skip']=1
                kwargs['Ncolors'] = None
                kwargs['color_type'] = "meth"
                kwargs['colorOption'] = "H3K9me3"
                kwargs['circles'] = [(31.0,[32.0,32.0,32.0])]
                kwargs['xlimits']=[32,36]
                kwargs['scalebar']=250/28.7
                kwargs['methFileName'] = baseName+"input/meth"
                kwargs['color_cohisn']=False
                kwargs['bindFileName']=None
                kwargs['color_palette']=None
                kwargs['ball_radius'] = 0.198
                kwargs['stick_radius']=0.05
                kwargs['cube']=None
                kwargs['view']="cube"
                kwargs['polymerLengthFile'] = None


                
            kwargs['OutName'] = baseName+"data/"+image+str(savept)+suffix+".png"
            
            r2pdb(xyzFileName, **kwargs)
            makepmlFile(**kwargs)

            os.system("pymol autoGen.pml")
           # try:
           #     r2pdb(xyzFileName, **kwargs)
           #     makepmlFile(**kwargs)

           #     os.system("pymol autoGen.pml")
           # except:
           #     print("cant display " + kwargs['OutName'])

