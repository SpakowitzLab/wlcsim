import textwrap
import seaborn as sns
from polymerMole.glob import *
def makepmlFile(OutName="out.png", colorOption="H3K9me3",
                Ncolors=default_Ncolors, closePymol=True, ball_radius=0.085,
                color_palette='hls', stick_radius=0.05, view = "cube",
                highlightPolymers = None, boundary_thickness=0.15, **kwargs):
    """ Make a pml comand file for pymol and run pymol to generate a png.

    Args:
        OutName (str): Output file name
        clorOptions (str): One of "H3K9me3", "Aseries"
        Ncolors (int): Nuber of colors
        ball_radius (float): radious of balls
        color_palette (str): seaborn color option.  e.t. "hls" or "coolwarm"
        boundary_thickness (float): boudnary wire radious
    """
    myfile = open("autoGen.pml","w")
    # General Head lines
    myfile.write(textwrap.dedent("""
    delete *
    load temp.pdb,snap
    set connect_mode,1

    hide all
    bg_color white
    """))
    if (colorOption == "highlight_homopoly"):
        myfile.write(textwrap.dedent("""
        set_color color1, [1.0, 0.0, 0.0]
        set_color color2, [0.0, 0.0, 1.0]

        show spheres, (name A1)
        color color1, (name A1)
        alter (name A1),vdw=""" +str(ball_radius) + """
        show sticks, (name A1)
        set_bond stick_radius,""" +str(stick_radius) + """, (name A1)
        hide lines, (name A1)

        show spheres, (name A2)
        color color2, (name A2)
        alter (name A2), vdw=""" +str(ball_radius) + """
        show sticks, (name A2)
        set_bond stick_radius, """ +str(stick_radius) + """, (name A2)
        hide lines, (name A2)

        set_bond stick_radius,""" +str(stick_radius) + """ , (name A2), (name A1)

        show spheres, (name C1)
        color black,(name C1)
        alter (name C1),vdw=0.15
        show sticks, (name C1)
        set_bond stick_radius, 0.05, (name C1)
        hide lines, (name C1)

        set sphere_transparency=0.0, (name A2)
        set_bond stick_transparency, 0.00, (name A2)
        """))
    elif (colorOption=="H3K9me3"):
        myfile.write(textwrap.dedent("""
        set_color low_H3K9me3, [0.5, 1.0, 0.917]
        set_color med_H3K9me3, [0.808, 0.596, 0.145]
        set_color high_H3K9me3, [0.417, 0.0, 0.5]

        show spheres, (name A5)
        color red,(name A5)
        alter (name A1),vdw=0.75

        show spheres, (name A1)
        color low_H3K9me3,(name A1)
        alter (name A1),vdw=0.174
        show sticks, (name A1)
        set_bond stick_radius, 0.05, (name A1)
        hide lines, (name A1)

        show spheres, (name A2)
        color med_H3K9me3,(name A2)
        alter (name A2),vdw=0.178
        show sticks, (name A2)
        set_bond stick_radius, 0.05, (name A2)
        hide lines, (name A2)

        show spheres, (name A3)
        color high_H3K9me3,(name A3)
        alter (name A3),vdw=0.178
        show sticks, (name A3)
        set_bond stick_radius, 0.05, (name A3)
        hide lines, (name A3)

        set_bond stick_radius, 0.05, (name A2), (name A1)
        set_bond stick_radius, 0.05, (name A2), (name A3)
        set_bond stick_radius, 0.05, (name A3), (name A1)

                                     
        show spheres, (name HL8)
        color black, (name HL8)
        alter (name HL8),vdw=0.5
        hide lines, (name HL8)

        show spheres, (name C1)
        color black,(name C1)
        alter (name C1),vdw=0.15
        show sticks, (name C1)
        set_bond stick_radius, 0.05, (name C1)
        hide lines, (name C1)
        set_bond stick_radius, 0.05, (name C1), (name A1)
        set_bond stick_radius, 0.05, (name C2), (name A1)
        set_bond stick_radius, 0.05, (name C3), (name A1)

        """))
    elif(colorOption=="Aseries"):
        # Choose from seaborn color palette
        colors = sns.color_palette(color_palette, Ncolors)
        for n in range(Ncolors):
            color = "["+ str(colors[n][0]) + ", " + str(colors[n][1]) + ", " +\
            str(colors[n][2]) + "] "
            color_name = "col"+str(n)
            name ="(name A"+str(n)+")"
            myfile.write("show spheres, "+name+"\n")
            myfile.write("set_color "+color_name+", "+color+"\n")
            myfile.write("color "+color_name+", "+name+"\n")
            myfile.write("alter "+name+", vdw="+str(ball_radius)+"\n")
            myfile.write("show sticks, "+name+"\n")
            myfile.write("hide lines, "+name+"\n")
            #myfile.write("set sphere_transparency, 0.20, "+name+"\n")
            myfile.write("set_bond stick_radius, " + str(stick_radius) + ", "+name+"\n")
            myfile.write("set_bond stick_radius, " + str(stick_radius) + ", "+name+", (name C1)\n")
            #myfile.write("set sphere_transparency=0.9, "+name+"\n")
            if highlightPolymers is not None:
                if n not in highlightPolymers:
                    myfile.write("set sphere_transparency=0.8, "+name+"\n")
                    myfile.write("set_bond stick_transparency, 0.80, "+name+"\n")
            myfile.write("\n\n")

        #for nn in range(Ncolors-1):
        #    for mm in range(1,Ncolors):
        #        if (False and nn != mm-1):
        #            continue
        #            # Use this if you want to shorten pml file and every
        #            # combination isn't needed
        #        name1 ="(name A"+str(nn)+")"
        #        name2 ="(name A"+str(mm)+")"
        #        myfile.write("set_bond stick_radius, " + str(stick_radius) + ", " + name +
        #                     ", " + name2 + "\n")

        myfile.write("set stick_radius, "+str(stick_radius))


        myfile.write(textwrap.dedent("""

        show spheres, (name C1)
        color black,(name C1)
        alter (name C1),vdw=0.15
        show sticks, (name C1)
        set_bond stick_radius, 0.05, (name C1)
        hide lines, (name C1)

        """))

    #Set view
    myfile.write(textwrap.dedent("""

    show (name BLCK)
    alter (name BLCK),vdw=0.02
    color black, (name BLCK)
    show sticks, (name BLCK)
    set_bond stick_radius, """+str(boundary_thickness)+""", (name BLCK)

    """))
    if (view == "single_chromosome_sphere"):
        myfile.write(textwrap.dedent("""
        set_view (\\
             0.000000000,    0.000000000,   -1.000000000,\\
             1.000000000,    0.000000000,   -0.000000000,\\
             0.000000000,   -1.000000000,   -0.000000000,\\
             0.000000000,    0.000000000, -186.786605835,\\
            28.395023346,   30.000000000,   32.000000000,\\
           -21.213415146,  415.786560059,  -20.000000000 )
        """))
    if (view == "chrom10"):
        # for several chromosomes coarse grained
        myfile.write(textwrap.dedent("""
        set_view (\\
             0.000000000,    0.000000000,   -1.000000000,\\
             1.000000000,    0.000000000,   -0.000000000,\\
             0.000000000,   -1.000000000,   -0.000000000,\\
             0.000000522,   -0.000000611, -128.224334717,\\
            28.395023346,   21.675159454,   21.032266617,\\
           -79.775764465,  357.224304199,  -20.000000000 )
        """))
    if (view == "cube"):
        #cube
        myfile.write(textwrap.dedent("""
        set_view (\\
             0.979412198,   -0.005551211,    0.201794311,\\
            -0.201868802,   -0.031240445,    0.978914082,\\
             0.000870125,   -0.999496341,   -0.031717874,\\
             0.000002684,    0.000009418, -260.142150879,\\
            30.736551285,   21.463550568,   30.930513382,\\
            52.141902924,  489.141906738,  -20.000000000 )
        """))

    if (view == "pablo"):
        #cube
        myfile.write(textwrap.dedent("""
        set_view (\
            -0.157233432,    0.314789474,   -0.936047614,\
             0.986173868,   -0.000170971,   -0.165710956,\
            -0.052324444,   -0.949160516,   -0.310410380,\
            -0.000010282,    0.000014623,  -58.981529236,\
            26.081180573,    8.462759018,   11.506986618,\
          -149.017959595,  287.981750488,  -20.000000000 )
        """))

    if (view == "chrom22"):
        #cube
        myfile.write(textwrap.dedent("""
        set_view (\
             0.000000000,    0.000000000,   -1.000000000,\
             1.000000000,    0.000000000,   -0.000000000,\
             0.000000000,   -1.000000000,   -0.000000000,\
            -0.000000700,    0.000000581, -117.193984985,\
            28.395023346,   17.435894012,   17.169824600,\
           -90.806114197,  346.193969727,  -20.000000000 )
        """))

    if (view == "Luke"):
        #cube
        myfile.write(textwrap.dedent("""
        set_view (\
             0.005956888,   -0.997197211,   -0.074564442,\
            -0.997922063,   -0.010712218,    0.063549034,\
            -0.064170003,    0.074030638,   -0.995188534,\
            -0.000077945,    0.000045210, -311.804016113,\
            41.645568848,   36.525054932,   29.839927673,\
           103.808670044,  540.808715820,  -20.000000000 )
        """))
    if (view =="zoom_boundary"):
        myfile.write(textwrap.dedent("""
        set_view (\
             0.979412198,   -0.005551211,    0.201794311,\
            -0.201868802,   -0.031240445,    0.978914082,\
             0.000870125,   -0.999496341,   -0.031717874,\
             0.000003330,    0.000002887,  -74.053474426,\
             9.356276512,   25.627843857,   23.458267212,\
          -133.947280884,  303.052459717,  -20.000000000 )
        """))


    myfile.write("ray\n")


    if (closePymol):
        myfile.write("png "+OutName+", dpi = 100, width=1000, height=1000, ray=1\n")
        myfile.write("quit")
        myfile.close()

