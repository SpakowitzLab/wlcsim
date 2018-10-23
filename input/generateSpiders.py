from spiderMove import *

def get_poly_lengths(file_name=None):
    if file_name == None:
        file_name = "polyLengths"
    lengthSet=[]
    with open(file_name) as lengthsFile:
        for line in lengthsFile:
            lengthSet.append(int(line))
    firstBeads=[1]
    npoly = len(lengthSet)
    for ii in range(0,npoly):
        firstBeads.append(firstBeads[-1]+lengthSet[ii])
    return [lengthSet,firstBeads,npoly]

def combine_bindpair_files(npoly,output="bindpairs"):
    all_bindpair_file = open(output,"w")
    lineNumber=0
    for ii in range(0,npoly):
        firstLine = lineNumber
        with open("bindpairs_chrom"+str(ii)) as fromfile:
            for line in fromfile:
                if int(line) == -1:
                    print(line,file=all_bindpair_file,end='')
                else:
                    print(int(line)+firstLine,file=all_bindpair_file)
                lineNumber=lineNumber+1
    all_bindpair_file.close()

def combine_spider_files(npoly,n_different_lengths,outName):
    # Combine files
    all_spider_file = open(outName,"w")
    print(totalSpiders,file=all_spider_file)
    for ii in range(0,npoly):
        for jj in range(0,n_different_lengths):
            with open("_temp_"+str(ii)+"_"+str(jj)) as fromfile:
                for line in fromfile:
                    print(line,file=all_spider_file,end='')
    all_spider_file.close()

#leglen=[10,25]
leglen=[2,5]
outName = "spiders"


#bindpairs = np.loadtxt(fname).astype(int)
#bindpairs[bindpairs != -1] = bindpairs[bindpairs != -1] - 1
#spiders = makeSpiders(bindpairs,leglength=leglen, moveable_loop=leglen*5)
#print_spiders_for_fortran(spiders,filename=outName,offset=1)


[lengthSet, firstBeads, npoly] = get_poly_lengths()
n_different_lengths = len(leglen)

# Generate spearate file of spiders for each polymer
totalSpiders = 0
for ii in range(0,npoly):
    for jj in range(0,n_different_lengths):
        bindpairs = np.loadtxt("bindpairs_chrom"+str(ii)).astype(int)
        bindpairs[bindpairs != -1] = bindpairs[bindpairs != -1] - 1
        spiders = makeSpiders(bindpairs, leglength=leglen[jj],
                              moveable_loop=leglen[jj]*5)
        totalSpiders=totalSpiders+len(spiders)
        print_spiders_for_fortran(spiders,
                                  filename="_temp_"+str(ii)+"_"+str(jj),
                              offset=firstBeads[ii], printNumber=False)

combine_bindpair_files(npoly)
combine_spider_files(npoly,n_different_lengths,outName)


