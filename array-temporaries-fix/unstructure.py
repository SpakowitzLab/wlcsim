import sys
import numpy as np
import re
if len(sys.argv) != 2:
      sys.exit('wrong number of arguments')
fname= sys.argv[1]
#fnameMeth= sys.argv[2]


subRutName='topStuff'
subRuts = {subRutName:set()}
whiteSpaceValues = {}
with open(fname) as f:
    for line in f:
        nWhiteSpace=len(line)
        line=line.lstrip()
        nWhiteSpace=nWhiteSpace-len(line)
        if (line[0:10].lower() == "subroutine" or line[0:8].lower() ==
            "function"):
            subRutName = re.findall(r"[\w']+",line)[1]
            subRuts[subRutName]=set()
            whiteSpaceValues[subRutName]=nWhiteSpace

        words = re.findall(r"[\w'%]+",line)
        for word in words:
            if str(word).lower()[:6] == "wlc_d%":
                subRuts[subRutName].add(word[6:])

#for key in subRuts.keys():
#    print('in '+key+' you have')
#    print(subRuts[key])


f = open(fname,"r")
contents = f.readlines()
f.close()


f = open(fname,"w")

for line in contents:
    if (line.lstrip().lower()[0:17] == "type(wlcsim_data)"):
        continue
    print(line, end="",file=f)
    line=line.lstrip()
    if (line[0:10].lower() == "subroutine" or line[0:8].lower() ==
            "function"):
        subRutName = re.findall(r"[\w']+",line)[1]
        if len(subRuts[subRutName]) == 0:
            continue
        for ii in range(whiteSpaceValues[subRutName]):
            print(" ",end="",file=f)
        print("! values from wlcsim_data",file=f)

        for ii in range(whiteSpaceValues[subRutName]):
            print(" ",end="",file=f)
        print("use params, only: ",end="",file=f)

        count = 0
        camma = False
        for name in subRuts[subRutName]:
            if count>4:
                print("&\n",end="",file=f)
                for ii in range(whiteSpaceValues[subRutName]+4):
                    print(" ",end="",file=f)
                count=0
            else:
                count=count+1
            if camma:
                print(", ",end="",file=f)
            camma=True
            print("wlc__"+name,end="",file=f)
        print("",file=f)

f.close()
