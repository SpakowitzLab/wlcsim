
fname = 'names'
print('#/bin/bash')
files=set()
with open(fname) as f:
    for line in f:
        line = line.rstrip()
        line = line.lstrip()
        files.add(line)

for line in files:
    print('python unstructure.py '+str(line))
    print('sed -i -e "s/wlc_d%/wlc__/g" '+str(line))
    print('sed -i -e "s/, wlc_d//g" '+str(line))
    print('sed -i -e "s/,wlc_d//g" '+str(line))
    print('sed -i -e "s/wlc_d,//g" '+str(line))
    print('sed -i -e "s/wlc_d ,//g" '+str(line))
    print('sed -i -e "s/wlc_d//g" '+str(line))
    print('sed -i -e "s/wlc__/wlc_/g" '+str(line))
    print('sed -i -e "s/wlcsim_data,//g" '+str(line))
    print('sed -i -e "s/wlcsim_data ,//g" '+str(line))
    print('sed -i -e "s/,wlcsim_data//g" '+str(line))
    print('sed -i -e "s/, wlcsim_data//g" '+str(line))
    print('sed -i -e "s/use params, only: wlcsim_data//g" '+str(line))


