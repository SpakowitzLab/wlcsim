#/bin/bash
grep -ri "wlcsim_data" src/mc/*.f* | awk -F':' '{print $1}' > names
python commands.py > allcomands.sh
chmod +x allcomands.sh
./allcomands.sh
grep -ri "wlcsim_data" src/wlcsim/*.f* | awk -F':' '{print $1}' > names
python commands.py > allcomands.sh
chmod +x allcomands.sh
./allcomands.sh
sed -i -e "s/max_sselta/maxWlcDelta/g" src/wlcsim/initcond.f03
sed -i -e "s/max_sselta/maxWlcDelta/g" src/wlcsim/params.f03
