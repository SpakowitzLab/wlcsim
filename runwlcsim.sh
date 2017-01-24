#/bin/bash

rm data/*
rm wlcsim.exe
make
./wlcsim.exe -i input/params
