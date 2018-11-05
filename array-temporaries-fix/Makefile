all: clean
	 gfortran -ggdb -cpp -fcheck=all -Wall -pedantic -fall-intrinsics -Warray-temporaries -Wno-surprising -Waliasing -Wconversion-extra -Wcompare-reals -Warray-bounds -Wattributes -Wextra mainQuinn.f03

clean: arrayf.f03 mainQuinn.f03
	rm -rf mymod.mod a.out
