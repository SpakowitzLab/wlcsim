import sys

from setuptools import setup, find_packages
import versioneer


if __name__ == "__main__":
    if sys.version_info < (3, 5):
        raise ImportError("requires Python>=3.5")

    setup(name="wlcsim",
          author="Bruno Beltran",
          version=versioneer.get_version(),
          cmdclass=versioneer.get_cmdclass(),
          url="https://github.com/SpakowitzLab/BasicWLC",
          license="",
          long_description=open("README.md").read(),
          classifiers=["Environment :: Console",
                       "Intended Audience :: Science/Research",
                       "License :: Fre for non-commercial use",
                       "Natural Language :: English",
                       "Programming Language :: Python :: 3 :: Only",
                       "Programming Language :: Fortran",
                       "Topic :: Scientific/Engineering :: Chemistry",
                       "Topic :: Scientific/Engineering :: Physics",
                       "Topic :: Scientific/Engineering :: Vizualization",
                      ],
          packages=find_packages(include=["wlcsim", "wlcsim.*"]),
          install_requires=["pandas>=0.18.1", "numpy>=1.11.1"],
          )
