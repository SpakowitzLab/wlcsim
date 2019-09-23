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
          classifiers=["Intended Audience :: Science/Research",
                       'Intended Audience :: Developers',
                       'Development Status :: 2 - Pre-Alpha',
                       'License :: OSI Approved :: MIT License',
                       "Natural Language :: English",
                       "Programming Language :: Python :: 3 :: Only",
                       "Programming Language :: Fortran",
                       "Topic :: Scientific/Engineering :: Chemistry",
                       "Topic :: Scientific/Engineering :: Physics",
                      ],
          packages=find_packages(include=["wlcsim", "wlcsim.*"]),
          package_data={
              'wlcsim.tabulation': ['*.csv'],
          },
          install_requires=["scipy", "statsmodels", "matplotlib", "seaborn",
                            "mpmath", "pandas", "numpy", "bruno_util",
                            "PyQt5"],
          )
