from setuptools import setup
import numpy

setup(name='VNtyper',
      version='1.0.0',
      description='A tool for genotyping MUC1 coding-VNTR from short-read sequencing data',
      author='Hassan Saei',
      author_email='hassan.saei@inserm.fr',
      license='X',
      url='https://github.com/mehrdadbakhtiari/adVNTR',
      test_suite='tests',
      packages=['advntr', 'pomegranate'],
      package_dir={'advntr': 'advntr', 'advntr.pomegranate': 'pomegranate'},
      install_requires=['networkx==1.11', 'scipy', 'biopython', 'cython', 'scikit-learn'],
      provides=["advntr"],
      entry_points={
            'console_scripts': ['advntr=advntr.__main__:main']
      },
      ext_modules=cythonize(["pomegranate/*.pyx"]),
      include_dirs=[numpy.get_include()],
      classifiers=["Environment :: Console",
                   "Intended Audience :: Developers",
                   "Intended Audience :: Science/Research",
                   "Operating System :: Unix",
                   "Programming Language :: Python",
                   "Programming Language :: Python :: 2",
                   "Programming Language :: Python :: 3",
                   "Topic :: Scientific/Engineering :: Bio-Informatics"],
      )
