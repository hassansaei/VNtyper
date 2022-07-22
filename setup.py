from setuptools import setup


setup(name='VNtyper',
      version='1.0.0',
      description='A tool for genotyping MUC1 coding-VNTR from short-read sequencing data',
      author='Hassan Saei',
      author_email='hassan.saei@inserm.fr',
      license='X',
      url='https://github.com/hassansaei/VNtyper',
      test_suite='tests',
      install_requires=['pandas', 'numpy', 'biopython', 'regex', 'argparse', 'PyVCF', ''])

