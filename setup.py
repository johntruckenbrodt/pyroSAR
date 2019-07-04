from setuptools import setup, find_packages
import os
import sys

# Create .pyrosar in HOME - Directory
directory = os.path.join(os.path.expanduser('~'), '.pyrosar')

if not os.path.exists(directory):
    os.makedirs(directory)

directory = os.path.abspath(os.path.dirname(__file__))
if sys.version_info >= (3, 0):
    with open(os.path.join(directory, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()
else:
    with open(os.path.join(directory, 'README.md')) as f:
        long_description = f.read()

setup(name='pyroSAR',
      packages=find_packages(),
      include_package_data=True,
      version='0.9.1',
      description='a framework for large-scale SAR satellite data processing',
      classifiers=[
          'License :: OSI Approved :: MIT License',
          'Operating System :: Microsoft :: Windows',
          'Operating System :: POSIX :: Linux',
          'Programming Language :: Python',
      ],
      install_requires=['progressbar2',
                        'pathos>=0.2',
                        'numpy',
                        'scoop',
                        'spatialist==0.2.9',
                        'pyyaml'],
      extras_require={
          'docs': ['sphinx'],
      },
      url='https://github.com/johntruckenbrodt/pyroSAR.git',
      author='John Truckenbrodt',
      author_email='john.truckenbrodt@uni-jena.de',
      license='MIT',
      zip_safe=False,
      long_description=long_description,
      long_description_content_type='text/markdown')
