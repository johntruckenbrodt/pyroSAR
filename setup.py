from setuptools import setup, find_packages
import platform
import zipfile as zf
import os

# Create .pyrosar in HOME - Directory
directory = os.path.join(os.path.expanduser('~'), '.pyrosar')

if not os.path.exists(directory):
    os.makedirs(directory)

if platform.system() is 'Windows':
    package_data = {'pyroSAR': ['pkgs/mod_spatialite/*']}
else:
    package_data = {}

setup(name='pyroSAR',
      packages=find_packages(),
      include_package_data=True,
      package_data=package_data,
      version='0.4',
      description='a framework for large-scale SAR satellite data processing',
      classifiers=[
          'Programming Language :: Python',
      ],
      install_requires=['progressbar2',
                        'pathos>=0.2',
                        'numpy',
                        'scoop'],
      url='https://github.com/johntruckenbrodt/pyroSAR.git',
      author='John Truckenbrodt',
      author_email='john.truckenbrodt@uni-jena.de',
      license='MIT',
      zip_safe=False)

if platform.system() is 'Windows':
    subdir = os.path.join(directory, 'mod_spatialite')
    if not os.path.isdir(subdir):
        os.makedirs(subdir)
    mod_spatialite = os.path.join(subdir, 'mod_spatialite.dll')
    if not os.path.isfile(mod_spatialite):
        source_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'pyroSAR', 'pkgs', 'mod_spatialite')
        suffix = 'amd64' if platform.machine().endswith('64') else 'x86'
        source = os.path.join(source_dir, 'mod_spatialite-4.3.0a-win-{}.zip'.format(suffix))
        archive = zf.ZipFile(source, 'r')
        archive.extractall(subdir)
        archive.close()
