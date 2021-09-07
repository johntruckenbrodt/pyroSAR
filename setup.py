from setuptools import setup, find_packages
import os

# Create .pyrosar in HOME - Directory
directory = os.path.join(os.path.expanduser('~'), '.pyrosar')
os.makedirs(directory, exist_ok=True)

directory = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='pyroSAR',
      packages=find_packages(),
      include_package_data=True,
      python_requires='>=3',
      setup_requires=['setuptools_scm'],
      use_scm_version=True,
      description='a framework for large-scale SAR satellite data processing',
      classifiers=[
          'License :: OSI Approved :: MIT License',
          'Operating System :: Microsoft :: Windows',
          'Operating System :: POSIX :: Linux',
          'Programming Language :: Python :: 3'
      ],
      install_requires=['progressbar2',
                        'pathos>=0.2',
                        'numpy',
                        'spatialist>=0.8',
                        'pyyaml',
                        'requests',
                        'psycopg2',
                        'SQLAlchemy<1.4',
                        'SQLAlchemy-Utils',
                        'GeoAlchemy2'],
      extras_require={
          'docs': ['sphinx', 'sphinxcontrib-bibtex'],
      },
      url='https://github.com/johntruckenbrodt/pyroSAR.git',
      author='John Truckenbrodt',
      author_email='john.truckenbrodt@uni-jena.de',
      license='MIT',
      zip_safe=False,
      long_description=long_description,
      long_description_content_type='text/markdown')
