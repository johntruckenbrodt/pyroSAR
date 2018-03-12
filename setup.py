from setuptools import setup, find_packages

setup(name='pyroSAR',
      packages=find_packages(),
      include_package_data=True,
      version='0.2',
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
