from setuptools import setup, find_packages

setup(name='pyroSAR',
      packages=find_packages(),
      version='0.1',
      description='a framework for large-scale SAR satellite data processing',
      classifiers=[
          'Programming Language :: Python :: 2.7',
      ],
      url='https://github.com/johntruckenbrodt/pyroSAR.git',
      author='John Truckenbrodt',
      author_email='john.truckenbrodt@uni-jena.de',
      license='MIT',
      zip_safe=False)
