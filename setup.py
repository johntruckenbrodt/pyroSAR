from setuptools import setup, find_packages

setup(name='pyroSAR',
      packages=find_packages(),
      include_package_data=True,
      version='0.1',
      description='a framework for large-scale SAR satellite data processing',
      classifiers=[
          'Programming Language :: Python :: 2.7',
      ],
      install_requires=['progressbar==2.3'],
      url='https://github.com/johntruckenbrodt/pyroSAR.git',
      author='John Truckenbrodt',
      author_email='john.truckenbrodt@uni-jena.de',
      license='MIT',
      zip_safe=False)
