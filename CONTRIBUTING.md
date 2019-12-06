# Contributing to pyroSAR

First off, thanks for considering a contribution to pyroSAR. Any contribution, may it be a feature suggestion, a pull 
request or s simple bug report, is valuable to the project and very welcome.
This document is intended as a guideline on best practices.

## How to open an issue
The easiest way to contribute to pyroSAR is by opening an issue. This is intended for reporting software bugs and 
suggesting new features. Before you do, please read through the list of 
[open issues](https://github.com/johntruckenbrodt/pyroSAR/issues) to see whether this issue has already been raised.
This way, duplicates can be reduced and it is easier for the developers to address them.
If you are not sure whether your issue is a duplicate of an existing one, just open a new issue. It is easier to link 
two existing similar issues than separating two different ones contained in one.
For reporting bugs please fill out the template, which is available once you open it. For suggesting new features you
can just delete the template text.  
The following questions need to be answered so that is is possible for the developers to start fixing the software:
- which operating system are you using?  
e.g. Windows 10, Ubuntu 18.4, etc.
- which environment is pyroSAR running in?  
e.g. system-wide Python installation, Anaconda environment, virtual environment, etc.
- which version of pyroSAR are you using?  
one installed via pip or a clone of the GitHub repository?
-  which function of pyroSAR did you call with which parameters?  
- if applicable, which version of SNAP or GAMMA are you using in pyroSAR?
- the full error message

This way the error is reproducible and can quickly be fixed.

## Checking pyroSAR's version
The used version can be obtained like this:
```python
import pyroSAR
print(pyroSAR.__version__)
```
Depending on how you installed pyroSAR the version might look differently. 
If installed via pip with `pip install pyroSAR`, the package is downloaded from 
[PyPI](https://pypi.org/project/pyroSAR/), 
where only the main releases are stored and versions are named e.g. `0.9.1`. 
These can also be found on GitHub [here](https://github.com/johntruckenbrodt/pyroSAR/releases).
If you have installed pyroSAR directly from GitHub like so:
```shell script
python3 -m pip install git+https://github.com/johntruckenbrodt/pyroSAR
```
or have directly cloned a branch from GitHub, your version might look like this:
`0.9.2.dev103+g57eeb30`, in which this naming pattern is used:  
`{next_version}.dev{distance}+{scm letter}{revision hash}`.
In this case we can see that git is used as scm and the latest commit of the software was 
[57eeb30](https://github.com/johntruckenbrodt/pyroSAR/commit/57eeb30970dc6adfee62ca12fd8c8818ecaf3a14), 
which, at the time of checking the version, had a distance of 103 commits to the latest commit.
See [here](https://www.diycode.cc/projects/pypa/setuptools_scm) for more details.
