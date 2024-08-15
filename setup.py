from setuptools import setup

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
   name='orcatools',
   version='1.0',
   description='Python tools for ORCA ',
   author='Patrick Melix',
   author_email='patrick.melix@uni-leipzig.de',
   url="https://github.com/Tonner-Zech-Group/ORCA-tools",
   packages=['orcatools'],
   install_requires=['ase', 'numpy', 'matplotlib']
)