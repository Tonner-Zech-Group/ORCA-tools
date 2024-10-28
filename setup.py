from setuptools import setup

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
   name='orcatools',
   version='1.0',
   description='Python tools for ORCA',
   long_description=open('README.md').read(),
   long_description_content_type='text/markdown',
   author='Patrick Melix',
   author_email='patrick.melix@uni-leipzig.de',
   url="https://github.com/Tonner-Zech-Group/ORCA-tools",
   packages=['orcatools'],
   install_requires=['ase', 'numpy', 'matplotlib'],
   dependency_links = ['https://github.com/virtualzx-nad/geodesic-interpolate.git']
)