from setuptools import setup, find_packages
from codecs import open
from os import path

__version__ = '0.0.1'

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# get the dependencies and installs
with open(path.join(here, 'requirements.txt'), encoding='utf-8') as f:
    all_reqs = f.read().split('\n')
with open(path.join(here, 'dev-requirements.txt'), encoding='utf-8') as f:
    dev_reqs = f.read().split('\n')

install_requires = [x.strip() for x in all_reqs if 'git+' not in x]
dependency_links = [x.strip().replace('git+', '') for x in all_reqs if x.startswith('git+')]
dev_reqs = [x.strip().replace('git+', '') for x in dev_reqs if x.startswith('git+')]

setup(
    name='mapping-graph',
    version=__version__,
    description='Utilities for creating and evaluating mapping graphs',
    long_description=long_description,
    url='https://github.com/ur-whitelab/mapping-graphs',
    download_url='https://github.com/ur-whitelab/mapping-graphs/tarball/' + __version__,
    license='GPL3',
    classifiers=[
      'Development Status :: 3 - Alpha',
      'Intended Audience :: Developers',
      'Programming Language :: Python :: 3',
    ],
    keywords='',
    packages=find_packages(exclude=['docs', 'tests*']),
    include_package_data=True,
    author='Andrew White',
    install_requires=install_requires,
    dependency_links=dependency_links,
    entry_points=
      {
        'console_scripts': ['mapping-graph=mapg.run:start']
      },
    tests_require=dev_reqs,
    author_email='white.d.andrew@gmail.com'
)
