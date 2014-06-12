import os

from setuptools import setup, find_packages

def read(*paths):
    """Build a file path from *paths* and return the contents."""
    with open(os.path.join(*paths), 'r') as f:
        return f.read()

setup(
    name='smashbenchmarking',
    version='1.0.1',
    packages=find_packages(exclude=["tests*","scripts*"]),
    description='Check the accuracy of one VCF callset against another',
    long_description=(read('README.md')),
    url='http://github.com/amplab/smash/',
    license='BSD',
    author='AMPlab, UC Berkeley',
    author_email='smashbenchmarking@googlegroups.com',
    py_modules=['smashbenchmarking'],
    install_requires=["pyvcf","pyfasta","numpy"],
    include_package_data=True,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)