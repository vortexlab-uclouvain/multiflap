#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""multiflap - periodic orbit detection and stability assessment"""

from setuptools import setup, find_packages

version = "1.1"

setup(
    name='multiflap',
    version=version,
    author="Gianmarco Ducci",
    author_email="gianmarco.ducci@uclouvain.be",
    description="periodic orbit detection and stability assessment",
    license='LICENSE',
    keywords="limit cycle, stability, Floquet, multiple shooting",
    url='https://github.com/vortexlab-uclouvain/multiflap',
    packages=find_packages(exclude=['img', 'amin', 'pdf']),
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    install_requires=[
        'numpy>=1.14,<2.0',
        'matplotlib',
        'scipy'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: Apache License, Version 2.0, January 2004',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering'
    ]
)
