from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.readlines()

setup(
    name='MassiveQC',
    version='1.0.0',
    author='shimw6828',
    author_email='shimw6828@qq.com',
    #url='https://github.com/shimw6828/MassiveQC',
    description='Tools for QC massive RNA-seq samples',
    license='MIT',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'MultiQC = MassiveQC.MultiProcess:main',
            'SingleQC = MassiveQC.SingleProcess:main',
            'IsoDetect = MassiveQC.IsoDetect:main'
        ]
    },
    classifiers=(
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX",
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ),
    install_requires=requirements
)