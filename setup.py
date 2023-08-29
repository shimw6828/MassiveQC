from setuptools import setup, find_packages

# with open('requirements.txt') as f:
#     requirements = f.readlines()
requirements = [
    'scikit-learn',
    'shap',
    'xopen',
    'NumPy',
    'Pandas >=1.3.2',
    'fastparquet',
    'more-itertools',
    'tqdm'
]

with open('README.md') as f:
    long_description = f.read()


setup(
    name='MassiveQC',
    version='0.1.2',
    author='shimw6828',
    author_email='shimw6828@qq.com',
    url='https://github.com/shimw6828/MassiveQC',
    description='Tools for QC massive RNA-seq samples',
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='MIT',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'MultiQC = MassiveQC.MultiProcess:main',
            'SingleQC = MassiveQC.SingleProcess:main',
            'IsoDetect = MassiveQC.IsoDetect:main'
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX",
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    install_requires=requirements,
    python_requires='>=3.7'
)