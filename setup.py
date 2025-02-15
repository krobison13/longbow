from io import open

from setuptools import find_packages, setup

# following src dir layout according to
# https://blog.ionelmc.ro/2014/05/25/python-packaging/#the-structure
version = "0.5.46"
setup(
    name="maslongbow",
    version=version,
    description="Annotation and segmentation of MAS-seq data",
    url="https://broadinstitute.github.io/longbow/",
    author="Kiran V Garimella, Jonn Smith",
    author_email="kiran@broadinstitute.org, jonn@broadinstitute.org",
    license="BSD 3-Clause",
    long_description=open("README.rst").read(),
    install_requires=[
        "cython",
        "click",
        "click-log",
        "editdistance",
        "pomegranate @ git+https://github.com/jonn-smith/pomegranate.git@d80ec9304bba4c42ee3943c887754fd78885cd55",
        "pandas",
        "numpy",
        "matplotlib",
        "ordered-set",
        "pysam==0.19.1",
        "progress",
        "construct",
        "pathos",
        "ssw==0.4.1",
        "tqdm",
        "symspellpy"
    ],
    tests_require=["coverage", "pytest"],
    python_requires=">=3.6, <3.8",
    packages=find_packages("src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: Implementation :: CPython",
    ],
    entry_points={"console_scripts": ["longbow=longbow.__main__:main_entry"]},
    include_package_data=True,
)
