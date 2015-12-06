import os
from setuptools import setup

setup(
    name = "Xylem",
    version = "0.0.1",
    author = "Pranjal Vachaspati",
    author_email = "pr@nj.al",
    description = ("phylogenetic pipelines"),
    license = "GPLv3",
    keywords = "phylogenetics pipeline",
    url = "",
    package_dir={'xylem':'src'},
    packages=['xylem', 'xylem.Tasks']
)
