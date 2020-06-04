from setuptools import setup


setup(
    name = "OPR Finder",
    version = "0.1",
    author = "Hans-Joachim Ruscheweyh",
    author_email = "hansr@ethz.ch",
    description = ("Bioinformatic toolkit for finding OPRs in sequencing data"),
    license = "GPLv3",
    keywords = "bioinformatics metagenomics ngs OPR ",
    url = "https://github.com/SushiLab/OPR_Finder",
    packages=['oprfinder'],
    entry_points = {
        'console_scripts': ['oprfinder=oprfinder.oprfinder:main'],
    }
)