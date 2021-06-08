from setuptools import setup


setup(
    name = "mVIRs",
    version = "1.1.0",
    author = "Hans-Joachim Ruscheweyh",
    author_email = "hansr@ethz.ch",
    description = ("Bioinformatic toolkit for finding prophages in sequencing data"),
    license = "GPLv3",
    keywords = "bioinformatics metagenomics NGS alignment OPRs Prophages",
    url = "https://github.com/SushiLab/mVIRs",
    packages=['mVIRs'],
    entry_points = {
        'console_scripts': ['mvirs=mVIRs.mvirs:main'],
    }
)
