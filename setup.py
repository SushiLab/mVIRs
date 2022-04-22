from setuptools import setup
import os

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()
install_requires = ['pysam']
long_description = read('README.md')
setup(
    name = "mVIRs",
    version = "1.1.1",
    author = "Hans-Joachim Ruscheweyh",
    author_email = "hansr@ethz.ch",
    description = ("mVIRs: Localisation of inducible prophages using NGS data"),
    license = "GPLv3",
    include_package_data=True,
    install_requires=install_requires,
    long_description=long_description,
    long_description_content_type='text/markdown',
    keywords = "bioinformatics metagenomics NGS alignment OPRs Prophages",
    url = "https://github.com/SushiLab/mVIRs",
    packages=['mVIRs'],
    download_url = "https://github.com/SushiLab/mVIRs/archive/refs/tags/1.1.1.tar.gz",
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    entry_points = {
        'console_scripts': ['mvirs=mVIRs.mvirs:main'],
    }
)
