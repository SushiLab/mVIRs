import os

from setuptools import setup, Extension
# from setuptools.command.build_ext import build_ext as _build_ext
from Cython.Build import build_ext

try:
    from Cython.Build import cythonize
except ImportError as err:
    cythonize = err

SRC_DIR = "mVIRs"
PACKAGES = [SRC_DIR]

EXTENSIONS = [
    Extension("mVIRs.oprs_c", [SRC_DIR + "/oprs_c.pyx"], language="c++"),
    Extension("mVIRs.extract_regions", [SRC_DIR + "/extract_regions.pyx"]),
]


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()
install_requires = ['pysam']
long_description = read('README.md')
setup(
    ext_modules=EXTENSIONS,
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
    },
    cmdclass={'build_ext': build_ext}
)
