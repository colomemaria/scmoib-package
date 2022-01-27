import sys, os
if sys.version_info < (3,):
    sys.exit('scmoib requires Python >= 3.6')

from setuptools import setup, find_packages
from pathlib import Path
from versioneer import get_version, get_cmdclass
try:
    from scmoib import __author__, __email__
except ImportError:  # Deps not yet installed
    __author__ = __email__ = ''

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="scmoib",
    version=get_version(),
    cmdclass=get_cmdclass(),
    author=__author__,
    author_email=__email__,
    license='MIT',
    description="Single Cell Multi Omic Integration Benchmarking",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/colomemaria/scmoib-package",
    project_urls={
        "Bug Tracker": "https://github.com/colomemaria/scmoib-package/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    python_requires=">=3.6",
    install_requires=[
        l.strip() for l in
        Path('requirements.txt').read_text('utf-8').splitlines()
    ],
)
