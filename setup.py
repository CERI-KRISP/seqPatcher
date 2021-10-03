from sys import version
from setuptools import setup

setup(
    name="SeqPatcher",
    version="0.0.1",
    scripts=["seqpatcher.py"],
    install_requires=[
        "pandas==1.2.4",
        "biopython==1.78",
        "click==7.1.2",
    ],
    # packages=[],
    # package_dir={},
    url="https://github.com/krisp-kwazulu-natal/seqPatcher",
    license="GPLv3",
    author="Anmol Kiran; San James Emmanuel",
    author_email="anmol.kiran@gmail.com;sanemmanueljames@gmail.com",
    description="Inserts sars-cov2 S-gene sanger sequencing data in assemblies generated using HTS.",
)
