from pathlib import Path
from setuptools import setup, find_packages

author = "Quan Xu"
author_email = "qxuchn@gmail.com"
description = "A method to map new scRNA-seq data to HEOCA"

long_description = Path("README.md").read_text("utf-8")
requirements = [
    l.strip() for l in Path("requirements.txt").read_text("utf-8").splitlines()
]

setup(
    name="sc2heoca",
    version="0.3.0",
    author=author,
    author_email=author_email,
    description=description,
    license="MIT",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    include_package_data=True,
    install_requires=requirements,
    python_requires="==3.9.16",
)

