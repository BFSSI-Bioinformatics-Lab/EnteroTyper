import setuptools
from EnteroTyper.__init__ import __version__, __author__, __email__

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="EnteroTyper",
    install_requires=['click', 'pandas', 'xlsxwriter', 'tqdm', 'dataclasses', 'pytest'],
    tests_require=['pytest'],
    setup_requires=['pytest-runner'],
    python_requires='~=3.6',
    description="Scripts to type and compare assemblies based off of Enterobase schemes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BFSSI-Bioinformatics-Lab/EnteroTyper",
    package_data={'EnteroTyper': ['*']},
    install_package_data=True,
    packages=setuptools.find_packages(),
    version=__version__,
    author=__author__,
    author_email=__email__,
    entry_points={
        'console_scripts': [
            'enterotyper=EnteroTyper.enterotyper:enterotyper'
        ]
    }
)
