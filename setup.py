from setuptools import setup, find_packages

setup(
    name="research_scripts",  # or whatever name you want to use
    version="0.1",
    packages=find_packages(),
    install_requires=[
        'prody',  # list your dependencies here
        'numpy',
        'griddataformats',
        'mdanalysis'
    ],
)