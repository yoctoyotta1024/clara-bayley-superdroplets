from setuptools import setup, find_packages

setup(
    name='version2.0',
    version='2.0',
    packages=find_packages(),
    install_requires=[
        'pytest',
        'sphinx',
    ],
)
