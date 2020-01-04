from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="multiplex",
    version="0.0.1",
    author="Benedikt Groever",
    author_email="groever@seas.harvard.edu",
    description="A package for personal use",
    long_description=long_description,
    packages=['multiplex'],#find_packages('pyautodiff'),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)