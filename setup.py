from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="multiplex",
    version="2020.1.0",
    author="Benedikt Groever",
    author_email="groever@seas.harvard.edu",
    description="A package for personal use",
    long_description=long_description,
    include_package_data=True,
    #packages=['multiplex'],#find_packages('pyautodiff'),
    packages=find_packages('multiplex'),
    #install_requires=['fenics','mshr','pyyaml','matplotlib'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ]
)
