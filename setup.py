from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="chemreact",
    version="0.0.3",
    author="Ilia Dorokhov",
    author_email="i.dorokhov14@ic.ac.uk",
    description="Reaction simulation package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dorox/react",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Operating System :: OS Independent",
    ],
    python_requires = ">=3.6",
    install_requires = ["numpy", "scipy", "matplotlib"]
)