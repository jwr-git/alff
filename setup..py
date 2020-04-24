import setuptools

with open("README.md", "r") as f:
    long_description = f.read()
    
setuptools.setup(
    name="jwr-git",
    version="1.1",
    author="Jamie Robinson",
    author_email="jwr_git@outlook.com",
    description="Package to find allele frequencies using NCBI's variation service",
    long_description =long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jwr-git/alff",
    packages=setuptools.find_packages(),
    classifiers=[
        "License :: GNU v3 License",
        "Operating System :: OS Independent",
        "Programming Language :: Python"
        ],
    python_requires='>=3'
    )
