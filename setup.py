import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
        name="pyconsFold",
        version="0.1",
        author="John Lamb",
        author_email="john@biolamb.se",
        description="A python wrapper around CNS for modelling and docking, inspired by CONFOLD",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/johnlamb/pyconsFold",
        packages=setuptools.find_packages(),
        include_package_data=True,
        zip_safe=False,
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: GPLv3",
            "Operating System :: OS Independent",
            ],
        python_requires='>3.5',
        install_requires=[
            'biopython>=1.7',
            'numpy>=1.19',
            'scipy>=1.4.0'
            ],
        )
