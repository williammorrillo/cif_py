import setuptools

__version__ = '0.0.1'

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
        name="cif_py",
        version=__version__,
        author="William Morrillo",
        author_email="williammorrillo@gmail.com"
        description="A python package for reading and writing CIF files",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url={
            "Source Code": "https://github.com/williammorrillo/cif_py",
            }
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Operating System :: OS Independent",
            ],
        packages=setuptools.find_packages(),
        python_requires='>=3.6',
        install_requires=[
            'numpy',
            'h5py',
            ],
        entry_points={
            'console_scripts': [
                'cif_py=cif_py.cli:main',
                ],
            },
        )
