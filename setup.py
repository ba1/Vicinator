import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Vicinator", # Replace with your own username
    version="0.0.9",
    author="Ba1",
    author_email="djahanschiri@bio.uni-frankfurt.de",
    description="A small python package to trace orthology neighborhood across feature files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ba1/vicinator",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        "ansi2html>=1.5.2",
        "colorama>=0.4.4",
        "ete3>=3.1.2",
        "pandas>=1.1.3"
    ],
    entry_points={
        'console_scripts': [
            'vicinator=vicinator.vicinator:main'
        ]
    },
)
