from setuptools import setup, find_packages

setup(
    name="chronocci",
    version="0.1.1",
    author="Your Name",
    author_email="your.email@example.com",
    description="Chronological Inference of Multi-Lineage Cell-Cell Interactions Using CellRank and LIANA",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.9, <3.12", 
    install_requires=[
        "scanpy==1.10.2",            
        "cellrank==2.0.4",           
        "liana>=1.0.0,<=1.2.0",      
        "numpy>=1.22.0,<2.0.0",      
        "pandas>=1.4.0,<2.2.0",      
        "matplotlib>=3.5.0,<3.9.0",
    ],

)
