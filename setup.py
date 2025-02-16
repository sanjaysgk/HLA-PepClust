from setuptools import setup, find_packages

setup(
    name="HLA-PEPCLUST",
    version="0.1.0-dev",
    author='Sanjay Krishna', 
    author_email='sanjay.sondekoppagopalakrishna@mail.com',
    packages=find_packages(),
    long_description=open('README.md').read(),
    install_requires=[
        "rich==13.9.4",
        "click==8.1.8",
        "pandas==2.2.3",
        "seaborn==0.13.2",
        "matplotlib==3.9.4",
        "rich-argparse==1.7.0",
    ],
    entry_points={
        "console_scripts": [
            "clust-search=cli.main:main",
        ],
    },
)
