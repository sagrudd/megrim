import setuptools
import megrim

with open("README.md", "r") as fh:
    long_description = fh.read()
setuptools.setup(
    name='megrim',
    version = megrim.__version__,
    author="Stephen Rudd",
    author_email="stephen.rudd@nanoporetech.com",
    description="A bioinformatics tutorial framework for epi2me-labs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sagrudd/megrim",
    packages=setuptools.find_packages(),
    package_data = {
        'megrim': ['data/ONT_logo.png', 'data/FontAwesome5Free-Solid-900.otf']
    },
    include_package_data=True,
    install_requires=['datetime', 'pytz', 'progressbar', 'tqdm', 'ipython', 
                      'ipywidgets', 'matplotlib', 'pandas', 'numpy',
                      'pyranges', 'pysam', 'fontawesome', 'scipy', 
                      'dask[complete]', 'bamread', 'weightedcalcs',
                      'biopython', 'ont_fast5_api'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)",
        "Operating System :: OS Independent",
    ],
)
