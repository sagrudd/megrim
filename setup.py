import setuptools
import megrim

with open("README.md", "r") as fh:
    long_description = fh.read()
setuptools.setup(
    name='megrim',
    version = megrim.__version__,
    author="Stephen Rudd",
    author_email="stephen@mnemosyne.co.uk",
    description="Python implementation of the floundeR workflow",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sagrudd/megrim",
    packages=["megrim", "megrim.plugins"],
    package_data = {
        'megrim': ['data/ONT_logo.png', 'data/FontAwesome5Free-Solid-900.otf']
    },
    include_package_data=True,
    entry_points={'console_scripts': ["megrim=megrim.toolbox:main"]},
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
