import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
setuptools.setup(
    name='megrim',
    version='0.1.3.9',
    author="Stephen Rudd",
    author_email="stephen.rudd@nanoporetech.com",
    description="A bioinformatics tutorial framework for epi2me-labs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sagrudd/megrim",
    packages=setuptools.find_packages(),
    install_requires=['datetime', 'pytz', 'progressbar', 'tqdm', 'ipython', 'ipywidgets'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)",
        "Operating System :: OS Independent",
    ],
)
