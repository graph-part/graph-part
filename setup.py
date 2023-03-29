import os
from setuptools import setup, find_packages


with open(os.path.join("README.md")) as f:
    readme = f.read()

requirements = [
    "numpy>=1.19.5",
    "networkx>=2.5.1",
    "tqdm>=4.62.3",
    "pandas>=1.1.5"
]




setup(
    name="graph-part",
    version='0.1.2',
    description="Graph-based partitioning of biological sequence data",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://healthtech.dtu.dk",
    author="F. Teufel and M.H. Gislason",
    packages=['graph_part'],
    python_requires=">=3.6, <4",
    install_requires=requirements,
    scripts=['graph_part/mmseqs_fake_prefilter.sh'],
    entry_points = {"console_scripts":['graphpart=graph_part:run_graph_part']},
)
