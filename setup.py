import os
from setuptools import setup, find_packages


with open(os.path.join("README.md")) as f:
    readme = f.read()

with open(os.path.join("requirements.txt")) as f:
    requirements = f.read()




setup(
    name="graph-part",
    version='1.0',
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
