import os
from setuptools import setup, find_packages

PATH_ROOT = os.path.dirname(__file__)


with open(os.path.join(PATH_ROOT, "README.md")) as f:
    readme = f.read()

with open(os.path.join(PATH_ROOT, "requirements.txt")) as f:
    requirements = f.read()




setup(
    name="graph-part",
    version='1.0',
    description="Graph-based partitioning of biological sequence data",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="healthtech.dtu.dk",
    author="Magnús Halldór Gíslason",
    author_email="mhgislason@gmail.com",
    packages=['graph_part'],
    python_requires=">=3.6, <4",
    install_requires=requirements,
    entry_points = {"console_scripts":['graphpart=graph_part:run_graph_part']},
)
