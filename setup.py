from setuptools import setup

setup(
   name='ship_plasmid',
   version='0.2.2',
   description='A a tool to find antimicrobial resistant (AMR) DNA regions with evidence of lateral transfer between plasmids.',
   author='Marco Teixeira',
   author_email='mcarvalh@broadinstitute.org',
   license='GNU General Public License v3 (GPLv3)',
   packages=['ship_plasmid'],
   install_requires=[
        "python==3.8",
        "bcbio-gff==0.6.9",
        "biopython==1.79",
        "joblib==1.2.0",
        "markov_clustering==0.0.6.dev0",
        "matplotlib==3.6.1",
        "networkx==2.8.6",
        "numpy==1.23.3",
        "pandas==1.4.4",
        "pymc==4.2.2",
        "pyvis==0.2.1",
        "pyyaml==6.0",
        "scikit-learn==1.1.2",
        "scipy==1.9.1",
        "tqdm==4.64.1"
    ],
   scripts=[
        'ship_plasmid/ship.py'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'GNU General Public License v3 (GPLv3)',  
        'Programming Language :: Python :: 3.8',
    ]
)