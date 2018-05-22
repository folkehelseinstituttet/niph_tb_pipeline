from setuptools import setup
from distutils.util import convert_path

versionholder = {}
ver_path = convert_path('niph_tb_pipeline/__init__.py')
with open(ver_path) as init_file:
    exec(init_file.read(), versionholder)

setup(name='niph_tb_pipeline',
    version=versionholder['__version__'],
    description='Complete analysis pipeline for MTBC WGS',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.'
    ],
    keywords=['TB','lineage','classification','typing','mycobacterium','tuberculosis', 'AMR', 'reference'],
    author='Ola Brynildsrud',
    author_email='ola.brynildsrud@fhi.no',
    license='MIT',
    packages=['niph_tb_pipeline'],
    install_requires=[
        'biopython',
        'ete3'],
    entry_points={
        'console_scripts': ['niph_tb_pipeline=niph_tb_pipeline.niph_tb_pipeline:main',
                            'json_to_tsv=niph_tb_pipeline.json_to_tsv:main',
                            'tex_finalizer=niph_tb_pipeline.Tex_finalizer:CreateReport']
    },
    package_data={'niph_tb_pipeline': ['LICENSE', 'setup.py','README.md','docs/*','data/*']},
    include_package_data=True
    )
