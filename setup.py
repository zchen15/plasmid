from setuptools import setup, find_packages

def get_version():
    fname = 'src/plasmid/__init__.py'
    out = {'__version__':'',
           '__author__':'',
           '__name__':'',
           '__email__':''}
    with open(fname,'r') as f:
        x = f.read()
    x = x.split('\n')
    for i in x:
        for k in out.keys():
            if k in i:
                val = i.split(' = ')[-1]
                out[k] = val.replace('"','')
    return out

info = get_version()
__version__ = info['__version__']
__author__ = info['__author__']
__name__ = info['__name__']
__email__ = info['__email__']

def get_requirements(fname):
    with open(fname,'r') as f:
        out = f.read()
    out = out.split('\n')
    return out

def read(path):
    with open(path, 'r') as f:
        out = f.read()
    return out

long_description = read('README.md')

setup(
    install_requires=get_requirements('requirements.txt'),
    python_requires='>=3.6',
    name=__name__,
    version=__version__,
    license='GPL-3',
    author=__author__,
    author_email=__email__,
    maintainer=__author__,
    maintainer_email=__email__,
    description='Plasmid: A python package for the visualizing, editing, and assembling DNA',
    long_description=long_description,
    long_description_content_type='text/markdown',
    zip_safe=False,
    packages=find_packages(include=[__name__]),
    include_package_data=False,
    url='https://github.com/zchen15/plasmid',
    keywords='gene editor, plasmid, genbank, dna assembly, genomics, synthetic biology',
    classifiers=[
        'Environment :: Console',
        'Framework :: Jupyter',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GPL-3 License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Utilities'
    ]
)

