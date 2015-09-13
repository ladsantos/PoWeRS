try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'pwoogs: python wrapper for MOOG synth',
    'author': 'Leonardo dos Santos',
    'download_url': 'https://github.com/laugustogs/pwoogs',
    'author_email': 'leonardoags@usp.br',
    'version': '0.1.150912',
    'install_requires': ['numpy','matplotlib'],
    'packages': ['pwoogs'],
    'name': 'pwoogs'
}

setup(**config)