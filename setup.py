try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'PoWeRS: Python WRapper for MOOG Synth',
    'author': 'Leonardo dos Santos',
    'download_url': 'https://github.com/RogueAstro/PoWeRS',
    'author_email': 'leonardoags@usp.br',
    'version': '0.1.160118',
    'install_requires': ['numpy','matplotlib'],
    'packages': ['powers'],
    'name': 'powers'
}

setup(**config)
