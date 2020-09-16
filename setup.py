from setuptools import setup

setup(
    name='Mantis',
    version='1.0.0',
    author='Pedro Queirós',
    author_email='pdqueiros@gmail.com',
    packages=['mantis'],
    scripts=['bin/stowe-towels.py','bin/wash-towels.py'],
    url='https://github.com/PedroMTQ/mantis',
    license='LICENSE',
    description='A package to annotate protein sequences.',
    install_requires=[
        'hmmer == 3.2.1',
        'nltk == 3.4.4',
        'numpy == 1.18.1',
        'pip',
        'psutil == 5.6.7',
        'python == 3.7.6',
        'requests == 2.22.0',
        'setuptools == 45.1.0',
    ],
)
