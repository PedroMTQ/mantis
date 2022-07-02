import pathlib
from setuptools import setup,find_packages

package_name='Mantis'
def get_package_version():
    import os
    from sys import platform
    if platform.startswith('win'):
        SPLITTER = '\\'
    else:
        SPLITTER = '/'


    dir_name=os.path.dirname(os.path.abspath(__file__))
    init_path=f'{dir_name}{SPLITTER}{package_name.lower()}{SPLITTER}__init__.py'
    package_version=None
    with open(init_path) as file:
        for line in file:
            if '__version__' in line:
                package_version=line.replace('__version__','')
                package_version=package_version.strip('\n')
                package_version=package_version.strip()
                package_version=package_version.strip('=')
                package_version=package_version.strip()
                package_version=package_version.strip('"')
                package_version=package_version.strip('"')
    return package_version



# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text(encoding='utf-8')
LICENSE = (HERE / "LICENSE").read_text(encoding='utf-8')

long_description='Mantis is protein function annotation, that dynamically integrates multiple reference databases to produce consensus-driven annotations.'

setup(
    name=package_name,
    version=get_package_version(),
    author="Pedro Queirós",
    author_email="pdqueiros@gmail.com",
    description="Tool for protein function annotation.",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/PedroMTQ/Mantis",
    project_urls={
        "Bug Tracker": "https://github.com/PedroMTQ/Mantis/issues",
    },
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        'Operating System :: MacOS',
        'Operating System :: POSIX :: Linux',
    ],
    license=LICENSE,
    include_package_data=True,
    install_requires=['cython', 'nltk>=1.18.1', 'psutil>=5.6.7', 'requests>=2.22.0', 'diamond>=2.0.13', 'python>=3.7'],
    entry_points={
        "console_scripts": [
            "mantis=mantis.__main__:main",
        ],
    },
)
