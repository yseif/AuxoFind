import os
from setuptools import setup, Command, find_packages


class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')
setup(
    name='AuxoFind',
    version='0.0.1',
    description='A package for predicting auxotrophies using metabolic models.',
    license='MIT',
    packages=find_packages(exclude=['tests*']),
    long_description=open('README.md').read(),
    install_requires=["cobrapy" ],
    url='https://github.com/yseif/panMEM.git',
    author='Yara Seif',
    cmdclass={
        'clean': CleanCommand,
    }
)
