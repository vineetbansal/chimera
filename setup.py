import os
from setuptools import setup, find_namespace_packages


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def get_version():
    """
    Get package version (without import the package, which may or may not work)
    :return: version info in maj.min.bld.0 format
    """
    version_dict = {}
    exec(open("src/chimera/version.py").read(), version_dict)
    return version_dict['version']


setup(
    name='chimera',
    version=get_version(),
    description='Protein Domain Identification Package',
    url='https://github.com/vineetbansal/chimera',
    author='Vineet Bansal',
    author_email='vineetb@princeton.edu',
    license='MIT',

    package_dir={'': 'src'},
    packages=find_namespace_packages(where='src'),
    include_package_data=True,

    zip_safe=False,
    test_suite='tests'
)
