from setuptools import setup, find_namespace_packages

setup(
    name='chimera',
    version='0.1.4',
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
