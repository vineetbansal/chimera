from setuptools import setup, find_namespace_packages

setup(
    name='chimera',
    version='0.1.0',
    description='Nothing useful in particular',
    url='https://github.com/vineetbansal/chimera',
    author='Vineet Bansal',
    author_email='vineetb@princeton.edu',
    license='MIT',

    package_dir={'': 'src'},
    packages=find_namespace_packages(where='src'),
    package_data={'chimera': ['config.json'], 'chimera.data': ['*.tsv'], 'chimera.data.pfms': ['*.pfm']},

    zip_safe=True,
    test_suite='tests'
)
