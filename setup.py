from setuptools import setup, find_packages

setup(
    name='IMMerge',
    version='0.0.3',
    license='MIT',
    author="Wanying Zhu from Below lab",
    author_email='email@example.com',
    long_description='Merge large Variant Call Format (VCF) files',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    url='https://github.com/belowlab/IMMerge',
    keywords='Genetics imputation VCF_merger',
    install_requires=[
          'pandas', 'xopen'
      ],
    python_requires='>=3',

)
