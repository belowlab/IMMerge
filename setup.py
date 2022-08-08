from setuptools import setup, find_packages

setup(
    name='IMMerge',
    version='0.0.1',
    license='MIT',
    author="Wanying Zhu from Below lab",
    author_email='email@example.com',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    url='https://github.com/belowlab/IMMerge',
    keywords='Genetics imputation VCF_merger',
    install_requires=[
          'pandas', 'xopen'
      ],
    python_requires='>=3',

)
