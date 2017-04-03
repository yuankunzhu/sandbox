from setuptools import setup

setup(
    name='pgcloader',
    version='1.0.0.a1',
    packages=['pgcloader'],
    description='A sample Python API wrapper for load data from S3 to Cavatica',
    url='https://github.com/d3b-center/s3-sbg-loader',
    author='Yuankun Zhu',
    author_email='zhuyuank@gmail.com',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7'
    ],
    install_requires=[
        'requests',
        'ConfigParser',
        'json', 'csv'
    ],
)
