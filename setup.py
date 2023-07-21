from setuptools import setup

setup(
    name='spacells',
    version='0.0.1',
    packages=['spacells'],
    install_requires=[
        'requests',
        'numpy',
        'scipy',
        'scikit-learn',
        'importlib-metadata; python_version == "3.8"',
    ],
)
