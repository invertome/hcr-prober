from setuptools import setup, find_packages

with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='hcr-prober',
    version='1.9.6',
    packages=find_packages(),
    author='Jorge L. Perez-Moreno, Ph.D.',
    description='An advanced tool for HCR v3.0 probe design with a high-performance, iterative filtering pipeline.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=[
        'biopython>=1.79', 'pandas>=1.3.0', 'openpyxl>=3.0.0',
        'numpy>=1.20.0', 'matplotlib>=3.3.0', 'PyYAML>=5.4.0', 'loguru>=0.5.3'
    ],
    entry_points={'console_scripts': ['hcr-prober = hcr_prober.main:main']},
    package_data={'hcr_prober': ['config/amplifiers/*.json', 'config/hcr-prober.yaml']},
    include_package_data=True,
    python_requires='>=3.8',
)
