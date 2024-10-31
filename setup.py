from setuptools import setup, find_packages

setup(
    name='yatomat',
    version='0.1.2',
    author='Marina Popova, Aleksey Komissarov',
    author_email='ad3002@example.com',
    description='Yet Another Tool for Making Artificial Genomes',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/aglabx/yatomat',
    # Изменяем packages чтобы правильно находить модули
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3.6',
    install_requires=[
        'pandas',
        'numpy',
        'plotly',
        'matplotlib',
        'scipy',
    ],
    include_package_data=True,
    # Исправляем entry_points чтобы правильно указать путь к main функции
    entry_points={
        'console_scripts': [
            'yatomat=yatomat.yatomat:main',
        ],
    },
    project_urls={
        'Bug Reports': 'https://github.com/aglabx/yatomat/issues',
        'Source': 'https://github.com/aglabx/yatomat',
    },
)
