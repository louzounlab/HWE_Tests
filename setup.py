import setuptools

with open("README.md", "r") as fh:
    description = fh.read()

setuptools.setup(
    name="hwetests",
    version="0.6.9",
    author="Or Shkuri",
    author_email="orshkuri2000@gmail.com",
    packages=["src/hwetests"],
    description="A python package containing two statistical tests for HWE testing: Gibbs Sampling test and a modified Chi Squared test that handles ambiguity",
    long_description=description,
    long_description_content_type="text/markdown",
    url="https://github.com/ExtraFlash/HWE_tests_package",
    license='MIT',
    python_requires='>=3.8',
    install_requires=['numpy>=1.22.4',
                      'pandas>=1.4.3',
                      'scipy>=1.7.3',
                      'matplotlib>=3.7.4',
                      'seaborn>=0.12.0']
)
