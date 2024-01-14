import setuptools

with open("README.md", "r") as fh:
    description = fh.read()

setuptools.setup(
    name="hwetests",
    version="0.5.8",
    author="Or Shkuri",
    author_email="orshkuri2000@gmail.com",
    packages=["hwetests"],
    description="A python package containing two statistical tests for HWE testing: Gibbs Sampling test and a modified Chi Squared test that handles ambiguity",
    long_description=description,
    long_description_content_type="text/markdown",
    url="https://github.com/ExtraFlash/HWE_tests_package",
    license='MIT',
    python_requires='>=3.8',
    install_requires=[]
)