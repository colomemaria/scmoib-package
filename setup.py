import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scmoib",
    version="0.2.1",
    author="Atai Dobrynin",
    author_email="atay.dobrynin@gmail.com",
    description="SCMOIB package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/colomemaria/scmoib-package",
    project_urls={
        "Bug Tracker": "https://github.com/colomemaria/scmoib-package/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    python_requires=">=3.6",
)
