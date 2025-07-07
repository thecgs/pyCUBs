from setuptools import setup, find_packages

setup(name="pycubs",
      version="2.0",
      description='An integrated package for codon usage bias analysis',
      url='https://github.com/thecgs/pyCUBs/',
      author='Guisen Chen',
      author_email='thecgs001@foxmail.com',
      long_description="It provides not only the calculation of several popular codon bias indices, but also several conventional codon usage bias analysis methods, such as neutrality plot analysis, ENC-GC3s plot analysis, parity rule 2-bias plot analysis, optimal codon analysis, correspondence analysis, principal component analysis, and phylogenetic analysis. It provides an interactive analysis method for freely exploring codon bias and visualizing these results.",
      license='MIT License',
      packages=find_packages(exclude=["example", ".github"]),
      keywords=['bioinformatics', 'codon usage bias', 'codonW', 'cusp'],
      install_requires=["scipy>=1.10.1",
                        "scikit-bio>=0.6.3",
                        "numpy>=1.24.4",
                        "pandas>=2.0.3",
                        "prince>=0.13.0",
                        "seaborn>=0.13.1",
                        "matplotlib>=3.7.5"]
     )
