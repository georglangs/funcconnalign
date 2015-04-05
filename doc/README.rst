Documentation source is written using restructuredText and can be built into
html and latex/pdf formats using the Python tool Sphinx.

I found the following documents useful as intros to restructuredText

- http://sphinx-doc.org/rest.html
- https://pythonhosted.org/an_example_pypi_project/sphinx.html

The full guide can be found here: http://docutils.sourceforge.net/docs/ref/rst/restructuredtext

In order to write documentation, make sure you have the following installed:

1. Python 2.7
2. pip 1.3
3. sphinx 1.2
4. sphinxcontrib-matlabdomain 0.2.6
5. sphinxcontrib-bibtex
6. readthedocs theme for sphinx
7. LaTeX (only necessary for pdf generation)

I believe LaTeX is optional for html generation, but is required for pdf
generation. For OS X, I recommend MacTeX https://www.tug.org/mactex/.
For Windows, I've used MiKTeX in the past: http://miktex.org/about.

I like using conda to manage multiple Python environments. Anaconda is
an package that includes conda, Sphinx and lots of NumPu/SciPy related
extensions. You can download Anaconda here: https://store.continuum.io/cshop/anaconda/.

Once you have Python installed, it's easy to install Sphinx, its
matlab-domain extension and the readthedocs theme using pip::

    pip install sphinx
    pip install sphinxcontrib-matlabdomain
    pip install sphinx_rtd_theme
    pip install sphinxcontrib-bibtex

To generate the html files, go to doc/external or doc/internal and run::

    make html

and to generate the pdf files, run::

    make latexpdf
