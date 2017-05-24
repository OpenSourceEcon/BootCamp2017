# JRC Overlapping Generations Model Training by Open Research Group

This public repository contains the training materials, tutorials, and code for the two week training delivered to the Joint Research Centre of the European Commission in Seville, Spain, May 29 to June 9, 2017 by [Richard Evans](https://sites.google.com/site/rickecon/) and [Jason DeBacker](http://jasondebacker.com/) of [Open Research Group, Inc.](http://openrg.com/) (OpenRG).

We are planning on the following schedule each of the 10 workdays during that time. Let us know if we need to change this schedule in any way.

* 9am to noon: Lecture, theory, computational instruction
* noon to 1pm: Lunch
* 1pm to 4pm: Guided computational practice, implementation, problem sets

We have created a series of textbook chapters complete with theoretical exposition, computational description and tips, and problem sets. These chapters and their instruction are meant to progressively build upon each other. In the end, you will be building a computational implementation of an OG model for fiscal policy that is thousands of lines of code. We will train your research group to understand each code section and to write the code in a way that is accessible, modular, scalable, and amenable to collaboration among a large group of researchers.

This will be an intensive two weeks. We are providing your researchers 7 areas of tutorials that they will benefit from reading before the training. We will, of course, teach these things as we go through the material. But we will be able to proceed at a faster pace if the attendees are already familiar with most of the concepts below.

Pre-course Tutorial Areas

1. Instructions for installing the Anaconda distribution of Python
2. Text editor suggestions (Atom, Sublime Text 3, Vim)
3. PEP8, docstring commenting, and module structure
4. Git and GitHub tutorial
5. Jupyter Notebooks
6. Basic Python tutorials (data structures, logic, functions and modules, pandas, root finders and minimizers)


## 1. Instructions for installing the Anaconda distribution of Python

We will be using the [Python](https://www.python.org/) programming language and many of its powerful libraries for writing the code that will run the overlapping generations models in this training. Using an open source language, such as Python, has the advantage of being free and accessible for anyone who wishes to contribute to this project. Being open source also allows Python users to go into the source code of any function to modify it to suit one's needs.

We recommend that each participant download the Anaconda distribution of Python provided by [Continuum Analytics](https://www.continuum.io/). We recommend the most recent stable version of Python, which is currently Python 3.6. This can be done from the [Anaconda download page](https://www.continuum.io/downloads) for Windows, Mac OSX, and Linux machines. The code we will be writing uses common Python libraries such as `NumPy`, `SciPy`, `pickle`, `os`, `matplotlib`, and `time`.


## 2. Text editor suggestions

In our recommended Python development workflow, you will write Python scripts and modules (`*.py` files) in a text editor. Then you will run those scripts from your terminal. You will want a capable text editor for developing your code. Many capable text editors exist, but we recommend three.

1. [Atom](https://atom.io)
2. [Sublime Text 3](https://www.sublimetext.com)
3. [Vim](http://www.vim.org)

Atom and Vim are completely free. A trial version of Sublime Text 3 is available for free, but a licensed version is $70 (US dollars). In the following subsections, we give some of the details of each of the above three text editors.


### 2.1. Atom

[Atom](https://atom.io) is an open source text editor developed by people at GitHub.com. This editor has all the features of Sublime Text 3, but it also allows users full customizability. Further, it has been a while now that the users of Atom have surpassed the critical mass necessary to keep the editor progressing with the most cutting edge additions.

There are several packages you'll want to install with Atom.  Once Atom is installed, you can add packages by navigating Atom->Preferences->Install and then typing in the name of the package you would like to install.  

For work with Python, we recommend the following packages be installed:

* MagicPython
* python-indent
* tabs-to-spaces
* minimap
* open-recent
* linter-python-pep8

For development with GitHub we recommend:

* merge-conflict

If using LaTex in this editor, the following packages are helpful:

* atom-latex
* latextools
* autocomplete-bitex
* dictionary
* latexer
* pdf-view

In addition, you will also want to download the [Skim](http://skim-app.sourceforge.net) PDF viewer to aid in displaying PDF files compiled from TeX with Atom.


### 2.2 Sublime Text 3

[Sublime Text 3](https://www.sublimetext.com) is the most widely used and versatile private software text editor. It has tremendous flexibility, as well as the polish of a piece of professional software. Sublime Text 3 will cost $70 for a license, although you can use a trial version indefinitely without charge while only having to suffer through frequent reminders to buy the full version.


### 2.3 Vim

[Vim](http://www.vim.org) is free and very powerful. Vim is the hard core developer's text editor of choice. The learning curve for using vim is a little steeper than that of Atom and Sublime Text 3, but it also has some advantages for efficient programming. Vim has navigation that does not use a mouse or trackpad. Eventually, your fingers never leave your keyboard. Further, most terminals have Vim built in so you can use Vim to edit scripts and modules on the fly with your terminal session.  


## 3. PEP 8, docstring commenting, and module structure

Computer code executes some set of commands in an organized way. In every case, there are often many ways to execute a set of instructions--some ways more efficient than others. However, code has at least three functions.

1. Efficiently execute the task at hand.
2. Be accessible and usable to other programmers.
3. Be scalable and integrable with other projects and procedures.

Bill Gates is credited with the following plea for efficiency and parsimony in code writing.

> "Measuring programming progress by lines of code is like measuring aircraft building progress by weight."

Strong support for points (2) and (3) is Eagleson's Law.

> "Any code of your own that you haven't looked at for six or more months might as well have been written by someone else."

Because of the latter two characteristics, Python code has developed some conventions and best practices, some of which have been institutionalized in the [PEP 8--Style Guide for Python Code](https://www.python.org/dev/peps/pep-0008/) ("PEP" stands for Python Enhancement Proposals). Key examples PEP 8 Python coding conventions are the following.

* Indents should be 4 spaces (not tab)
* Limit all lines to a maximum of 79 characters long blocks of text being limited to 72 characters
* Use a space after a comma
* Use a space before and after arithmetic operators

In the text editors Atom, Sublime Text 3, and Vim, you can install Linter packages that highlight areas of your code that break PEP 8 rules and tell you what the violation is.

There are fewer conventions in docstring structure, but we have developed some of our own that are outlined in the [PythonFuncs.ipynb](https://github.com/OpenRG/OG-JRC/blob/master/Tutorials/PythonFuncs.ipynb) Jupyter notebook. See especially Sections 3 and 4.


## 4. Git and GitHub tutorial

We have included a tutorial on using [Git and GitHub.com](https://github.com/OpenRG/OG-JRC/blob/master/Tutorials/git_tutorial.pdf) in the [Tutorials](https://github.com/OpenRG/OG-JRC/tree/master/Tutorials) directory of this repository. Git is a powerful version control software that comes natively installed on many machines and is widely used. GitHub.com is the most widely used online platform for hosting open source projects and integrating with Git software. Git has a significant learning curve, but it is essential for large collaborations that involve software development.


## 5. Jupyter Notebooks

[Jupyter notebooks](http://jupyter.org/) are files that end with the `*.ipynb` suffix. These notebooks are opened in a browser environment and are an open source web application that combines instructional text with live executable and modifyable code for many different programming platforms (e.g., Python, R, Julia). Jupyter notebooks are an ideal tool for teaching programming as they provide the code for a user to execute and they also provide the context and explanation for the code. We have provided a number of Jupyter notebooks in the [Tutorials](https://github.com/OpenRG/OG-JRC/tree/master/Tutorials) folder of this repository.

These notebooks used to be Python-specific, and were therefore called iPython notebooks (hence the `*.ipynb` suffix). But Jupyter notebooks now support many programming languages, although the name still pays homage to Python with the vestigal "py" in "Jupyter". The notebooks execute code from the kernel of the specific programming language on your local machine.

Jupyter notebooks capability will be automatically installed with your download of the [Anaconda distribution](https://www.continuum.io/downloads) of Python. If you did not download the Anaconda distribution of Python, you can download Jupyter notebooks separately by following the instructions on the Jupyter [install page](http://jupyter.org/install.html).


### 5.1. Opening a Jupyter notebook

Once Jupyter is installed--whether through Anaconda or through the Jupyter website--you can open a Jupyter notebook by the following steps.

1. Navigate in your terminal to the folder in which the Jupyter notebook files reside. In the case of the Jupyter notebook tutorials in this repository, you would navigate to the `~/OG-JRC/Tutorials/` directory.
2. Type `jupyter notebook` at the terminal prompt.
3. A Jupyter notebook session will open in your browser, showing the available `*.ipynb` files in that directory.
4. Double click on the Jupyter notebook you would like to open.

It is worth noting that you can also simply navigate to the URL of the Jupyter notebook file in the GitHub repository on the web (e.g., [https://github.com/OpenRG/OG-JRC/blob/master/Tutorials/PythonReadIn.ipynb](https://github.com/OpenRG/OG-JRC/blob/master/Tutorials/PythonReadIn.ipynb)). You can read the Jupyter notebook on GitHub.com, but you cannot execute any of the cells. You can only execute the cells in the Jupyter notebook when you follow the steps above and open the file from your terminal.


### 5.2. Using an open Jupyter notebook

Once you have opened a Jupyter notebook, you will find the notebook has two main types of cells: Markdown cells and Code cells. Markdown cells have formatted Jupyter notebook markdown text, and serve primarily to present context for the coding cells. A reference for the markdown options in Jupyter notebooks is found in the [Jupyter markdown documentation page](http://jupyter-notebook.readthedocs.io/en/latest/examples/Notebook/Working%20With%20Markdown%20Cells.html).

You can edit a Markdown cell in a Jupyter notebook by double clicking on the cell and then making your changes. Make sure the cell-type box in the middle of the top menu bar is set to `Markdown`. To implement your changes in the Markdown cell, type `Shift-Enter`.

A Code cell will have a `In [ ]:` immediately to the left of the cell for input. The code in that cell can be executed by typing `Shift-Enter`. For a Code cell, the  cell-type box in the middle of the top menu bar says `Code`.


### 5.3. Closing a Jupyter notebook

When you are done with a Jupyter notebook, you first save any changes that you want to remain with the notebook. Then you close the browser windows associated with that Jupyter notebook session. You must then close the local server that was opened to run the Jupyter notebook in your terminal window. On a Mac or Windows, this is done by going to your terminal window and typing `Ctrl-C` and then selecting `y` for yes and hitting `Enter`.


## 6. Python tutorials

For this training, we have included in this repository five basic Python tutorials in the [Tutorials](https://github.com/OpenRG/OG-JRC/tree/master/Tutorials) directory.

1. [PythonReadIn.ipynb](https://github.com/OpenRG/OG-JRC/blob/master/Tutorials/PythonReadIn.ipynb). This Jupyter notebook provides instruction on basic Python I/O, reading data into Python, and saving data to disk.
2. [PythonNumpyPandas.ipynb](https://github.com/OpenRG/OG-JRC/blob/master/Tutorials/PythonNumpyPandas.ipynb). This Jupyter notebook provides instruction on working with data using `NumPy` as well as Python's powerful data library `pandas`.
3. [PythonDescribe.ipynb](https://github.com/OpenRG/OG-JRC/blob/master/Tutorials/PythonDescribe.ipynb). This Jupyter notebook provides instruction on describing, slicing, and manipulating data in Python.
4. [PythonFuncs.ipynb](https://github.com/OpenRG/OG-JRC/blob/master/Tutorials/PythonFuncs.ipynb). This Jupyter notebook provides instruction on working with and writing Python functions.
5. [PythonVisualize.ipynb](https://github.com/OpenRG/OG-JRC/blob/master/Tutorials/PythonVisualize.ipynb). This Jupyter notebook provides instruction on creating visualizations in Python.
6. [PythonRootMin.ipynb](https://github.com/OpenRG/OG-JRC/blob/master/Tutorials/PythonRootMin.ipynb). This Jupyter notebook provides instruction on implementing univariate and multivariate root finders and unconstrained and constrained minimizers using functions in the [`scipy.optimize`](https://docs.scipy.org/doc/scipy/reference/optimize.html) sub-library.

To further ones Python programming skills, a number of other great resources exist.

* The [official Python 3 tutorial site](https://docs.python.org/3/tutorial/)
* [QuantEcon.net](https://lectures.quantecon.org/py/index.html) is a site run by [Thomas Sargent](http://www.tomsargent.com/) (NYU Stern) and [John Stachurski](http://johnstachurski.net/) (Australia National University). QuantEcon has a very large number of high-quality economics focused computational tutorials in Python.
* [Python computational labs](http://www.acme.byu.edu/2016-2017-materials/) of the Applied and Computational Mathematics Emphasis at Brigham Young University
* [Code Academy's Python learning module](https://www.codecademy.com/learn/python)

In addition, a number of excellent textbooks and reference manuals are very helpful and may be available in your local library. Or you may just want to have these in your own library. Lutz (2013) is a giant 1,500-page reference manual that has an expansive collection of materials targeted at beginners. Beazley (2009) is a more concise reference but is targeted at readers with some experience using Python. Despite its focus on a particular set of tools in the Python programming language, McKinney (2013) has a great introductory section that can serve as a good starting tutorial. Further, its focus on Python's data analysis capabilities is truly one of the important features of Python. Rounding out the list is Langtangen (2010). This book's focus on scientists and engineers makes it a unique reference for optimization, wrapping C and Fortran and other scientific computing topics using Python.


## References

* Beazley, David M., *Python Essential Reference*, 4th edition, Addison-Wesley (2009).
* Langtangen, Hans Petter, *Python Scripting for Computational Science*, Texts in Computational Science and Engineering, 3rd edition, Springer (2010).
* Lutz, Mark, *Learning Python*, 5th edition, O'Reilly Media, Inc. (2013)
* McKinney, Wes, *Python for Data Analysis*, O'Reilly Median, Inc. (2013)

