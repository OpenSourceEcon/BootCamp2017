# OSM Lab: Computational Methods Section Materials

This directory in the repository contains all the materials for the computational methods section of the OSM Lab Boot Camp.


## Prerequisite and tutorial resources

We expect students in the Boot Camp to be jumping into Python at a level beyond an absolute beginner. Some Great resources for getting to that point include the following Jupyter Notebooks in the [Tutorials](https://github.com/OpenSourceMacro/BootCamp2017/tree/master/Tutorials) folder of this repository.

1. [PythonReadIn.ipynb](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Tutorials/PythonReadIn.ipynb). This Jupyter notebook provides instruction on basic Python I/O, reading data into Python, and saving data to disk.
2. [PythonNumpyPandas.ipynb](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Tutorials/PythonNumpyPandas.ipynb). This Jupyter notebook provides instruction on working with data using `NumPy` as well as Python's powerful data library `pandas`.
3. [PythonDescribe.ipynb](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Tutorials/PythonDescribe.ipynb). This Jupyter notebook provides instruction on describing, slicing, and manipulating data in Python.
4. [PythonFuncs.ipynb](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Tutorials/PythonFuncs.ipynb). This Jupyter notebook provides instruction on working with and writing Python functions.
5. [PythonVisualize.ipynb](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Tutorials/PythonVisualize.ipynb). This Jupyter notebook provides instruction on creating visualizations in Python.
6. [PythonRootMin.ipynb](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Tutorials/PythonRootMin.ipynb). This Jupyter notebook provides instruction on implementing univariate and multivariate root finders and unconstrained and constrained minimizers using functions in the [`scipy.optimize`](https://docs.scipy.org/doc/scipy/reference/optimize.html) sub-library.

We also recommend the ["Intro to Python" lab](http://www.acme.byu.edu/wp-content/uploads/2016/08/PythonIntro.pdf) from Brigham Young University's Math Department's Applied and Computational Math Emphasis (ACME) as well as the ["An Introductory Example"](https://lectures.quantecon.org/py/python_by_example.html) and ["Python Essentials"](https://lectures.quantecon.org/py/python_essentials.html) lectures from [QuantEcon](https://lectures.quantecon.org/py/).


## Schedule

The computational methods lab sessions for the OSM Lab will be held from 8:00-11:50am, Tuesday and Thursday in Saieh Hall, Room 247. The lab files in the schedule under the "Materials column" that start with "ACME" come from the Brigham Young University Math Department's [Applied and Computational Math Emphasis (ACME program)](http://www.acme.byu.edu/). These computational labs are open source. We cover only a subset of these excellent applied math Python labs, which are available in their entirety at [http://www.acme.byu.edu/2016-2017-materials/](http://www.acme.byu.edu/2016-2017-materials/). We highly recommend that you take time after the Boot Camp to work through some of the other labs that are available to you.

### Week 1

| Date | Day | Topic | Instructor | Materials | Problem Set |
|:---:|:---:|:--- |:--- |:--- | --- |
6-19  | M   |     |     |     |     |
6-20  | T   |     | Justin Gardiner | [ACME: Intro to NumPy](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk1_PyIntro/NumpyIntro.pdf) | [Comp Prob Set 1](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk1_PyIntro/PyIntro_probset.pdf) |
|     |     |     |                 | [ACME: Standard Library](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk1_PyIntro/StandardLibrary.pdf) | due T, 6-27, 8am |
|     |     |     |                 | [ACME: Unit Testing](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk1_PyIntro/Vol1B-Testing-2017.pdf) |   |
6-21  | W   |     |         |          |    |
6-22  | Th  |     | Justin Gardiner | [ACME: Object Oriented Programming](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk1_PyIntro/OOP.pdf) |   |
|     |     |     |                 | [ACME: Exceptions and File I/O](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk1_PyIntro/Exceptions.pdf) |   |
6-23  | F   |     |     |     |     |

### Week 2

| Date | Day | Topic | Instructor | Materials | Problem Set |
|:---:|:---:|:--- |:--- |:--- | --- |
6-26  | M   |     |     |     |     |
6-27  | T   | Visualizations | Justin Gardiner | [Visualizations Notebook](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Tutorials/PythonVisualize.ipynb) | [Comp Prob Set 2](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk2_VisPandas/VisPandas_probset.pdf) |
|     |     | and Pandas  |    | [ACME: Intro to Matplotlib](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk2_VisPandas/PlottingIntro.pdf) | due W, 7-5, 8am |
|     |     |     |    | [ACME: Data Visualization](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk2_VisPandas/Vol1A-DataVisualization-2016.pdf) |   |
|     |     |     |    | [ACME: Pandas 1](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk2_VisPandas/Vol3A-Pandas1-2016.pdf) |   |
6-28  | W   |     |         |          |     |
6-29  | Th  | Visualizations | Justin Gardiner | [Pandas Notebook](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Tutorials/PythonNumpyPandas.ipynb) |   |
|     |     | and Bokeh      |          | [ACME: Pandas 3](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk2_VisPandas/Vol3A-Pandas3-2016.pdf) |    |
|     |     |      |          | [ACME: Pandas 4](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk2_VisPandas/Vol3A-Pandas4-2016.pdf) |    |
|     |     |     |    | [ACME: Bokeh](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk2_VisPandas/Vol3A-Bokeh-2016.pdf) |   |
6-30  | F   |     |     |     |     |

### Week 3

| Date | Day | Topic | Instructor | Materials | Problem Set |
|:---:|:---:|:--- |:--- |:--- | --- |
7-3  | M    | NO CLASSES: HOLIDAY | NO CLASSES: HOLIDAY | NO CLASSES: HOLIDAY  |  |
7-4  | T    | NO CLASSES: HOLIDAY | NO CLASSES: HOLIDAY | NO CLASSES: HOLIDAY  |  |
7-5  | W    |     |         |          |   |
7-6  | Th   |  Matrix Decomposition | Justin Gardiner | [ACME: QR Decomp 1](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk3_Decomp/Vol1A-QR1-Decomposition-2016.pdf) | [Comp Prob Set 3](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk3_Decomp/Decomp_probset.pdf) |
|     |     | and Compute Speed |          | [ACME: QR Decomp 2](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk3_Decomp/Vol1A-QR2-LeastSqEigs-2016.pdf) | due T, 7-11, 8am |
|     |     |      |          | [ACME: SVD Image Compress](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk3_Decomp/Vol1A-SVD-2016.pdf) |  |
|     |     |      |          | [ACME: Drazin Inverse](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk3_Decomp/Vol1B-Drazin-2017.pdf) |  |
|     |     |      |          | [QuantEcon: Need For Speed](https://lectures.quantecon.org/py/need_for_speed.html) |  |
7-7  | F    |     |     |     |  |

### Week 4

| Date | Day | Topic | Instructor | Materials | Problem Set |
|:---:|:---:|:--- |:--- |:--- | --- |
7-10  | M   |     |     |     |     |
7-11  | T   | Matrix conditions | Justin Gardiner | [ACME: PageRank Algorithm](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk4_DifIntOpt/ACME_PageRank.pdf) | [Comp Prob Set 4](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk4_DifIntOpt/DifIntOpt_probset.pdf)  |
|     |     | Numerical Diff. and Int. |           | [ACME: Conditioning and Stability](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk4_DifIntOpt/ACME_CondStab.pdf) | due Th, 7-20, 8am |
|     |     |     |      | [ACME: Numerical Differentiation](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk4_DifIntOpt/ACME_NumDiff.pdf) |  |
|     |     |     |      | [Evans: Numerical Integration](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk4_DifIntOpt/NumIntegr_Evans.pdf) |  |
7-12  | W   |     |     |     |    |
7-13  | Th  | Root finding | Justin Gardiner | [ACME: Simplex Method](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk4_DifIntOpt/ACME_Simplex.pdf) |  |
|     |     |     |      | [ACME: Line Search Methods](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk4_DifIntOpt/ACME_LineSrch.pdf) |  |
|     |     |     |      | [ACME: Newton's Method](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk4_DifIntOpt/ACME_Newtons.pdf) |  |
|     |     |     |      | [ACME: Iterative Solvers](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk4_DifIntOpt/ACME_IterSolvers.pdf) |  |
7-14  | F   |     |     |     |     |

### Week 5

| Date | Day | Topic | Instructor | Materials | Problem Set |
|:---:|:---:|:--- |:--- |:--- | --- |
7-17  | M   |     |     |     |     |
7-18  | T   | Minimization | Justin Gardiner | [ACME: Interior Point, Linear Programs](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk4_DifIntOpt/ACME_IntPtLin.pdf) |  |
|     |     |     |      | [ACME: Interior Point, Quadratic Programs](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk4_DifIntOpt/ACME_IntPtQuad.pdf) |  |
|     |     |     |      | [ACME: Newton and Quasi Newton Methods](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk4_DifIntOpt/ACME_QuasNewt.pdf) |  |
|     |     |     |      | [ACME: Scipy.optimize](https://github.com/OpenSourceMacro/BootCamp2017/blob/master/Computation/Wk4_DifIntOpt/ACME_ScipyOpt.pdf) |  |
7-19  | W   |     |     |     |     |
7-20  | Th  | Adaptive sparse grids, Smolyak | Simon Scheidegger |  | Comp Prob Set 5 |
7-21  | F   |     |     |     | due T, 7-25, 8am    |

### Week 6

| Date | Day | Topic | Instructor | Materials | Problem Set |
|:---:|:---:|:--- |:--- |:--- | --- |
7-24  | M   |     |     |     |     |
7-25  | T   | HPC, parallel computing | Simon Scheidegger |  | Comp Prob Set 6 |
7-26  | W   |     |     |     | due T, 8-1, 8am |
7-27  | Th  | HPC, parallel computing | Simon Scheidegger |  |  |
7-28  | F   |     |     |     |     |

### Week 7

| Date | Day | Topic | Instructor | Materials | Problem Set |
|:---:|:---:|:--- |:--- |:--- | --- |
7-31 | M   |      |     |     |     |
8-1  | T   | HPC, parallel computing | Simon Scheidegger |  | Comp Prob Set 7  |
8-2  | W   |      |     |     | due F, 8-4, 8am |
8-3  | Th  | HPC, parallel computing | Simon Scheidegger |  |  |
8-4  | F   | Conclusion: Hwk due |  |  |  |


## References

* Humpherys, Jeffrey, "[Computational Labs for *Foundations of Applied Mathematics*](http://www.acme.byu.edu/2016-2017-materials/)" (2017).

