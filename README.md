# Numerics supporting the findings in "Robust self-testing of two-qubit states"

arXiv-link: ...

The code in this repository supports the finding in the paper *Robust self-testing of two-qubit states* by Tim Coopmans, Jędrzej Kaniewski and Christian Schaffner (2018, arXiv-link). Using this code one can numerically verify the positivity of the operator

<dl>
T<sub>&alpha;</sub>(a, b) = K<sub>&alpha;</sub>(a, b) - s<sub>&alpha;</sub> W<sub>&alpha;</sub>(a, b) - &mu;<sub>&alpha;</sub> I
</dl>

over a discretization of (a,b) &isin; \[0, &pi;/4\] &times; \[0, &pi;/2\] (see also the definition of T<sub>&alpha;</sub>(a, b) on page 14 of the article).


## Downloading the code

In order to download the code, click the green button on the right top of this web page with the text **Clone or download**. 

For downloading the code as a compressed file, click **Download ZIP**. 

An alternative that uses the command line is to copy the url [https://github.com/timcp/Self-Testing_Pure_TwoQubit_States.git](https://github.com/timcp/Self-Testing_Pure_TwoQubit_States.git). Then open a terminal and execute the command
```
git clone <copied_url>
```
to save the code on your computer.

## Running the code

The code has been written in [Python](https://www.python.org/) and can be run from any Python development environment (such as [PyCharm](https://en.wikipedia.org/wiki/PyCharm)). One can also run the code directly in a terminal by first navigating to the folder where you stored the code
```
cd C://path/where/you/downloaded/the/code
```
followed by executing the command
```
python3 numerical_evidence_CKS2018_bounds.py
```

Note: the code only works with Python 3 instead of Python 2. To check your Python version, open a terminal, type `python -V` and press enter. The output should be something like ```Python 3.x.x```.

## File overview

**Definitions**

 - *tiltedCHSH.py*: definition of the tilted CHSH operators and its quantum- and classical value.
 - *CKS2018.py*: the definitions from the paper *Robust self-testing of two-qubit states* as specified above (arXiv-link).

**Numerics**

 - *output*: folder where the plots and data produced by the code will be put.
 - *numerical_evidence_CKS2018_bounds.py*: code for numerically verifying the positivity of the operator T<sub>&alpha;</sub>(a, b) as defined above.

**Miscellaneous**

 - *README.md*: this file.
