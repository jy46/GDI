# Graphical Directed Information (GDI)
This code is provided in tandem with our paper on graphical directed information (GDI), which is a model-free technique that infers causal influences between pairs of time series and in particular captures unique influences between pairs by conditioning on other time series. 

The directed information (DI) from a time series X(t) to another time series Y(t) is the mutual information (MI) between the past of X(t) and the present of Y(t) conditioned on the past of Y(t), which means DI quantifies the causal influence exclusively from X(t) to Y(t). 

GDI is an extension of DI that conditions not only on the past of Y(t), but also conditions on the pasts of other time series W(t), Z(t), etc. Conditioning on other time series enables the identification and quantification of a direct connection from one time series to another as well as the possibility for eliminating indirect connections.

We note that one must choose a history parameter which controls the number of past samples that are considered when using the pasts of X(t) and the time series being conditioned on, which is referred to as M in our code.

We offer Python and MATLAB implementations of GDI, which both include five different example GDI analyses as well as instructions on how to install both. We also include the [SNNAP](https://med.uth.edu/nba/snnap/) files which were used to simulate networks of neurons.

For all relevant references, please see our paper.

## Python
### Installation
GDI was tested using Python 3.7, and requires the following packages:
  - NumPy
  - seaborn
  - SciPy
  - TensorFlow (<v2.0)

GDI also requires [CCMI](https://github.com/sudiptodip15/CCMI). Extract all files from the CIT folder of CCMI, and place them in the gdi_python directory of our code.

*Note*: Our GDI results on binned spike times are normalized using code from CTM-DI, which is only available for MATLAB. For the purposes of this repository, we included normalization factors that were precalculated using CTM-DI and are loaded for our two binned spike time examples. Please note that these normalizations are only appropriate for the already specified bin widths and M values in those examples.

### Usage
The `GDI.py` file contains all functions used for GDI. The core functions are:
  - `DI(X,M,B)`: Compute the pairwise (non-graphical) DI between columns of X with history parameter M (i.e. past number of samples relevant in estimation) and B bootstrap iterations.
  - `GDI(X,M,B)`:
  - `sign_inference(X,M)`:
  - `GDI_mask(X,M,B,mask)`:
For more detail, view the headers for each function in the `GDI.py` file.

### Minimal Working Example
The

## MATLAB
### Installation
The key file

GDI also requires [CCMI](https://github.com/sudiptodip15/CCMI) and [CTM-DI](http://www.ece.rice.edu/neuroengineering/).

### Usage
The

### Minimal Working Example
The

## SNNAP Files
We include 

## Examples
We included five different example analyses in our code which correspond to the five results figures in our paper:

### 1. Scaling (scaling.m, scaling.ipynb)
GDI's performance with regard to sample size, the number of dimensions being conditioned on, and the number of bootstrap iterations (see Method section below) is analyzed.

### 2. Gaussian Network (gaussian.m, gaussian.ipynb)
GDI is applied to a Gaussian network, which consists of causal influences between nodes that have their own Gaussian noise. The analytic solution for GDI is known for this network, and the accuracy of our GDI estimates is compared with the derived values.

### 3. Nonlinear Network (nonlinear.m, nonlinear.ipynb)
GDI is applied to a nonlinear network, which has the same structure as the Gaussian network however source nodes now follow uniform distributions and causal influence involves a squared relationship.

### 4. Arbitrary Network (arb.m, arb.ipynb)
GDI is applied to binned spike times produced by an abritrary model of a network of neurons.

### 5. CPG Network (cpg.m, cpg.ipynb)
GDI is applied to binned spike times produced by a model of the central pattern generator (CPG) in *Aplyisa*.
