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
  - `DI(X,M,B)`: Compute the pairwise (non-graphical) DI between columns of X with history parameter M (i.e. past number of samples relevant in estimation) and B bootstrap iterations. Output is DI matrix where each element is the DI from the row to the column.
  - `GDI(X,M,B)`: Compute the GDI between columns of X with history parameter M (i.e. past number of samples relevant in estimation) and B bootstrap iterations. When computing the GDI between two columns, all other columns are conditioned on. Output is GDI matrix where each element is the GDI from the row to the column.
  - `sign_inference(X,M)`: Determine the sign of the relationship between columns of X with history parameter M. First output is the partial sign matrix where each element is the sign of the relationship from the row to the column based on partial correlations, and second output is the same but based on regular correlations.
  - `GDI_mask(X,M,B,mask)`: Same as `GDI()`, however a mask in the form of a square matrix containing zeros and ones specifies which time series GDI is to be computed for as well as which time series are to be conditioned on in such GDI computations. For example, a mask with all zeros except for ones at elements (2,3), (4,3), and (5,4) would mean that GDI would only be computed from columns 2 to 3, 4 to 3, and 5 to 4 of X. Furthermore, the GDI computation from column 2 to 3 would only be conditioned on column 4, while the GDI from 4 to 3 would only be conditioned on 2, and finally the GDI from 5 to 4 would not be conditioned on any other column, making it equivalent to the DI from 5 to 4. Output is GDI matrix where each element is the GDI from the row to the column.
  
For more detail, view the header for each function in the `GDI.py` file.

### Minimal Working Example
The

```python
import GDI
import numpy as np

num_samples = 5000
mean = [0,0,0]
cov  = [[1, 0.75, 0],[0.75, 1, 0],[0, 0, 1]]

X = np.random.multivariate_normal(mean, cov, num_samples)
X[3:,1] = X[:-3,1]+X[:-3,2]

M = 3
B = 2
X_GDI = GDI(X,M,B)
X_DI = DI(X,M,B)
```

## MATLAB
### Installation
The key file

GDI also requires [CCMI](https://github.com/sudiptodip15/CCMI) and [CTM-DI](http://www.ece.rice.edu/neuroengineering/).

### Usage
The

### Minimal Working Example
The

```matlab
for ii=1:10
  help(ii) = 1
end
```

## Examples
We included five different example analyses in our code which correspond to the five results figures in our paper, and links in parentheses go to the files within our repository corresponding to those examples. 

### 1. Scaling ([MATLAB](gdi_matlab/scaling.m), [Jupyter Notebook](gdi_python/scaling.ipynb))
GDI's performance with regard to sample size, the number of dimensions being conditioned on, and the number of bootstrap iterations (see Method section below) is analyzed.

### 2. Gaussian Network ([MATLAB](gdi_matlab/gaussian.m), [Jupyter Notebook](gdi_python/gaussian.ipynb))
GDI is applied to a Gaussian network, which consists of causal influences between nodes that have their own Gaussian noise. The analytic solution for GDI is known for this network, and the accuracy of our GDI estimates is compared with the derived values.

### 3. Nonlinear Network ([MATLAB](gdi_matlab/nonlinear.m), [Jupyter Notebook](gdi_python/nonlinear.ipynb))
GDI is applied to a nonlinear network, which has the same structure as the Gaussian network however source nodes now follow uniform distributions and causal influence involves a squared relationship.

### 4. Arbitrary Network ([MATLAB](gdi_matlab/arb.m), [Jupyter Notebook](gdi_python/arb.ipynb))
GDI is applied to binned spike times produced by an abritrary model of a network of neurons.

### 5. CPG Network ([MATLAB](gdi_matlab/cpg.m), [Jupyter Notebook](gdi_python/cpg.ipynb))
GDI is applied to binned spike times produced by a model of the central pattern generator (CPG) in *Aplyisa*.
