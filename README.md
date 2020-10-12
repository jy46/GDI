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

(Optional) If you want to run example 4, the arbitrary network, then please download the simulation file [here](https://drive.google.com/file/d/1OOLpoqL5_SYA9FDRIifWdFgqsfy8w2Q5/view?usp=sharing) and add it to the `gdi_python` folder.

*Note*: Our GDI results on binned spike times are normalized using code from CTM-DI, which is only available for MATLAB. For the purposes of this repository, we included normalization factors that were precalculated using CTM-DI and are loaded for our two binned spike time examples. Please note that these normalizations are only appropriate for the already specified bin widths and M values in those examples.

### Usage
The `GDI.py` file contains all functions used for GDI. The core functions are:
  - `DI(X,M,B)`: Compute the pairwise (non-graphical) DI between columns of X with history parameter M (i.e. past number of samples relevant in estimation) and B bootstrap iterations. Output is DI matrix where each element is the DI from the row to the column.
  - `GDI(X,M,B)`: Compute the GDI between columns of X with history parameter M (i.e. past number of samples relevant in estimation) and B bootstrap iterations. When computing the GDI between two columns, all other columns are conditioned on. Output is GDI matrix where each element is the GDI from the row to the column.
  - `sign_inference(X,M)`: Determine the sign of the relationship between columns of X with history parameter M. First output is the partial sign matrix where each element is the sign of the relationship from the row to the column based on partial correlations, and second output is the same but based on regular correlations.
  - `GDI_mask(X,M,B,mask)`: Same as `GDI()`, however a mask in the form of a square matrix containing zeros and ones specifies which time series GDI is to be computed for as well as which time series are to be conditioned on in such GDI computations. For example, a mask with all zeros except for ones at elements (2,3), (4,3), and (5,4) would mean that GDI would only be computed from columns 2 to 3, 4 to 3, and 5 to 4 of X. Furthermore, the GDI computation from column 2 to 3 would only be conditioned on column 4, while the GDI from 4 to 3 would only be conditioned on 2, and finally the GDI from 5 to 4 would not be conditioned on any other column, making it equivalent to the DI from 5 to 4. Output is GDI matrix where each element is the GDI from the row to the column.
  
For more detail, view the header for each function in the `GDI.py` file.

### Minimal Working Example ([Jupyter Notebook](gdi_python/minimal_working_example.ipynb))
Here we construct an example where node 0 causally influences node 1 and node 1 causally influences node 2. This results in direct connections from 0 to 1 and 1 to 2 as well as an indirect connection from 0 to 2 which GDI correctly eliminates.

```python
import GDI
import numpy as np

num_samples = 5000
mean = [0,0,0]
cov  = np.eye(3)

X = np.random.multivariate_normal(mean, cov, num_samples)

X[2:,1] = X[2:,1] + X[:-2,0]
X[2:,2] = X[2:,2] + X[:-2,1]

M = 4
B = 10

X_DI = GDI.DI(X,M,B)
X_GDI = GDI.GDI(X,M,B)

print(X_DI)
print(X_GDI)
```

An example run produces a DI matrix of:
<pre>
[[ 0.          0.36758208  <b>0.15418205</b>]
 [ 0.04643047  0.          0.46776479]
 [-0.00124419  0.01545924  0.        ]]
</pre>
where the incorrectly identified indirect connection from 0 to 2 is bolded. 

Then GDI matrix then looks like:
<pre>
[[ 0.          0.29086685 <b>-0.00878638</b>]
 [-0.03832781  0.          0.21460593]
 [ 0.01501656 -0.00237674  0.        ]]
</pre>
which shows that GDI eliminated the indirect connection from node 0 to 2 that was incorrectly identified by DI.

## MATLAB
### Installation
These steps are required (although all steps regarding CTM-DI can be ignored if not running arb.m or cpg.m):
1. Download [CCMI](https://github.com/sudiptodip15/CCMI) and [CTM-DI](http://www.ece.rice.edu/neuroengineering/). Extract all files from the CIT folder of CCMI, and place them in the ccdi_mat directory of our code. The CTM_DI_package folder of CTM-DI should be placed at the same directory level as the ccdi_mat directory of our code. Ensure that you have the necessary Python packages installed as listed in the prior Python section.

2. Modify a few lines of the CCMI code to include the CMI estimate for each bootstrap iteration before averaing across bootstrap iterations occurs. The specific lines that we want to modify are right around the return statement at the end of the definition for `get_cmi_est()` in the file `CCMI.py`. When we downloaded CCMI, the lines to be considered were 102 to 106, which appear as:
```python
            cmi_est = I_xyz - I_xz
        else:
            raise NotImplementedError

        return cmi_est
```

Change those lines to include the CMI estimates before averaging over bootstrap iterations:
```python
            cmi_est_list = I_xyz_list - I_xz_list
            cmi_est = I_xyz - I_xz
        else:
            raise NotImplementedError

        return cmi_est, cmi_est_list
```

3. Replace the contents of the file `CTM_DI_package/Supporting_functions/CTM_spiketime_wrapper.m` with the following:
```matlab
    function [CMatrix, w_i, HMatrix] = CTM_spiketime_wrapper(sig,M, IRI, causal)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % note this function takes vectors of spike times per channel as input.
    % When a channel has less spike than others, zero pad the rest of the
    % entries
    %
    % sig is an array where rows are spike times in ms and columns are neurons
    % M is the maximum depth that the tree can be (maximum number of past bins used)
    % IRI is the bin_width to use in ms
    % causal should be 1 (causal definition of DI)
    % SPIKE TIMES SHOULD BE IN MS!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        channel = size(sig,2);
        CMatrix = zeros(channel,channel,length(IRI));

        parfor cnt = 1:length(IRI)
            cnt
            [CMatrix(:,:,cnt) w_i(:,:,:,cnt) HMatrix(:,:,cnt)]= ...
                Connect_CTM(spiketime2train(sig,IRI(cnt)),M,causal);
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    return
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
```

4. Modify `CTM_DI_package/CTM_DI/Connect_CTM.m` as follows:
- Change the function definition to:
```matlab
  function [CMatrix,w_i,HMatrix] = Connect_CTM(Sig,D,causal,sumspike,sig_label)
```
- Add this line in between lines 68 and 69 (which should be between `weight_matrix = zeros(channel);` and `for cnt = 1:length(row)`):
```matlab
  w_i = nan(channel,channel,D+1);
```
- Just after what should now be line 75 (`weight_matrix(col(cnt),row(cnt))=sum(weights2(1+causal:end));`), add these two lines just before `end`:
```matlab
  w_i(row(cnt),col(cnt),1:length(weights1)) = weights1;
  w_i(col(cnt),row(cnt),1:length(weights2)) = weights2;    
```
- Just after the definition of CMatrix (`CMatrix = zeros(channel,channel);`) which should be at line 83 now, add this right after:
```matlab
  HMatrix = zeros(channel,channel);
```
- Just after `CMatrix(col(cnt), row(cnt)) = DI21/H1;` which should now be at line 117, add these two lines:
```matlab
  HMatrix(row(cnt), col(cnt)) = H2;
  HMatrix(col(cnt), row(cnt)) = H1;
```

5. Since the core of this toolbox relies on the CCMI implementation which is written in Python, you must insert your system path and python path in the `gdi_matlab/python_path_script.m` file. This script is called by deeper functions to access python. This means copying the terminal output for the command `echo $PATH` and putting it in between the '' for system_path in the `gdi_matlab/python_path_script.m` file, and then also copying the output for the command `which python` and putting it in between the '' for the python_path in the `gdi_matlab/python_path_script.m` file.

6. (Optional) If you want to run example 4, the arbitrary network, then please download the simulation file [here](https://drive.google.com/file/d/1OOLpoqL5_SYA9FDRIifWdFgqsfy8w2Q5/view?usp=sharing) and place it in the `gdi_matlab` folder.

### Usage
The `ccdi_mat` folder contains all of the files/functions for GDI. The core functions are:
  - `di_compute(X,M,C,B)`: Computes the DI or GDI between each column of X using the history parameter M (i.e. past number of samples relevant in estimation) and B bootstrap iterations. C is 0 (compute DI) or 1 (compute GDI).
  - `di_compute_pair(X,M,C,B,pairs)`: Same as `di_compute()`, but computes DI/GDI only for the specified pairs.
  - `di_compute_post(DI_uncond,thresh,M,X,B)`: Computes GDI based on thresholding the (non-grapical) DI values. GDI will only be computed between channels with DI values >= thresh, and will only be conditioned on channels with DI values >= thresh. If thresholding means there are no other channels to condition on for a particular GDI analysis, then that GDI analysis will not be performed and the DI value will be taken to be the GDI value.
  - `sign_inference(X,M)`: Determine the sign of the relationship between columns of X with history parameter M. First output is the partial sign matrix where each element is the sign of the relationship from the row to the column based on partial correlations, and second output is the same but based on regular correlations.
  
For more detail, view the header for each function/file in the `ccdi_mat` folder.

### Minimal Working Example ([MATLAB](gdi_matlab/minimal_working_example.m))
Here we construct an example where node 1 causally influences node 2 and node 2 causally influences node 3. This results in direct connections from 1 to 2 and 2 to 3 as well as an indirect connection from 1 to 3 which GDI correctly eliminates.

```matlab
addpath ccdi_mat

num_samples = 5000;
mean_vec = [0,0,0];
cov_mat  = eye(3);

X = mvnrnd(mean_vec, cov_mat, num_samples);

X(3:end,2) = X(3:end,2) + X(1:(end-2),1);
X(3:end,3) = X(3:end,3) + X(1:(end-2),2);

M = 4;
B = 10;

X_DI = di_compute(X,M,0,B)
X_GDI = di_compute(X,M,1,B)
```

An example run produces a DI matrix of:
<pre>
       NaN    0.3569    <b>0.1730</b>
   -0.0298       NaN    0.4008
   -0.0060   -0.1157       NaN
</pre>
where the incorrectly identified indirect connection from 1 to 3 is bolded. 

Then GDI matrix then looks like:
<pre>
       NaN    0.2392    <b>-0.0166</b>
    0.0490       NaN    0.3186
   -0.0438   -0.1290       NaN
</pre>
which shows that GDI eliminated the indirect connection from node 1 to 3 that was incorrectly identified by DI.

## Examples Used in Paper
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

## Licensing
We have adopted the GPLv2 license for this toolbox (see LICENSE and GPLv2_note.txt files).
