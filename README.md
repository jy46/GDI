# Graphical Directed Information (GDI)
Graphical directed information (GDI) is a model-free technique that infers causal influences between pairs of time series and in particular captures unique influences between pairs by conditioning on other time series . 

The directed information (DI) from a time series X(t) to another time series Y(t) is the mutual information (MI) between the past of X(t) and the present of Y(t) conditioned on the past of Y(t), which means DI quantifies the causal influence exclusively from X(t) to Y(t). 

GDI is an extension of DI that conditions not only on the past of Y(t), but also conditions on the pasts of other time series W(t), Z(t), etc. Conditioning on other time series enables the identification and quantification of a direct connection from one time series to another as well as the possibility for eliminating indirect connections.

We offer Python and MATLAB implementations of GDI, which both include five different example GDI analyses as well as instructions on how to install both.

## Python
### Installation
GDI was tested using Python 3.7, and requires the following packages:
  - NumPy
  - seaborn
  - SciPy
  - TensorFlow (<v2.0)

### Usage
The

## MATLAB
### Installation
The key file

### Usage
The

https://github.com/sudiptodip15/CCMI

## Examples
We included five different example anlyses in our code.

### 1. Scaling
Scaling

### 2. Gaussian Network
Gaussian

### 3. Nonlinear Network
Nonlinear

### 4. Arbitrary Network
Arbitrary

### 5. CPG Network
CPG

## Method
GDI takes advantage 

## References
