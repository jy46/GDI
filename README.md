# Graphical Directed Information (GDI)
Graphical directed information (GDI) is a model-free technique that infers causal influences between pairs of time series that captures unique influences between pairs by conditioning on other time series . 

The directed information (DI) from a time series X(t) to another time series Y(t) is the mutual information (MI) between the past of X(t) and the present of Y(t) conditioned on the past of Y(t), which means DI quantifies the causal influence exclusively from X(t) to Y(t). 

GDI is an extension of DI that conditions not only on the past of Y(t), but also conditions on the pasts of other time series W(t), Z(t), etc. Conditioning on other time series enables the identification and quantification of a direct connection from one time series to another as well as the possibility for eliminating indirect connections.

We offer both Python and MATLAB implementations of GDI, along 


## Python Implementation

## MATLAB Implementation

https://github.com/sudiptodip15/CCMI

