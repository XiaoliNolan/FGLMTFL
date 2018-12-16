# Fused Group Lasso Regularized Multi-task Feature Learning (FGL-MTFL) #

This repository contains a MATLAB implementation of the FGL-MTFL algorithm proposed in the paper [Fused Group Lasso Regularized Multi-Task Feature Learning and Its Application to the Cognitive Performance Prediction of Alzheimer’s Disease](https://link.springer.com/article/10.1007/s12021-018-9398-5).

## Overview ##

The Fused Group Lasso Regularized Multi-task Feature Learning model assumes unequal correlation of the tasks and effects of different correlated tasks on different brain regions. FGL-MTFL models the underlying structures with a graph structure within tasks and a group structure among the image features. FGL-MTFL model includes two regularization processes: (1) all tasks are
regularized by the l21-norm regularizer to capture global relationships. (2) Two local structures (graph structures within tasks and group structures within features) are considered to capture the specific task-ROI structures. For optimization, the alternating direction method of multipliers (ADMM) is employed to efficiently solve the proposed non-smooth formulation.

This code has been tested only in MATLAB in both Linux and Mac.

## How to run? ##

We created the file `fglMtfl.m` to show how to run FGL-MTFL code. FGL-MTFL first constructs a graph to model the task correlation, then optimized based on ADMM algorithm.

## Structure of the input data files ##

In order to run the code the input data files containing the training and test data must follow a specific format. 
The `FLADMM()` function, which is the core algorithm of GFL-SGL, receives 1) two matrices, **X** (covariate matrix) *n x p* with the number of *n* samples and *p* covariates, and **Y** (response matrix) *n x k* with *k* tasks; 2) group information vector with *p* covariates divided *q* disjoint groups. Note that, the number of features in each group can be different; 3) task undirected graph **G** with **V** representing the set of vertices and **E** the set of edges, which is created using only the training data.

## How to cite it? ##

If you like it and want to cite it in your papers, you can use the following:

```
#!latex

@article{liu2018fused,
  title={Fused Group Lasso Regularized Multi-Task Feature Learning and Its Application to the Cognitive Performance Prediction of Alzheimer’s Disease},
  author={Liu, Xiaoli and Cao, Peng and Wang, Jianzhong and Kong, Jun and Zhao, Dazhe},
  journal={Neuroinformatics},
  pages={1--24},
  year={2018},
  publisher={Springer}
}
```

## Have a question? ##

If you found any bug or have a question, don't hesitate to contact me:

[Xiaoli Liu]
email: `neuxiaoliliu -at- gmail -dot- com`
