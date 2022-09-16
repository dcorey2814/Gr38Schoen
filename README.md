# Gr38Schoen
Code used in the paper <a href="https://arxiv.org/abs/2206.14993"> The Grassmannian of 3-planes in C^8</a> is sch\"on by Daniel Corey and Dante Luber. This code works with OSCAR version 0.10.1. Note: it may be essential to use this version of OSCAR.  

**Convention** $\operatorname{TGr}_0(3,8)$ is the support of a fan in $\wedge^{3}\mathbb{R}^{8} \cong \mathbb{R}^{\binom{[8]}{3}}$. We order the coordinates of $\mathbb{R}^{\binom{[8]}{3}} \cong \mathbb{R}$ lexicographically.  In the paper, we use the MIN convention, however, the data of $\operatorname{TGr}_{0}(3,8)$ is given to us with respect to the MAX convention. Thus, when forming the associated subdivision of $\Delta(3,8)$, one must *negate* the vector. 

## The seconday fan structure $\operatorname{TGr}_{0}(3,8)$: Raw data

The data of $\operatorname{TGr}_{0}(3,8)$ was shared with us by Benjamin Schr\"oter, which was computed in the paper <a href="https://arxiv.org/abs/2003.13752"> Parallel Computation of tropical varieties, their positive part, and tropical Grassmannians</a> by Dominik Bendle, Janko Boehm, Yue Ren, and Benjamin Schr\"oter. This is in the files ```group38```, ```GrRays.data```, and ```ConesDrOfGr.data```. 

```group38``` is a polymake data file containing an ```Array<Array<Int>>```. Denote by $S_{n}$ the symmetric group on $[n]$.  This file records the subgroup of $S_{56}$  isomorphic to  $S_8$ induced by the action of $S_8$ on $\binom{[8]}{3}$ given by

$\sigma \cdot (i,j,k ) = (\sigma(i),\sigma(j),\sigma(k))$

Each of these permutation in $S_{56}$ is recorded in the standard ``one-line'' notation of a permutation, i.e., $\sigma = (\sigma(1), \sigma(2), \ldots, \sigma(56))$.

```GrRays.data``` is a polymake data file containing a ```Matrix<Rational,NonSymmetric>``` recording the rays of $\operatorname{TGr}_{0}(3,8)$. Each row of this matrix is a ray of $\operatorname{TGr}_{0}(3,8)$ with its *Gr\"obner* fan structure; only the first 12 are needed for the secondary fan structure. 

```ConesDrOfGr.data``` is a text file recording the maximal cones of $\operatorname{TGr}_{0}(3,8)$ with its secondary fan structure. Each line represents a cone as a space separated list of symbols ```r#s```. Here, ```r``` represents a row of ```GrRays.data```, and ```s``` a row of ```group38```. Thus, the symbol ```r#s``` mean ``the r-th ray whose coordinates are permuted by the s-th permutation in ```group38```.'' 

## The seconday fan structure $\operatorname{TGr}_{0}(3,8)$: The intermediate cones


## Verifying propositions in Section 6

The notebooks ```G2.ipynb```, ```G3.ipynb```, ```G4.ipynb```, ```G5.ipynb```, ```G6.ipynb``` contain the code used in the proof of Propositions 6.7, 6.13, 6.14, 6.16, and 6.19, respectively. These rely on the functions in the files contained in the ```src``` directory. The documentation for these functions is in the notebook ```functionDocumentation.ipynb```. Instructions on full verifications and examples are also provided. 



