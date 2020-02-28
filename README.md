
# Nonlinear Dimensionality Reduction for Clustering

## Introduction
Clusters defined in low dimensional manifolds can have highly nonlinear
structure, which can cause linear dimensionality reduction methods to fail. A
number of established nonlinear dimensionality reduction methods use a graph
representation of the data to define the transformation (embedding) to a low
dimensional space. We employ the Isometric mapping (Isomap) methodology
to investigate conditions under which the separability of clusters defined in
manifolds is guaranteed in the low dimensional representation.

The proposed algorithm uses the acronym "i-DivClu" (Isometric mapping for Divisive Clustering).
It is based on the idea that a suitably defined one-dimensional
representation is sufficient for constructing cluster boundaries that split the data
without breaking any of the clusters. Repeating the procedure recursively provides a
theoretically justified and efficient non-linear clustering technique.

We provide two variations of this methodology ,"i-DivClu-M" for maximum margin clustering and "i-DivClu-D"
for density based clustering.

The provided scripts can be used to compare experimentaly "i-DivClu" with popular methods also considered
able to handle non linearity in a number of real world dataset with respect to clustering efficiency and
computational cost.
![Image of]
(https://github.com/usersotiris/nonlinearclustering/blob/master/myexampleplot4.pdf)

## License
This project is licensed under the BSD-3-Clause License - see the LICENSE file for details
