~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KIOPS : Krylov with Incomplete Orthogonalization Process Solver
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KIOPS is the Matlab package described in "KIOPS: A fast adaptive Krylov
subspace solver for exponential integrators" by Stéphane Gaudreault, Greg
Rainwater and Mayya Tokman.

https://doi.org/10.1016/j.jcp.2018.06.026

To install, copy the files kiops.m and expm_kiops.m
to a directory where Matlab can find them.

KIOPS is an algorithm for computing linear combinations of φ-functions that
appear in exponential integrators. This algorithm is suitable for large-scale
problems in computational physics where little or no information about the
spectrum or norm of the Jacobian matrix is known a priori. 