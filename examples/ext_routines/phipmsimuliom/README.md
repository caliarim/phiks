# Phipm_simul_iom

Evaluates a linear combination of the phi functions evaluated at the scaled matrix, t*A, acting on a set of input vectors.

The evaluation expresses eveything in terms of the highest order phi function and evaluates the action of this on a vector using a Krylov technique and then computes w using the recurrence relation. 

The size of the Krylov subspace is changed dynamically during the integration. The Krylov subspace is computed using the Arnoldi process.

This file is a modified version of phipm_iom.m.

The modifications were performed by Daniel R. Reynolds and Vu Luan, Mathematics Department, Southern Methodist University, Fall 2017.  This routine was used within the manuscript, "Further development of efficient and accurate time integration schemes for meterological models," arXiv:1805.02144 [math.NA].
