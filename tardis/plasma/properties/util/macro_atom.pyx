# module for fast macro_atom calculations

# cython: profile=False
# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False
# cython: language_level=3

import numpy as np

cimport numpy as np

ctypedef np.int64_t int_type_t

from astropy import constants


cdef extern from "math.h":
    double exp(double)


cdef double h_cgs = constants.h.cgs.value
cdef double c = constants.c.cgs.value
cdef double kb = constants.k_B.cgs.value
cdef double inv_c2 = 1 / (c ** 2)


def calculate_transition_probabilities(
        double [:] transition_probability_coef,
        double [:, ::1] beta_sobolev, double [:, ::1] j_blues,
        double [:, ::1] stimulated_emission_factor,
        int_type_t [:] transition_type,
        int_type_t [:] lines_idx,
        int_type_t [:] block_references,
        double [:, ::1] transition_probabilities,
        normalize):

    cdef int i, j, k, line_idx
    cdef np.ndarray[double, ndim=1] norm_factor = np.zeros(transition_probabilities.shape[1])

    for i in range(transition_probabilities.shape[0]):
        line_idx = lines_idx[i]
        for j in range(transition_probabilities.shape[1]):
            transition_probabilities[i, j] = (transition_probability_coef[i] * beta_sobolev[line_idx, j])
        if transition_type[i] == 1:
            for j in range(transition_probabilities.shape[1]):
                transition_probabilities[i, j] *= stimulated_emission_factor[line_idx, j] * j_blues[line_idx, j]

    if normalize:
        for i in range(block_references.shape[0] - 1):
            for k in range(transition_probabilities.shape[1]):
                norm_factor[k] = 0.0
            for j in range(block_references[i], block_references[i + 1]):
                for k in range(transition_probabilities.shape[1]):
                    norm_factor[k] += transition_probabilities[j, k]
            for k in range(transition_probabilities.shape[1]):
                if norm_factor[k] != 0.0:
                    norm_factor[k] = 1 / norm_factor[k]
                else:
                    norm_factor[k] = 1.0
            for j in range(block_references[i], block_references[i + 1]):
                for k in range(0, transition_probabilities.shape[1]):
                    transition_probabilities[j, k] *= norm_factor[k]
