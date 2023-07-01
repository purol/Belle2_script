#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# usage: basf2 MakeNtuple_multi.py "./20210402/evt-0.mdst"
# last: 2021-10-29-00

import numpy as np

Eff_matrix = np.array([[0.000358195,0.000010985,0.000103152], [0.000021071,0.000226903,0.000016230], [0.000000000,0.000006966,0.000033843]])

Eff_matrix_inv = np.linalg.pinv(Eff_matrix)
#print(Eff_matrix_inv)

#print(Eff_matrix_inv@Eff_matrix)

Nevt_data = np.array([[12.676286],[11.529586],[4.390777]])

print(Eff_matrix_inv@Nevt_data)
