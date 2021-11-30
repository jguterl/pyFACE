#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 11:53:29 2021

@author: guterlj
"""

import pyface
state=True
pyface.solver.allocate()
pyface.solver.initialize_wrapper()
state=True
pyface.solver.do_step(state)
pyface.input.verbose_compute=True
pyface.solver.update_dt()
state=True
for i in range(100):
    pyface.solver.do_step(state)
    pyface.solver.update_dt()
   