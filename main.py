#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Kroto Membrane Form-finding.
#  Copyright © 2015, Thinkshell & Laboratoire Navier, ENPC.
#
#  This file is part of Kroto.
#
#  Kroto is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
#
#  Kroto is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with Kroto. If not, see http://www.gnu.org/licenses/.
#

import meshminimize as mm
import rhinoscriptsyntax as rs

mm.DEBUG = 0
mm.GRAPHIC = 0
mm.SHOW_RESULT = 1
mm.SAVE_RESULTS = 0
mm.MAX_DISP = 0.001
mm.MAX_ITER = 100
mm.MAX_ITER_QS = 10
mm.MAX_DEV_SIGMA = 1
mm.SPEED = 1
mm.METHOD = 'fixed-point'
mm.FIXED_CABLE_ENDS = True


mesh = rs.GetObject("Select mesh", rs.filter.mesh)
cables = rs.GetObjects("Select cables", rs.filter.curve)

q_cables = None
if cables:
    q_cables = [1000 for i in cables]

pressure = 100

objective = mm.minimize_mesh(mesh, cables, fixed=None, qs=None,
                             q_cables=q_cables, reference=None, P=pressure)
