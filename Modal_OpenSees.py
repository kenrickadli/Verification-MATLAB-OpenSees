# =========================================================================
# 3D PORTAL FRAME - Modal Analysis (Eigenvalue Analysis)
# =========================================================================

import os
import numpy as np
import openseespy.opensees as ops
import opsvis as opsv
import matplotlib.pyplot as plt
import pandas as pd

ops.wipe()
ops.model('BasicBuilder', '-ndm', 3, '-ndf', 6) 
ops.defaultUnits("-force", "N", "-length", "m", "-time", "sec") 
g = 9.80665  # m/s²

# =========================================================================
# STRUCTURE GEOMETRY
# =========================================================================
StructureArray = [
    7.2,    # Lx - bay width (m)
    7.2,    # Ly - bay depth (m)
    3.6,    # Lz - story height (m)
    1,      # nx - number of bays in X
    1,      # ny - number of bays in Y
    1,      # nz - number of stories
]

Lx, Ly, Lz, nx, ny, nz = StructureArray
ndt = (nx + 1) * (ny + 1) * (nz + 1)  # Total nodes
nds = (nx + 1) * (ny + 1)             # Nodes per floor

# =========================================================================
# LOADING
# =========================================================================
LoadingArray = [
    2.0,    # SIDL (kN/m²)
    2.4,    # Live (kN/m²)
    180,    # h_slab (mm)
    2400,   # Concrete density (kg/m³)
]

SIDL, Live, h_slab, rho_c = LoadingArray

floor_area = (Lx * nx) * (Ly * ny) # (m²)
m_Dead = rho_c * (h_slab * 1e-3) * floor_area # (kg)
m_SIDL = (SIDL * 1e3) * floor_area / g  # (kg)
m_Live = (Live * 1e3) * floor_area / g  # (kg)
mass_source = (m_Dead + m_SIDL + 0.5 * m_Live)
mz = mass_source * ((Lx * nx)**2 + (Ly * ny)**2) / 12  # (kg⋅m²)
 
print(f"Total mass per floor: {mass_source:.0f} kg")

# =========================================================================
# MATERIAL PROPERTIES
# =========================================================================
MaterialArray = [
    345,    # Fy (MPa)
    1.1,    # Ry
    200000, # Es (MPa)
    0.3     # v
]

Fy, Ry, Es, v = MaterialArray
Es = Es * 1e6
G = Es / (2 * (1 + v))
J = 1e-10

print(Es,G)

# =========================================================================
# BEAM SECTION PROPERTIES
# =========================================================================
BeamArray = [
    312.4,  # d_beam
    304.8,  # bf_beam
    17,     # tf_beam
    10.9,   # tw_beam
]

d_beam, bf_beam, tf_beam, tw_beam = BeamArray

A_beam = (bf_beam * tf_beam * 2 + (d_beam - 2 * tf_beam) * tw_beam) * 1e-6 # (m²)
Ix_beam = ((bf_beam * d_beam**3 / 12) - ((bf_beam - tw_beam) * (d_beam - 2 * tf_beam)**3 / 12)) * 1e-12 # (m⁴)
Iy_beam = ((d_beam - 2 * tf_beam) * tw_beam**3 / 12 + 2 * (tf_beam * bf_beam**3 / 12)) * 1e-12 # (m⁴)
ry_beam = ((Iy_beam / A_beam) ** 0.5) # (m)
Sx_beam = Ix_beam / (d_beam / 2) # (m³)
Sy_beam = Iy_beam / (bf_beam / 2) # (m³)
Zx_beam = (bf_beam * d_beam**2 / 4) - ((bf_beam - tw_beam) * (d_beam - 2 * tf_beam)**2 / 4) # (m³)
Zy_beam = 4 * tf_beam * (bf_beam**2 / 8) + (d_beam - 2 * tf_beam) * (tw_beam**2 / 4) # (m³)
L_beam = Lx

print("\n=== BEAM PROPERTIES ===")
print(f"Area: {A_beam:.5e} m²")
print(f"Ix: {Ix_beam:.5e} m⁴, Iy: {Iy_beam:.5e} m⁴")
print(f"ry: {ry_beam:.5e} m")
print(f"Sx: {Sx_beam:.5e} m³, Sy: {Sy_beam:.5e} m³")
print(f"Zx: {Zx_beam:.5e} m³, Zy: {Zy_beam:.5e} m³")

# =========================================================================
# DEFINE NODES
# =========================================================================
n = 1
node_coords = {}
for k in range(nz + 1):
    for j in range(ny + 1): 
        for i in range(nx + 1): 
            x_coord = i * Lx
            y_coord = j * Ly  
            z_coord = k * Lz
            ops.node(n, x_coord, y_coord, z_coord)
            node_coords[n] = (x_coord, y_coord, z_coord, i, j, k)
            n += 1 

# =========================================================================
# BASE RESTRAINTS
# =========================================================================
for i in range(1, nds + 1): 
    ops.fix(i, 1, 1, 1, 1, 1, 1) # Fixed base

# =========================================================================
# RIGID DIAPHRAGM
# =========================================================================
rDtag = []
for k in range(1, nz + 1):
    master_node = ndt + k
    rDtag.append(master_node)
    ops.node(master_node, Lx * nx / 2, Ly * ny / 2, k * Lz)

for k in range(1, nz + 1): 
    floor_nodes = [int(n) for n in np.arange(k*nds + 1, (k + 1)*nds + 1)] 
    ops.rigidDiaphragm(3, rDtag[k-1], *floor_nodes)

for i in range(len(rDtag)): 
    ops.fix(rDtag[i], 0, 0, 1, 1, 1, 0)

for i in range(len(rDtag)):
    ops.mass(rDtag[i], mass_source, mass_source, 0, 0, 0, mz)

# =========================================================================
# SECTIONS AND TRANSFORMATIONS
# =========================================================================
ops.section('Elastic', 1, Es, A_beam, Ix_beam, Iy_beam, G, J) # Beam
ops.section('Elastic', 2, Es, A_beam, Ix_beam, Iy_beam, G, J) # Column

ops.geomTransf('Linear', 1, 0, 1, 0)  # Beam X-direction (local y = global z)
ops.geomTransf('Linear', 2, 1, 0, 0)  # Beam Y-direction (local x = global y)  
ops.geomTransf('Linear', 3, 1, 0, 0)  # Column (local x = global y)

# =========================================================================
# DEFINE ELEMENTS
# =========================================================================
eletag = 1
beam_count_x = 0
beam_count_y = 0
column_count = 0

# Beam
nn = nds + 1
for k in range(1, nz + 1):
    for j in range(ny + 1):
        for i in range(nx + 1):
            if i < nx:
                node1 = nn
                node2 = nn + 1
                ops.element('elasticBeamColumn', eletag, node1, node2, 1, 1)
                beam_count_x += 1
                eletag += 1 
            nn += 1

# Beam Y Direction
nn = nds + 1
for k in range(1, nz + 1):
    for j in range(ny + 1):
        for i in range(nx + 1):
            if j < ny:
                node1 = nn
                node2 = nn + (nx + 1)
                ops.element('elasticBeamColumn', eletag, node1, node2, 1, 2)
                beam_count_y += 1
                eletag += 1 
            nn += 1      

# Column
nn = 1
for k in range(nz):
    for j in range(ny + 1):
        for i in range(nx + 1):
            node1 = nn
            node2 = nn + nds
            ops.element('elasticBeamColumn', eletag, node1, node2, 2, 3)
            column_count += 1
            nn += 1 
            eletag += 1

# =========================================================================
# LOAD PATTERN
# =========================================================================
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)

# =========================================================================
# EIGENVALUE ANALYSIS
# =========================================================================
n_modes = 3
eigenValues = ops.eigen('-fullGenLapack', n_modes)  # Use different solver

# =========================================================================
# RESULTS
# =========================================================================
periods = []
print("\n=== EIGENVALUE ANALYSIS ===")

for i, lamb in enumerate(eigenValues):
    omega = lamb**0.5
    period = 2 * np.pi / omega
    periods.append(period)
    print(f"Mode {i+1}: λ = {lamb:.6f}, ω = {omega:.6f}, T = {period:.6f} sec")