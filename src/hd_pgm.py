#!/usr/bin/env python

from matplotlib import rc
from daft import PGM, Node, Plate
rc("font", family="serif", size=10)
rc("text", usetex=True)


pgm = PGM([7.5, 6.5], origin=[0., 0.2], observed_style='inner')

pgm.add_node(Node('HD',r"$\vec{\mu},\vec{z}$", 1,6))
pgm.add_node(Node('TypeHost',r"Type, Host", 3, 6, scale=1.8))
pgm.add_node(Node('Flux',r"$n(t,\lambda)$", 2, 5))
pgm.add_node(Node('Flux_g',r"$n_g(\lambda)$", 3, 5))
pgm.add_node(Node('Transmission',r"$T(\lambda)$", 6, 6))
pgm.add_node(Node('Spectrum',r"$\hat{\mathit{Spec}}$.", 2, 2, observed=True))
pgm.add_node(Node('Counts',r"$\mathit{ADU}$", 3, 4))
pgm.add_node(Node('Counts_g',r"$\mathit{ADU}_g$", 4, 4))
pgm.add_node(Node('Zeropoints',r"$\hat{Z}$", 6, 5, observed=True))
pgm.add_node(Node('^Counts',r"$\hat{f}$", 3, 3, observed=True))
pgm.add_node(Node('Spars',r"$\hat{z_S},\hat{\theta_S}$", 3, 2, observed=True))
pgm.add_node(Node('^Host',r"$\hat{z_H},\hat{\theta_H}$", 3, 1, observed=True))
pgm.add_node(Node('Galaxies',r"Gals.", 1, 1, observed=True))
pgm.add_node(Node('^z',r"$\hat{z}$", 5, 2, observed=True))
pgm.add_node(Node('^mu',r"$\hat{\mu}$", 5, 1, observed=True))
pgm.add_node(Node('^Type',r"$\hat{\mathit{Type}}$", 5, 3, observed=True))
pgm.add_node(Node('Detected',r"$\hat{D}$", 5, 4, observed=True))

pgm.add_edge("HD","TypeHost")
pgm.add_edge("HD","Flux")
pgm.add_edge("HD","Flux_g")

pgm.add_edge("TypeHost","Flux")
pgm.add_edge("TypeHost","Flux_g")


pgm.add_edge("Flux","Counts")
pgm.add_edge("Flux_g","Counts_g")

pgm.add_edge("Counts","^Counts")
pgm.add_edge("Counts_g","^Counts")
pgm.add_edge("Zeropoints","^Counts")


pgm.add_edge("Transmission","Counts")
pgm.add_edge("Transmission","Counts_g")
pgm.add_edge("Transmission","Zeropoints")

pgm.add_edge("Flux","Spectrum")
pgm.add_edge("Flux_g","Spectrum")
pgm.add_edge("Spectrum","Spars")

pgm.add_edge("Galaxies","^Host")


pgm.add_edge("^Counts","Detected")
pgm.add_edge("^Counts","^Type")
pgm.add_edge("^Counts","^z")
pgm.add_edge("^Counts","^mu")

pgm.add_edge("Spectrum","^Type")

pgm.add_edge("Spars","^z")
pgm.add_edge("Spars","^mu")

pgm.add_edge("^Host","^Type")
pgm.add_edge("^Host","^z")
pgm.add_edge("^Host","^mu")

# pgm.add_edge("Galaxies", "gal_zi")


# Per-Galaxy parameters: third line in the plate


# Observed photometry

# pgm.add_edge("gal_ni","^gal_zi")
# pgm.add_edge("gal_ni","^gal_thetai")
# pgm.add_edge("gal_ni","^gal_ni")


# pgm.add_edge("gal_D","^gal_zi")
# pgm.add_edge("gal_D","^gal_thetai")
# pgm.add_edge("gal_D","^gal_ni")

# pgm.add_edge("^ni","^Type")


# pgm.add_edge("Standard","Z")


#pgm.add_node(Node("t0true", r"$t_0^{\mathrm{true}}$", 7, 1))


# Big Plate: Galaxy
pgm.add_plate(Plate([1.5, 0.5, 4, 6.],
                    label=r"SNe $\i = 1, \cdots, N_{SN}$",
                    shift=-0.1))

# # Big Plate: SNe
# pgm.add_plate(Plate([5.5, 0.5, 4, 6.],
#                     label=r"SNe $i = 1, \cdots, N_{SN}$",
#                     shift=-0.1))


# pgm.add_plate(Plate([2.6, 1.6, 1.8, 3.8],
#                     label=r"Photo-$z$",
#                     shift=-0.1))

# pgm.add_plate(Plate([3.6, 5.6, 2.8, 1.8],
#                     label=r"Host match",
#                     shift=-0.1))

# pgm.add_plate(Plate([5.6, 1.6, 1.8, 3.8],
#                     label=r"LC Fit",
#                     shift=-0.1))

# pgm.add_plate(Plate([5.6, 2.6, 2.8, 0.8],
#                     label=r"Efficiency",
#                     shift=-0.1))


# Render and save.
pgm.render()
pgm.figure.savefig("../results/hdpgm.eps")
