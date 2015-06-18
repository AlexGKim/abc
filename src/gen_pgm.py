#!/usr/bin/env python

from matplotlib import rc
from daft import PGM, Node, Plate
rc("font", family="serif", size=10)
rc("text", usetex=True)


pgm = PGM([10.5, 7.5], origin=[0., 0.2], observed_style='inner')

pgm.add_node(Node("Standard", r"Stands.", 5, 1, scale=1.5))
pgm.add_node(Node("^thetai", r"$\hat{\theta}_{S i}$", 8, 1, observed=True))
pgm.add_node(Node("^zi", r"$\hat{z}_{Si}$", 9, 1, observed=True))

pgm.add_node(Node("^ni", r"$\hat{ADU}_{i}$", 6, 2, observed=True,scale=1.5))
pgm.add_node(Node("gal_D", r"$\hat{D}$", 1, 2, observed=True))
pgm.add_node(Node("^gal_ni", r"$\hat{ADU}_{\alpha}$", 4, 2, observed=True,scale=1.5))
pgm.add_node(Node("Z", r"$\hat{Z}$", 5, 2, observed=True))


pgm.add_node(Node("^gal_thetai", r"$\hat{\theta}_{S \alpha}$", 2, 3, observed=True))
pgm.add_node(Node("Optics", r"Trans.", 5, 3, scale=1.5))
pgm.add_node(Node("gal_ni", r"$ADU_{\alpha}$", 3, 3,scale=1.5))
pgm.add_node(Node("ni", r"$ADU_{i}$", 6, 3,scale=1.5))

pgm.add_node(Node("^gal_zi", r"$\hat{z}_{S \alpha}$", 1, 4, observed=True))
pgm.add_node(Node("D_L", r"$D_L$", 5, 4))
pgm.add_node(Node("D", r"$\hat{D}$", 8, 4, observed=True))

pgm.add_node(Node("gal_zi", r"$z_{\alpha}$", 4, 5))
pgm.add_node(Node("zi", r"$z_{i}$", 6, 5))
pgm.add_node(Node("pi", r"$p^{SN}_{i}$", 7, 5))

pgm.add_node(Node("gal_Li", r"$L_{\alpha},p^{gal}_\alpha$", 2, 6,scale=1.5))
pgm.add_node(Node("coord", r"$RA/Dec_i$", 6, 6, scale=1.5))
pgm.add_node(Node("^Type", r"$\hat{Type}_{i}$", 8, 6, observed=True,scale=1.5))
pgm.add_node(Node("SN pdf", r"SN $pdf_{\alpha}$", 4, 6,scale=1.5))
pgm.add_node(Node("^Host", r"$\hat{Host}_{i}$", 7, 6, observed=True,scale=1.5))

pgm.add_node(Node("Galaxies", r"Gal. $pdf$",2, 7, scale=1.5))
pgm.add_node(Node("gal_coord", r"$\{RA/Dec\}$", 4, 7, scale=1.7))

pgm.add_edge("Galaxies", "gal_zi")
pgm.add_edge("Galaxies", "gal_Li")
pgm.add_edge("Galaxies", "gal_coord")

pgm.add_edge("gal_coord","^Host")

pgm.add_edge("^Host","^Type")

pgm.add_edge("SN pdf","pi")
pgm.add_edge("SN pdf", "coord")

pgm.add_edge("gal_Li","gal_ni")
pgm.add_edge("gal_Li", "SN pdf")

pgm.add_edge("coord","^Host")

pgm.add_edge("gal_zi","D_L")
pgm.add_edge("gal_zi","zi")
pgm.add_edge("gal_zi","gal_ni")

pgm.add_edge("zi","D_L")
pgm.add_edge("zi","ni")

pgm.add_edge("pi","ni")

pgm.add_edge("D_L","gal_ni")
pgm.add_edge("D_L","ni")

pgm.add_edge("gal_ni", "^ni")

pgm.add_edge("Optics","gal_ni")
pgm.add_edge("Optics","ni")
pgm.add_edge("Optics","Z")

# pgm.add_edge("ni","D")
pgm.add_edge("ni","^ni")
pgm.add_edge("ni","^zi")
pgm.add_edge("ni","^thetai")

pgm.add_edge("D","^thetai")
pgm.add_edge("D","^Type")
pgm.add_edge("D","^zi")

pgm.add_edge("^gal_ni", "D")
pgm.add_edge("^gal_ni", "gal_D")

pgm.add_edge("^ni","D")

# Per-Galaxy parameters: third line in the plate


# Observed photometry

pgm.add_edge("gal_ni","^gal_zi")
pgm.add_edge("gal_ni","^gal_thetai")
pgm.add_edge("gal_ni","^gal_ni")


pgm.add_edge("gal_D","^gal_zi")
pgm.add_edge("gal_D","^gal_thetai")
pgm.add_edge("gal_D","^gal_ni")

pgm.add_edge("^ni","^Type")


pgm.add_edge("Standard","Z")


#pgm.add_node(Node("t0true", r"$t_0^{\mathrm{true}}$", 7, 1))


# Big Plate: Galaxy
pgm.add_plate(Plate([0.5, 0.5, 4, 6.],
                    label=r"Galaxies $\alpha = 1, \cdots, N_{gal}$",
                    shift=-0.1))

# Big Plate: SNe
pgm.add_plate(Plate([5.5, 0.5, 4, 6.],
                    label=r"SNe $i = 1, \cdots, N_{SN}$",
                    shift=-0.1))


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
pgm.figure.savefig("../results/snpgm.eps")
