#!/usr/bin/env python

from matplotlib import rc
from daft import PGM, Node, Plate
rc("font", family="serif", size=10)
rc("text", usetex=True)


pgm = PGM([4.5, 4.5], origin=[0., 0.2], observed_style='inner')



pgm.add_node(Node('Pars',r"$\theta_G$", 1, 3))
pgm.add_node(Node('mu',r"$\mu_G$", 1, 1))
pgm.add_node(Node('SNpars',r"$\theta^{dist}$", 1, 4))
pgm.add_node(Node('^z',r"$\hat{z_i}$", 3,1))

pgm.add_node(Node('Spars',r"$\hat{T}_{Si}$,$\hat{z_{Si}},\hat{\theta}_{Si}$", 3, 2, scale=1.8,observed=True))
pgm.add_node(Node('^Counts',r"$\hat{f_i}$", 2, 2, observed=True))
pgm.add_node(Node('^Host',r"$\hat{z_{Hi}},\hat{\theta}_{Hi}$", 3, 3, observed=True,scale=1.5))

pgm.add_node(Node('^mu',r"$\hat{\mu_i}$", 2, 1))

pgm.add_node(Node('SNpars_i',r"$\theta^{true}_i$", 2, 3))




pgm.add_node(Node('^Type',r"$\hat{T_i}$", 2, 4))




pgm.add_edge("^Type","SNpars_i")


pgm.add_edge("SNpars","SNpars_i")

pgm.add_edge("SNpars_i","^Counts")
pgm.add_edge("SNpars_i","Spars")
pgm.add_edge("SNpars_i","^Host")
pgm.add_edge("Pars","^Counts")
pgm.add_edge("^mu","^Counts")

pgm.add_edge("mu","^mu")
pgm.add_edge("^z","Spars")

#pgm.add_node(Node("t0true", r"$t_0^{\mathrm{true}}$", 7, 1))


# Big Plate: Galaxy
# pgm.add_plate(Plate([0.5, 0.5, 4, 6.],
#                     label=r"Galaxies $\alpha = 1, \cdots, N_{gal}$",
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
pgm.figure.savefig("../results/snpgm.eps")
