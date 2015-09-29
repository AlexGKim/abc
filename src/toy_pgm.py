#!/usr/bin/env python

from matplotlib import rc
from daft import PGM, Node, Plate
rc("font", family="serif", size=12)
rc("text", usetex=True)


pgm = PGM([6.5, 4.5], origin=[0., 0.2], observed_style='inner')

#pgm.add_node(Node('G',r"$G$", 3,1))
#pgm.add_node(Node('Coords',r"${RA}_i$/${Dec}_i$", 2,1,fixed=True))
pgm.add_node(Node('G_i',r"$z_i, z_{Ni}$", 2,2, scale=1.2))

pgm.add_node(Node('mu',r"$\Omega_M$, $w$", 3,4, scale=1.4))
pgm.add_node(Node('rate',r"$\theta_{r1}$, $\theta_{r2}$", 5,4,scale=1.4))
pgm.add_node(Node('HD',r"$\mu_i$", 3,2,fixed=True,offset=(0,-22)))
pgm.add_node(Node('theta_T',r"\center{$\alpha_{Ia}$, $\alpha_{non-Ia}$ \newline $\sigma_{Ia}$, $\sigma_{non-Ia}$}", 4,4,scale=1.65,aspect=1.5))
#pgm.add_node(Node('theta_Ti',r"$\theta_{Ti}^{Ia}$, $\theta_{Ti}^{non-Ia}$", 2,6,scale=1.5,aspect=1.5))

#pgm.add_node(Node('Host',r"$\theta_{gi}$", 2, 3, fixed=True,offset=(-10,-5)))
#pgm.add_node(Node('z',r"$z_i$", 2, 2, fixed=True,offset=(-10,-5)))

pgm.add_node(Node('Type',r"$T_i$", 5,3))
pgm.add_node(Node('Luminosity',r"$L_i(t,\lambda)$", 4,3, scale=1.4))
#pgm.add_node(Node('Flux',r"$n_i(t,\lambda)$", 5, 5, scale=1.2,fixed=True,offset=(15,0)))
#pgm.add_node(Node('Flux_g',r"$n_{gi}(\lambda)$", 5, 2,fixed=True,offset=(0,-20)))
#pgm.add_node(Node('Transmission',r"$\phi(\lambda)$", 7, 7))
#pgm.add_node(Node('Counts',r"$\overline{f}_i$", 8, 5,scale=1.2,fixed=True,offset=(15,0)))
#pgm.add_node(Node('Counts_g',r"$\overline{f}_{gi}$", 8, 4,scale=1.2,fixed=True,offset=(10,-25)))
#pgm.add_node(Node('Zeropoints',r"${Z}$", 9, 7, observed=True))
pgm.add_node(Node('^Counts',r"${f_i}$", 4, 2, observed=True))

pgm.add_node(Node('Spars',r"${T}_{Si}$, ${z}_{Si}$", 4,1, scale=1.4,aspect=1.2, observed=True))
#pgm.add_node(Node('^Host',r"${z}_{Hi},{\theta}_{Hi}$", 7, 2, scale=1.5,observed=True))

#pgm.add_node(Node('Detected',r"Detected$_i$", 8, 3, fixed=True,offset=(-10,-20)))
#pgm.add_node(Node('^Type',r"$\tau_i$", 9, 2, fixed=True,offset=(10,-10)))

pgm.add_node(Node('Gals',r"$\theta_{G1}$, $\theta_{G2}$", 1,1, scale=1.4, aspect=1.2, observed=True))


pgm.add_edge("G_i","Gals")

pgm.add_edge("rate","Type")


pgm.add_edge("mu","HD")
#pgm.add_edge("^Counts","^Type")
#pgm.add_edge("Spars","^Type")
#pgm.add_edge("Host","Luminosity")


pgm.add_edge("G_i","HD")


# pgm.add_edge("Coords","G_i")
#pgm.add_edge("G_i","Host")

#pgm.add_edge("G_i","z")
pgm.add_edge("HD","^Counts")
#pgm.add_edge("G_i","^Counts")
#pgm.add_edge("G_i","Flux_g")
#pgm.add_edge("HD","Flux_g")
#pgm.add_edge("z","Flux_g")
#pgm.add_edge("HD","^Host")


#pgm.add_edge("Host","Type")

pgm.add_edge("Type","Luminosity")
pgm.add_edge("Luminosity","^Counts")


pgm.add_edge("theta_T","Luminosity")
#pgm.add_edge("Galaxies","G")

#pgm.add_edge("theta_Ti","Luminosity")


#pgm.add_edge("Flux","Counts")
#pgm.add_edge("Flux_g","Counts_g")

#pgm.add_edge("Counts","^Counts")
#pgm.add_edge("Counts_g","^Counts")


#pgm.add_edge("Transmission","Counts")
#pgm.add_edge("Transmission","Counts_g")
#pgm.add_edge("Transmission","Zeropoints")

pgm.add_edge("Type","Spars")
pgm.add_edge("G_i","Spars")

#pgm.add_edge("Counts_g","^Counts")


#pgm.add_edge("Detected","^Counts")




# Big Plate: Galaxy
pgm.add_plate(Plate([1.5, 0.5, 4.2, 3],
                    label=r"SNe $i = 1, \cdots, N_{SN}$",
                    shift=-0.1))


# Render and save.
pgm.render()
# pgm.figure.text(0.2,0.98,r'\underline{UNIVERSE}',size='large')
# pgm.figure.text(0.45,0.98,r'\underline{OBSERVATORY}',size='large')
# pgm.figure.text(0.72,0.98,r'\underline{DATA}',size='large')

pgm.figure.savefig("../results/toy_pgm.pdf")
