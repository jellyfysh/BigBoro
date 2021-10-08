[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)

# BigBoro

The BigBoro software package accompanies the work "Sparse hard-disk packings and local Markov chains" (see 
[[Hoellmer2021]](https://arxiv.org/abs/2109.13343) in the [References.bib](References.bib) file), where we propose 
locally stable sparse hard-disk packings, as introduced by Böröczky, as a model for the analysis and benchmarking of 
Markov-chain Monte Carlo (MCMC) algorithms. Here, we focus on such Böröczky packings in a square box with periodic 
boundary conditions. The work studies how local MCMC algorithms, namely the Metropolis algorithm and several versions of 
event-chain Monte Carlo (ECMC), escape from configurations that are obtained by slightly reducing all disk radii by a 
relaxation parameter (epsilon-relaxed Böröczky configurations). A scaling analysis is confirmed by simulation results. 
We obtain two classes of ECMC, one in which the escape time varies algebraically with the relaxation parameter (as for 
the local Metropolis algorithm) and another in which the escape time scales as the logarithm of the relaxation 
parameter. 

This BigBoro software package was used to construct Böröczky packings (see, for example, Fig. 1 and Table 1 of 
[[Hoellmer2021]](https://arxiv.org/abs/2109.13343)). It also computed collective escape modes from such locally stable 
packings, where all disks are infinitesimally displaced simultaneously (see Fig. 2). Finally, its arbitrary-precision 
ECMC implementation confirmed the different scalings of escape times from epsilon-relaxed Böröczky configurations (see 
Figs 4, 5, and 6). To be precise, the BigBoro software package consists of three parts:

1. The arbitrary-precision Python script [`construct_packing.py`](construct_packing/construct_packing.py)
constructs finite-N Böröczky packings of hard disks in a periodic square box (see 
[`construct_packing/`](construct_packing) directory). 
2. The Python script [`collective_escape_modes.py`](collective_escape_modes/collective_escape_modes.py) computes 
collective infinitesimal displacements of hard disks in a packing that result in an escape (see
[`collective_escape_modes/`](collective_escape_modes) directory). 
3. The arbitrary-precision Go application [`go-hard-disks`](go-hard-disks) performs
hard-disk ECMC simulations that may start from epsilon-relaxed Böröczky configurations derived from Böröczky packings
(see [`go-hard-disks/`](go-hard-disks) directory).

## Installing and using

In order to use (parts of) the BigBoro software package, it is necessary to clone this repository. We refer to the 
`README.md` files within the directories of each part for further details on the installation and usage. 

## Contributing

If you find a bug, please raise an Issue here on GitHub to let us know.

As an open-source project, the BigBoro software package solicits contributions from the community. We would be happy to 
receive your fixes, extensions for the BigBoro software packages, or solutions to open issues, and are looking forward 
to your pull requests. After a successful review of your code, we will merge your code into the
master branch and add your name to the [AUTHORS.md](AUTHORS.md) file. In addition, your name will appear in the commit
history of the repository.

For successful pull requests make sure you follow the following points:

- Follow the coding style of the existing code.
- Test your code.
- Read the [code of conduct](CODE_OF_CONDUCT.md).
- At the top of any new file, include the license notice that you can find in the existing files.

We will be happy to assist contributors with their first few pull requests.

Please note that this project is released with the Contributor Covenant [code of conduct](CODE_OF_CONDUCT.md). By 
participating in this project you agree to abide by its terms. Report unacceptable behavior to 
[werner.krauth@ens.fr](mailto:werner.krauth@ens.fr).

## Versioning

Versioning of the BigBoro software project adopts two-to-four-field version numbers defined as 
Milestone.Feature.AddOn.Patch. The current version 1.0 represents the first development milestone which reproduces 
published data in [[Hoellmer2021]](https://arxiv.org/abs/2109.13343). Patches and bugfixes of this version will be given 
number 1.0.0.1, 1.0.0.2 etc. New extensions are expected to lead to versions 1.0.1, 1.0.2 etc. In the development of the 
BigBoro software package, two-field versions (2.0, 3.0, etc.) may introduce incompatible code, while three- and 
four-field version numbers are intended to be backward compatible.

## Authors

Check the [AUTHORS.md](AUTHORS.md) file to see who participated in this project.

## License

This project is licensed under the GNU General Public License, version 3 (see the [LICENSE](LICENSE) file).

## Contact

If you have questions regarding the BigBoro software package, just raise an issue here on GitHub or contact us via mail
(see the [AUTHORS.md](AUTHORS.md) file). We are happy to help you!

## Citation

If you use (parts of the) BigBoro software package in published work, please cite the following reference (see
[[Hoellmer2021]](https://arxiv.org/abs/2109.13343) in [References.bib](References.bib)):

Philipp Höllmer, Nicolas Noirault, Botao Li, A. C. Maggs, Werner Krauth,\
Sparse hard-disk packings and local Markov chains,\
arXiv e-prints: 2109.13343 (2021), https://arxiv.org/abs/2109.13343.
