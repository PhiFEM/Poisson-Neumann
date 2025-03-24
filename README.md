# $\varphi$-FEM for Poisson Neumann

Repository containing the codes used in M. Duprez, V. Lleras, and A. Lozinski. "A new $\phi$-FEM approach for problems with natural boundary conditions", *Numer. Methods Partial Differ. Eq. 39* (2023), 281-303. [https://doi.org/10.1002/num.22878](https://onlinelibrary.wiley.com/doi/10.1002/num.22878)

## This repository is for reproducibility purposes only

It is "frozen in time" and not maintained.
To use our latest $\varphi$-FEM code please refer to the [phiFEM repository](https://github.com/PhiFEM/Poisson-Dirichlet-fenicsx).

## Generate the results

The results can be generated using two container images:
- the (legacy) FEniCS image is based on [ghcr.io/scientificcomputing/fenics-gmsh:2024-05-30](https://github.com/scientificcomputing/packages/pkgs/container/fenics-gmsh), with the additional python libraries [`pyamg`](https://github.com/pyamg/pyamg) and [`multiphenics`](https://github.com/multiphenics/multiphenics/tree/7b23c85c070a092775666c7dad84c8d6471c0b0c) (legacy),
- the [FreeFem++](https://doc.freefem.org/introduction/index.html) image is [https://hub.docker.com/r/freefem/freefem](https://hub.docker.com/r/freefem/freefem).

### Prerequisites

- [Git](https://git-scm.com/),
- [Docker](https://www.docker.com/)/[podman](https://podman.io/).

### Install the image and launch the container

1) Clone this repository in a dedicated directory:
   
   ```bash
   mkdir poisson-neumann-phifem/
   git clone https://github.com/PhiFEM/Poisson-Neumann.git poisson-neumann-phifem
   ```

2) Download the images from the docker.io registry, in the main directory:
   
   ```bash
   export CONTAINER_ENGINE=docker
   cd poisson-neumann-phifem
   sudo -E bash pull-fenics-image.sh
   sudo -E bash pull-freefem-image.sh
   ```

3) Launch the FEniCS container:

   ```bash
   sudo -E bash run-fenics-image.sh
   ```

   or the FreeFem++ container:

   ```bash
   sudo -E bash run-freefem-image.sh
   ```

### Example of usage

#### FEniCS

From the main directory `poisson-neumann-phifem`, launch the 1st test case:

```bash
python3 phiFEM_test_case1_and_2_block.py
```

#### FreeFem

From the main directory `poisson-neumann-phifem`, launch the 1st test case:

```bash
FreeFem++ testcas1.edp
```

## Issues and support

Please use the issue tracker to report any issues.

## Authors (alphabetical)

[Michel Duprez](https://michelduprez.fr/), Inria Nancy Grand-Est  
[Vanessa Lleras](https://vanessalleras.wixsite.com/lleras), Université de Montpellier  
[Alexei Lozinski](https://orcid.org/0000-0003-0745-0365), Université de Franche-Comté  