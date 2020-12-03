metabarcoding with dada2: environment installation
================

# update VM configuration

``` bash
sudo apt-get update -y 
sudo apt-get install -y libbz2-dev
sudo apt-get install -y liblzma-dev
```

    ## sudo: unable to resolve host d8eb10e9e619: Name or service not known
    ## Hit:1 http://archive.ubuntu.com/ubuntu focal InRelease
    ## Hit:2 http://security.ubuntu.com/ubuntu focal-security InRelease
    ## Hit:3 http://archive.ubuntu.com/ubuntu focal-updates InRelease
    ## Hit:4 http://archive.ubuntu.com/ubuntu focal-backports InRelease
    ## Reading package lists...
    ## sudo: unable to resolve host d8eb10e9e619: Name or service not known
    ## Reading package lists...
    ## Building dependency tree...
    ## Reading state information...
    ## libbz2-dev is already the newest version (1.0.8-2).
    ## 0 upgraded, 0 newly installed, 0 to remove and 40 not upgraded.
    ## sudo: unable to resolve host d8eb10e9e619: Name or service not known
    ## Reading package lists...
    ## Building dependency tree...
    ## Reading state information...
    ## liblzma-dev is already the newest version (5.2.4-1ubuntu1).
    ## 0 upgraded, 0 newly installed, 0 to remove and 40 not upgraded.

\#package install

Following instruction on
<https://benjjneb.github.io/dada2/dada-installation.html> \#précède a
installé une première version de dada2 et à obtenir la version actuelle
de dada2 en mettant à jour vers Bioconductor 3.11

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = '3.11')
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

``` r
BiocManager::install("dada2", version = "3.11")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'dada2'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

\#phangorn: package pour la reconstruction et l’analyse phylogénétique.
phangorn offre la possibilité de reconstruire des phylogénies avec des
méthodes basées sur la distance, la parcimonie maximale ou le maximum
de vraisemblance (ML) et d’effectuer une conjugaison Hadamard..

``` r
BiocManager::install("phangorn")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'phangorn'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

\#DECIPHER :ensemble d’outils logiciels qui peut être utilisé pour
déchiffrer et gérer efficacement des séquences biologiques

``` r
BiocManager::install("DECIPHER")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'DECIPHER'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

\#Fournit un certain nombre de fonctions au niveau de l’utilisateur pour
travailler avec des Grid graphiques, notamment pour organiser plusieurs
tracés basés sur une grille sur une page et dessiner des tableaux

``` r
install.packages("gridExtra")
```

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

\#Installer quelques packages disponibles pour ces analyses
complémentaires

``` r
.cran_packages <- c( "shiny","miniUI", "caret", "pls", "e1071", "ggplot2", "randomForest", "dplyr", "ggrepel", "nlme", "devtools",
                  "reshape2", "PMA", "structSSI", "ade4",
                  "ggnetwork", "intergraph", "scales")
.github_packages <- c("jfukuyama/phyloseqGraphTest")
.bioc_packages <- c("genefilter", "impute")
```

\#les packages CRAN ne font pas la distinction entre les branches devel
et release. afin d’installer à partir de miroirs CRAN préconfigurés
plutôt que de saisir une URL et de s’approvisionner sur le Web.

``` r
install.packages(.cran_packages)
```

    ## Installing packages into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

    ## Warning: package 'structSSI' is not available for this version of R
    ## 
    ## A version of this package for your version of R might be available elsewhere,
    ## see the ideas at
    ## https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages

``` r
devtools::install_github(.github_packages)
```

    ## Skipping install of 'phyloseqGraphTest' from a github remote, the SHA1 (3fb6c274) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
BiocManager::install(.bioc_packages)
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'genefilter', 'impute'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

``` r
library(vegan)
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.5-7

\#installer ce package qui fait une estimation modérée du changement de
pli et de la dispersion pour les données RNA-seq.

``` r
BiocManager::install("DESeq2")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'DESeq2'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

\#\#Ce package met en œuvre deux techniques récemment développées dans
le domaine de l’inférence sélective et simultanée (SSI).

``` r
library(devtools)
```

    ## Loading required package: usethis

    ## 
    ## Attaching package: 'devtools'

    ## The following object is masked from 'package:permute':
    ## 
    ##     check

``` r
install_local("./structSSI_1.1.1.tar.gz")
```

    ## Skipping 1 packages not available: multtest

    ##      checking for file ‘/tmp/Rtmp2HMOq7/remotes2f291778755a/structSSI/DESCRIPTION’ ...  ✓  checking for file ‘/tmp/Rtmp2HMOq7/remotes2f291778755a/structSSI/DESCRIPTION’ (338ms)
    ##   ─  preparing ‘structSSI’:
    ##      checking DESCRIPTION meta-information ...  ✓  checking DESCRIPTION meta-information
    ##   ─  checking for LF line-endings in source and make files and shell scripts
    ##   ─  checking for empty or unneeded directories
    ## ─  looking to see if a ‘data/datalist’ file should be added
    ##   ─  building ‘structSSI_1.1.1.tar.gz’
    ##      
    ## 

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

\#CRAN est un réseau de serveurs FTP et Web à travers le monde qui
stockent des versions identiques et à jour du code et de la
documentation pour R.

``` bash
wget https://cran.r-project.org/src/contrib/Archive/structSSI/structSSI_1.1.1.tar.gz
```

    ## --2020-12-03 23:23:06--  https://cran.r-project.org/src/contrib/Archive/structSSI/structSSI_1.1.1.tar.gz
    ## Resolving cran.r-project.org (cran.r-project.org)... 137.208.57.37
    ## Connecting to cran.r-project.org (cran.r-project.org)|137.208.57.37|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 25591 (25K) [application/x-gzip]
    ## Saving to: ‘structSSI_1.1.1.tar.gz.2’
    ## 
    ##      0K .......... .......... ....                            100% 1.05M=0.02s
    ## 
    ## 2020-12-03 23:23:06 (1.05 MB/s) - ‘structSSI_1.1.1.tar.gz.2’ saved [25591/25591]

\#usethis(=use\_release\_issue) pour créer une liste de contrôle des
tâches d’étude

``` r
install.packages("usethis")
```

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

``` r
library(usethis)
use_git_config(user.name = "Masmoudi-FZ", user.email = "fatizohra.masmoudi@gmail.com")
```

``` r
save.image (file = "my_work_space.RData")
```
