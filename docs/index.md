## sc2heoca

### What is sc2heoca?
A method for querying newly generated organoid single-cell RNA sequencing (scRNA-seq) data to the Human Endoderm Organoid Cell Atlas ([HEOCA](https://cellxgene.cziscience.com/e/6725ee8e-ef5b-4e68-8901-61bd14a1fe73.cxg)). HEOCA is an integrated single-cell transcriptome collection that includes 218 samples from organoids derived from various endoderm tissues, including the lung, pancreas, intestine, liver, biliary system, stomach, and prostate. This current version of HEOCA includes nearly one million cells from a variety of conditions, data sources, and protocols.

!!! note "Citation"
> An integrated transcriptomic cell atlas of human endoderm-derived organoids. Quan Xu, Lennard Halle, Soroor Hediyeh-zadeh, Merel Kuijs, Umut Kilik, Qianhui Yu, Tristan Frum, Lukas Adam, Shrey Parikh, Manuel Gander, Raphael Kfuri-Rubens, Dominik Klein, Zhisong He, Jonas Simon Fleck, Koen Oost, Maurice Kahnwald, Silvia Barbiero, Olga Mitrofanova, Grzegorz Maciag, Kim B. Jensen, Matthias Lutolf, Prisca Liberali, Joep Beumer, Jason R. Spence, Barbara Treutlein, Fabian J. Theis, J. Gray Camp. bioRxiv 2023.11.20.567825; doi: https://doi.org/10.1101/2023.11.20.567825 

### Getting started

* Install sc2heoca on Linux or Mac, see the [Installation](installation.md) page for details.
* Have a look at these simple [examples](examples.md) to get a taste of what is possible.
* Check out the [command-line reference](command-line_reference.md) to get going with your own data.

### Get help

* First, check the [FAQ](faq.md) for common issues.
* The preferred way to get support is through the [Github issues page](https://github.com/devsystemslab/sc2heoca/issues).
* Finally, you can reach us by email to <a href="mailto:qxuchn@gmail.com" target="_blank">Quan Xu</a>.

### Full contents

* [Model description and download](model_description.md)
    - [Overview of sc2heoca model](model_description/#overview_of_sc2heoca)
    - [Prepare reference atlas and model](model_description/#Prepare_reference_atlas_and_model)
* [installation instructions](installation.md)
    - [Creat conda environment](installation/#Creat_conda_environment)
    - [Install sc2heoca from GitHub](installation/#Install_sc2heoca_from_GitHub)
* [Input data and examples](examples.md)
    - [sc2heoca setup](examples/#prepare-code-and-dataset)
    - [Data](examples/#data)
    - [Example of sc2heoca query and annotation](examples/#query)
    - [Example of sc2heoca on/off target influence](examples/#target)
    - [Example of sc2heoca maturation influence](examples/#maturation)
    - [Example of sc2heoca de genes influence](examples/#de_genes)
* [FAQ](faq.md)
* [Acknowledgments](acknowledgments.md)
