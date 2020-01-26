<!--
author:   Alexander Botzki
email:    Alexander.Botzki@vib.be
version:  0.1.1
language: en
narrator: US English Female

comment:  slides Docker Introduction

logo: img/Logo.png

link:     https://cdnjs.cloudflare.com/ajax/libs/animate.css/3.7.2/animate.min.css
link:     https://raw.githubusercontent.com/vibbits/material-liascript/master/img/org.css

debug: true

-->


# Docker Introduction

25.02.2020 UGent

## what options are available

- for ad hoc workflows and tools: docker, singularity, make
- for lightweight application: jupyter notebooks
- for ad hoc pipelines: shell scripts in github
- for (routine) pipelines: snakemake, Galaxy, nextflow (with tools from bioconda)     
- for software tools: bioconda/biocontainers
- ...


## So why does everyone love containers and Docker?

                 {{0-1}}
************************************************

- VM hypervisors are fat in terms of system requirements
- small, neat capsule containing your application
- enables CI, CD
- containers gives you instant application portability and easy to deploy in a cloud
- make applications and workloads more portable in an effective, standardized, and repeatable way [^1](https://www.zdnet.com/article/what-is-docker-and-why-is-it-so-darn-popular/)

************************************************

                 {{1}}
************************************************

![Containers VMs](https://zdnet2.cbsistatic.com/hub/i/r/2017/05/08/af178c5a-64dd-4900-8447-3abd739757e3/resize/770xauto/78abd09a8d41c182a28118ac0465c914/docker-vm-container.png)<!-- width="100%"" -->

************************************************


### virtualisation

pros and cons

- ++ very similar to a full OS
- ++ high OS diversity
- -- need of more space and resources
- -- slower than containers
- -- not as good automation

### containers

pros and cons

- ++ faster
- ++ no need for full OS
- ++ easy solutions for distribution of recipes. high portability
- ++ easy to automate
- -- still OS dependant solutions
- -- not real OS in some cases

## Docker

                 {{0-1}}
************************************************

![Docker](https://msdnshared.blob.core.windows.net/media/2017/10/docker.png)

- platform for developing, shipping, and running applications
- infrastructure as application/code
- Open Container Initiative
- Docker community edition

************************************************

                 {{1}}
************************************************

![Docker components](https://docs.docker.com/engine/images/architecture.svg)

************************************************

## Docker image

- read-only templates
- containers are run from them
- images are not run

### Docker image - building

- can be built from existing images
  - ubuntu, alpine
- base images can be created with tools such as Debootstrap
- any modification from base image is a new layer ( tip: use && )
- images have several layers

### Docker image - instructions

- Recipe: Dockerfile
- Instructions
- FROM
- ADD, COPY
- RUN
- ENV PATH, ARG
- USER, WORKDIR, LABEL
- VOLUME, EXPOSE
- CMD, (ENTRYPOINT)

[Reference](https://docs.docker.com/engine/reference/builder/)

### One tool, one image

- start from packages e.g. [pip/PyPI](https://pypi.org/), [CPAN](https://www.cpan.org/), or [CRAN](https://www.cran.org/)
- use versions for tools and images
- use ENV PATH instead of ENTRYPOINT
- reduce size as much as possible
- keep data outside the image/container
- check the license
- make your container discoverable e.g. biocontainers, quay.io, docker hub

## Example from BioContainers

                 {{0-1}}
************************************************

```bash
################## BASE IMAGE ######################

FROM biocontainers/biocontainers:v1.0.0_cv4

################## METADATA ######################

LABEL base_image="biocontainers:v1.0.0_cv4"
LABEL version="2"
LABEL software="NCBI BLAST+"
LABEL software.version="2.2.31"
LABEL about.summary="basic local alignment search tool"
LABEL about.home="http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastHome"
LABEL about.documentation="http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastHome"
LABEL about.license_file="https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/scripts/projects/blast/LICENSE"
LABEL about.license="SPDX:MIT"
LABEL extra.identifiers.biotools="BLAST"
LABEL about.tags="Genomics"

################## MAINTAINER ######################

MAINTAINER Saulo Alves Aflitos <sauloal@gmail.com>

################## INSTALLATION ######################

RUN conda install blast=2.2.31

WORKDIR /data/
```

************************************************

                 {{1}}
************************************************

[repo for biocontainers base image](https://github.com/BioContainers/containers/blob/master/biocontainers/1.0.0/Dockerfile)

```
# Base image
FROM ubuntu:16.04

################## METADATA ######################

LABEL base_image="ubuntu:16.04"
LABEL version="4"
LABEL software="Biocontainers base Image"
LABEL software.version="1.0.0"
LABEL about.summary="Base image for BioDocker"
LABEL about.home="http://biocontainers.pro"
LABEL about.documentation="https://github.com/BioContainers/specs/wiki"
LABEL about.license_file="https://github.com/BioContainers/containers/blob/master/LICENSE"
LABEL about.license="SPDX:Apache-2.0"
LABEL about.tags="Genomics,Proteomics,Transcriptomics,General,Metabolomics"

################## MAINTAINER ######################
MAINTAINER Felipe da Veiga Leprevost <felipe@leprevost.com.br>

ENV DEBIAN_FRONTEND noninteractive

RUN mv /etc/apt/sources.list /etc/apt/sources.list.bkp && \
    bash -c 'echo -e "deb mirror://mirrors.ubuntu.com/mirrors.txt xenial main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-updates main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-backports main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-security main restricted universe multiverse\n\n" > /etc/apt/sources.list' && \
    cat /etc/apt/sources.list.bkp >> /etc/apt/sources.list && \
    cat /etc/apt/sources.list

RUN apt-get clean all && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y  \
        autotools-dev   \
        automake        \
        cmake           \
        curl            \
        grep            \
        sed             \
        dpkg            \
        fuse            \
        git             \
        wget            \
        zip             \
        openjdk-8-jre   \
        build-essential \
        pkg-config      \
        python          \
	python-dev      \
        python-pip      \
        bzip2           \
        ca-certificates \
        libglib2.0-0    \
        libxext6        \
        libsm6          \
        libxrender1     \
        git             \
        mercurial       \
        subversion      \
        zlib1g-dev &&   \
        apt-get clean && \
        apt-get purge && \
        rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.0.5-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

RUN TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

RUN mkdir /data /config

# Add user biodocker with password biodocker
RUN groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse biodocker && \
    echo `echo "biodocker\nbiodocker\n" | passwd biodocker` && \
    chown biodocker:biodocker /data && \
    chown biodocker:biodocker /config

# give write permissions to conda folder
RUN chmod 777 -R /opt/conda/

# Change user
USER biodocker

ENV PATH=$PATH:/opt/conda/bin
ENV PATH=$PATH:/home/biodocker/bin
ENV HOME=/home/biodocker

RUN mkdir /home/biodocker/bin

RUN conda config --add channels r
RUN conda config --add channels bioconda

RUN conda upgrade conda

VOLUME ["/data", "/config"]

# Overwrite this with 'CMD []' in a dependent Dockerfile
CMD ["/bin/bash"]

WORKDIR /data
```

************************************************

### How to run the docker image

```
 $ cd /home/user/workplace
 $ docker pull biocontainers/blast
 $ docker run biocontainers/blast blastp -help
 $ wget http://www.uniprot.org/uniprot/P04156.fasta    
 $ curl -O ftp://ftp.ncbi.nih.gov/refseq/D_rerio/mRNA_Prot/zebrafish.1.protein.faa.gz
 $ gunzip zebrafish.1.protein.faa.gz
 $ docker run -v /Users/yperez/workplace:/data/ biocontainers/blast makeblastdb -in zebrafish.1.protein.faa -dbtype prot
 $ docker run -v /Users/yperez/workplace:/data/ biocontainers/blast blastp -query P04156.fasta -db zebrafish.1.protein.faa -out results.txt
```

## Other examples from BioContainers

                 {{0}}
************************************************

![BioContainers Architecture](https://onlinelibrary.wiley.com/cms/attachment/d7a972cb-8975-4ce4-ba57-5ba89f000939/pmic13223-fig-0003-m.jpg)

[^1](https://onlinelibrary.wiley.com/doi/full/10.1002/pmic.201900147)
[^2](https://hub.docker.com/u/biocontainers)
[^3](https://quay.io/organization/biocontainers)
[^4](http://biocontainers.pro/regitry)

************************************************

                --{{0}}--
BioContainers architecture from the container request by the user in GitHub to the final container deposited in DockerHub and Quay.io. The BioContainers community in collaboration with the BioConda community defines a set of guidelines and protocols to create a Conda and Docker container including mandatory metadata, tests, and trusted images. The proposed architecture uses a continuous integration system to test and build the final containers and deposit them into public registries. All the Containers and tools can be searched from the BioContainers registry.

## Reproducibility stack

![reproducibility stack](http://data.bits.vib.be/pub/trainingen/CommunityMeeting/20180627-reproducibility-stack.png)

[Reference](https://www.biorxiv.org/content/early/2017/10/11/200683)

## Recommendations

- Carefully define a set of tools for a given analysis
- Use tools from the Bioconda registry
- Adopt containers to guarantee consistency of results
- Use virtualization to make analyses “resistant to time”.

## Further reading

- [impact of docker containers on performance](https://peerj.com/articles/1273/)
- [container-based virtualization for HPC environments](https://arxiv.org/abs/1709.10140)
- recommendations on [containers](https://f1000research.com/articles/7-742/v1)
- [practical computational reproducibility in life sciences](https://www.biorxiv.org/content/early/2017/10/11/200683)

# Thanks

- Bioinfo Core at CRG  [slides](https://github.com/biocorecrg/C4LWG-2018/tree/master/slides)
- based on recommendation from [F1000](https://f1000research.com/articles/7-742/v1)

## devops

- software development + software operations
- automate and monitor

![Devops Explained](https://hostadvice.com/wp-content/uploads/2018/03/devopsext.jpg)
