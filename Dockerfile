FROM continuumio/anaconda3

## INSTALL system dependencies
RUN apt-get update
RUN apt-get -y install apt-utils
RUN apt-get -y install make unzip gcc g++ zlib1g-dev figtree imagemagick locales-all

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

## INSTALL TeX stuff
RUN apt-get -y install texlive-full

## INSTALL python dependencies
RUN pip install ete3 Biopython

## SWITCH TO MAMBA
RUN conda install -c conda-forge mamba

## INSTALL conda stuff
RUN mamba config --add channels defaults
RUN mamba config --add channels conda-forge
RUN mamba config --add channels bioconda
RUN mamba update -y -n base conda
RUN mamba install -y snp-dists snippy fastqc snp-sites fasttree trimal mash rapidnj
RUN mamba install -y mykrobe mykatlas

## UPDATE MYKROBE
RUN mamba upgrade -c bioconda -n mykrobe mykrobe

## REINSTALL SNIPPY - FAULTY BUILD
RUN mamba create -y -n snippy -c bioconda snippy
## FORCE SNPEFF TO VERSION 4, VERSION 5.1 THROWS ERROR
RUN mamba install -y -n snippy -c bioconda snpEff=4

## ADD MASH DATABASE
RUN mkdir -p /mnt/mash
WORKDIR /media/mash
RUN wget https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh
WORKDIR /

## REINSTALL COLLTYPER TO CURRENT PYTHON
## ALSO NEED TO DOWNGRADE SETUPTOOLS FOR THIS
RUN pip install "setuptools<58.0.0"
RUN pip install --upgrade git+https://github.com/admiralenola/colltyper

WORKDIR /opt/niph_tb_pipeline
COPY niph_tb_pipeline .
RUN pip install --upgrade .
WORKDIR /

## MAKE new user with no root privileges
RUN useradd -ms /bin/bash -g sudo tbuser
USER tbuser

MAINTAINER Ola Brynildsrud "olbb@fhi.no"
