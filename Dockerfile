FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
        build-essential \
        curl \
        lbzip2 \
        libboost-all-dev \
        libcurl3-dev \
        libgsl-dev \
        libhdf5-serial-dev \
        openjdk-8-jdk \
        python3 \
        python3-pip \
        unzip \
        vim-common \
        wget \
        zlib1g-dev \
		r-base \
		libxml2-dev \
		libssl-dev \
		vcftools \
		-y git \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/
	
# workaround for PEER, see https://github.com/mz2/peer/issues/4
RUN apt-get update \
    && apt-get install -y \
        gcc-5 \
        g++-5 \
        gfortran-5 \
        cmake \
    && rm -rf /var/lib/apt/lists/* \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 50 --slave /usr/bin/g++ g++ /usr/bin/g++-5
	

# R
RUN wget https://github.com/downloads/PMBio/peer/R_peer_source_1.3.tgz && \
    R CMD INSTALL R_peer_source_1.3.tgz && \
    rm R_peer_source_1.3.tgz && \
    echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "install.packages(c('argparser'), dependencies=TRUE)" && \
    Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("qvalue"); biocLite("sva"); biocLite("edgeR");'
	
# standalone R math library
RUN wget https://cran.r-project.org/src/base/R-3/R-3.4.4.tar.gz && tar -xf R-3.4.4.tar.gz && rm R-3.4.4.tar.gz \
    && cd R-3.4.4 && ./configure --with-x=no && cd src/nmath/standalone && make
ENV RMATH /R-3.4.4/src

# Python
RUN pip3 install --upgrade pip setuptools
RUN pip3 install numpy tables pandas scipy matplotlib h5py feather-format pysam statsmodels scikits.bootstrap xlrd
# numpy dependencies:
RUN pip3 install pyBigWig

# htslib
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar -xf htslib-1.9.tar.bz2 && rm htslib-1.9.tar.bz2 && cd htslib-1.9 && make && make install && make clean

# bcftools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
    tar -xf bcftools-1.9.tar.bz2 && rm bcftools-1.9.tar.bz2 && cd bcftools-1.9 && \
    ./configure --with-htslib=/opt/htslib-1.9 && make && make install && make clean

# samtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -xf samtools-1.9.tar.bz2 && rm samtools-1.9.tar.bz2 && cd samtools-1.9 && \
    ./configure --with-htslib=/opt/htslib-1.9 && make && make install && make clean

# PLINK
RUN mkdir /opt/plink && cd /opt/plink && \
    wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200107.zip && \
    unzip plink_linux_x86_64_20200107.zip && rm plink_linux_x86_64_20200107.zip
ENV PATH $PATH:/opt/plink


ENV PATH /opt/fastqtl/bin:$PATH
ENV PATH /usr/bin:$PATH

RUN Rscript -e "install.packages(c('dplyr','httr'), dependencies=TRUE)"
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("GEOquery");'
RUN Rscript -e "install.packages(c('devtools', 'MatrixEQTL'), dependencies=TRUE)"
RUN Rscript -e "devtools::install_github('foster-gabe/PFExpTools', ref = 'backcompat')"


RUN wget --no-check-certificate https://qtltools.github.io/qtltools/binaries/QTLtools_1.2_source.tar.gz && \
    tar -xf QTLtools_1.2_source.tar.gz && rm QTLtools_1.2_source.tar.gz
COPY /pfalsrc/Makefile QTLtools_1.2_source/
RUN cd QTLtools_1.2_source && make
RUN cd QTLtools_1.2_source/script $$ sed -i 's/Ph\$7/Ph\$V19/' runFDR_ftrans.R && sed -i 's/Nh\$7/Nh\$V12/' runFDR_ftrans.R
ENV PATH "$PATH:/QTLtools_1.2_source/bin"

RUN git clone https://github.com/francois-a/fastqtl
RUN cd fastqtl && mkdir obj && sed -i 's/RMATH=/#RMATH/' Makefile && make
ENV PATH "$PATH:/fastqtl/bin"

RUN pip install multiprocess


# copy scripts
COPY pfalsrc pfalsrc/
