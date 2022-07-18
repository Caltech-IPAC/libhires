FROM debian:bullseye

RUN echo "\
# Base repository \n\
deb http://deb.debian.org/debian bullseye main contrib non-free \n\
deb-src http://deb.debian.org/debian bullseye main contrib non-free \n\
deb http://deb.debian.org/debian-security/ bullseye-security main contrib non-free \n\
deb-src http://deb.debian.org/debian-security/ bullseye-security main contrib non-free \n\
deb http://deb.debian.org/debian bullseye-updates main contrib non-free \n\
deb-src http://deb.debian.org/debian bullseye-updates main contrib non-free \n\
" > sources.list \
    && mv -f sources.list /etc/apt/sources.list
RUN apt-get update

RUN apt-get install -y \
    ant \
    autoconf \
    automake \
    binutils \
    bison \
    cabextract \
    checkstyle \
    coreutils \
    cscope \
    csh \
    curl \
    diffutils \
    dos2unix \
    doxygen \
    flex \
    fonts-freefont-ttf \
    g++ \
    gawk \
    gcc \
    gdb \
    gearman \
    gettext \
    gfortran \
    gimp \
    git \
    git-cola \
    git-doc \
    git-gui \
    gitg \
    gnuplot \
    gnuplot-nox \
    golang \
    graphviz \
    grep \
    gsoap \
    gzip \
    hdf5-tools \
    html2ps \
    htmldoc \
    htop \
    iftop \
    imagemagick \
    imagemagick-doc \
    indent \
    iotop \
    iperf \
    iperf3 \
    javahelp2 \
    jq \
    ksh \
    libguice-java \
    libjs-sphinxdoc \
    libpqxx-doc \
    libprotobuf-java \
    libtool \
    libvsqlitepp-doc \
    libxerces2-java \
    libxml2-utils \
    libzookeeper-java \
    m4 \
    magit \
    make \
    memcached \
    nicstat \
    odbc-postgresql \
    patch \
    pdftk \
    pdl \
    perl \
    pgplot5 \
    pkg-config \
    plotutils \
    pwgen \
    qgit \
    rcs \
    sed \
    sloccount \
    source-highlight \
    sqlite3 \
    strace \
    swig \
    tar \
    tcpdump \
    testng \
    texinfo \
    tidy \
    unzip \
    vim-doc \
    vim-nox \
    vim-scripts \
    wcslib-doc \
    wget \
    zip \
    libaio-dev \
    libarmadillo-dev \
    libboost-all-dev \
    libboost-dev \
    libcairo2-dev \
    libccfits-dev \
    libctemplate-dev \
    libcurl4-gnutls-dev \
    libeigen3-dev \
    libelf-dev \
    libevent-dev \
    libexpat1-dev \
    libexpat1-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libglib2.0-dev \
    libgmp-dev \
    libgraphviz-dev \
    libgsl0-dev \
    libgtest-dev \
    libhdf5-dev \
    libjpeg-dev \
    libjsoncpp-dev \
    liblapack-dev \
    liblzo2-dev \
    libmemcached-dev \
    libmlpack-dev \
    libmpfr-dev \
    libncurses5-dev \
    libpango1.0-dev \
    libpcre3-dev \
    libpixman-1-dev \
    libplplot-dev \
    libpqxx-dev \
    libprotobuf-dev \
    libreadline-dev \
    libsqlite3-dev \
    libssl-dev \
    libvsqlitepp-dev \
    libxml2-dev \
    libxslt1-dev \
    python-dev \
    python3-dev \
    unixodbc-dev \
    uuid-dev \
    wcslib-dev \
    zlib1g-dev \
    libgsoap-dev \
    libwcstools-dev \
    gedit 

RUN apt-get install -y \
    libgd-dev \
    libcgicc-dev 

RUN ln -sf /usr/bin/python3 /usr/bin/python

ARG token
RUN git clone -b debian11 https://$token:@github.com/IPAC-SW/ipac-base.git ipac-base
RUN git clone -b debian11 https://$token:@github.com/IPAC-SW/irsa-libcells.git irsa-libcells
RUN git clone -b debian11 https://$token:@github.com/IPAC-SW/irsa-libjbcvtc.git irsa-libjbcvtc
RUN git clone -b debian11 https://$token:@github.com/IPAC-SW/irsa-libned.git irsa-libned
RUN git clone -b master   https://$token:@github.com/IPAC-SW/irsa-libscew.git irsa-libscew
RUN git clone -b debian11 https://$token:@github.com/IPAC-SW/irsa-libsim-1.0.git irsa-libsim
RUN git clone -b debian11 https://$token:@github.com/IPAC-SW/irsa-libim.git irsa-libim
RUN git clone -b debian11 https://$token:@github.com/IPAC-SW/irsa-hcompress.git irsa-hcompress
RUN git clone -b debian11 https://$token:@github.com/IPAC-SW/irsa-cspice.git irsa-cspice

RUN cd /irsa-libscew && ./configure --prefix=/usr && make && make install
RUN cd /irsa-libjbcvtc && ./waf configure --prefix=/usr && ./waf build install
RUN cd /irsa-libcells && ./waf configure --prefix=/usr && ./waf build install
RUN cd /irsa-libned && ./waf configure --prefix=/usr && ./waf build install
RUN cd /irsa-libsim && ./waf configure --prefix=/usr && ./waf build install
RUN cd /irsa-libim && ./waf configure --prefix=/usr && ./waf build install
RUN cd /irsa-hcompress && ./waf configure --prefix=/usr && ./waf build install
RUN cd /irsa-cspice && ./waf configure --prefix=/usr && ./waf build install

ARG user
ENV USER=$user

ENV CM_ENV_DIR=/usr
ENV INFORMIXDIR=/usr/local/informix_client
ENV INFORMIXSERVER=rmt_dbms20

ENV ORACLE_HOME=/usr/lib/oracle/19.3/client64
ENV ODBC_HOME=/usr
ENV ODBCINI=/etc/odbc.ini
ENV TNS_ADMIN=/etc/oracle

ENV CM_TPS_DIR=/usr
ENV CM_STK_DIR=/irsa/cm/website/stk
ENV CM_MONTAGE_DIR=/irsa/cm/website/montage
ENV CM_OASIS_DIR=/irsa/cm/website/oasis

ENV DEFDBMS=ORACLE
ENV PYTHONPATH=$CM_STK_DIR/python

ENV LD_LIBRARY_PATH=/usr/lib/slurm-drmaa/lib/:${ORACLE_HOME}/lib

ENV JAVA_HOME=/usr/lib/jvm/java-1.11.0-openjdk-amd64
ENV JDKVERSION="11.0.5"

ENV CM_BASE_DIR=/irsa/cm/ws/$USER/base
RUN mkdir -p $CM_BASE_DIR
ENV CM_BASE_APP_DIR=/irsa/cm/ws/$USER/irsa
RUN mkdir -p $CM_BASE_APP_DIR
ENV CM_IRSA_DIR=/irsa/cm/ws/$USER/irsa

ENV PATH=/usr/bin:${PATH}
ENV PATH=${CM_BASE_DIR}/bin:${PATH}
ENV PATH=${CM_IRSA_DIR}/bin:${PATH}

ARG CACHE_DATE=2016-01-01
RUN cd /ipac-base && git pull
RUN cd /ipac-base/src && make 
RUN cp -rf /ipac-base/bin /ipac-base/lib /ipac-base/include /ipac-base/share /ipac-base/docs/svc /ipac-base/docs/lib ${CM_BASE_DIR}/

#ARG CACHE_DATE=2016-01-01
RUN mkdir /irsa-cgi
ADD . /irsa-cgi/
RUN cd /irsa-cgi &&./waf configure && ./waf clean && ./waf build
#RUN ./waf configure --libxml2-libdir=/usr/lib/aarch64-linux-gnu/  --libxml2-incdir=/usr/include/libxml2 && ./waf -j1 
#BUILD
	# docker build --network=host --build-arg CACHE_DATE="$(date)" --build-arg token=${PERSONAL_GITHUB_TOKEN} --build-arg user=$USER --tag irsa-cgi:latest .

#RUN
	#LINUX X11 server
		# docker run --privileged -it --net=host -e DISPLAY=:0 -v /tmp/.X11-unix:/tmp/.X11-unix irsa-cgi:latest bash
	
	#MAC
	
		# .bashrc entry for initializing display and xquartz
			#export DISPLAY_MAC=`ifconfig en0 | grep "inet " | cut -d " " -f2`:0
			
			#function startx() 
			#{
			#	if [ -z "$(ps -ef|grep Xquartz|grep -v grep)" ] ; then
			#         socat TCP-LISTEN:6000,reuseaddr,fork UNIX-CLIENT:\"$DISPLAY\" &
			#	fi
			#}
		# xhost +
		# docker run --privileged -it --net=host -e DISPLAY=$DISPLAY_MAC sandbox:latest
