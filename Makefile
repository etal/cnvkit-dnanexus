#SHELL=/bin/bash -e -x -o pipefail
all:
	sudo apt update
	sudo apt --yes install \
		cython3 \
		libblas-dev \
		libfreetype6-dev \
		libgfortran-4.8-dev \
		liblapack-dev \
		liblzma-dev \
		libpng-dev \
		pkg-config \
		python3-dev \
		python3-matplotlib \
		python3-numpy \
		python3-pip
	sudo -H pip3 install --upgrade pip
	sudo -H pip3 install -U pandas reportlab pomegranate
	sudo -H pip3 install cnvkit

