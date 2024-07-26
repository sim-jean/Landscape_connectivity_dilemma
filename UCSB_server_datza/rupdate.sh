#!/bin/bash

# Add the CRAN repository to your sources list
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

# Update the package list
sudo apt-get update

# Install the latest version of R
sudo apt-get install -y r-base

# Verify the installation
R --version
