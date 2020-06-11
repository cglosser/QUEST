#!/bin/bash
cd /tmp
wget https://bitbucket.org/eigen/eigen/get/3.3.2.tar.gz -O eigen.tar.gz
mkdir -p eigen/build && tar -xvzf eigen.tar.gz -C eigen --strip-components 1
cd eigen/build && cmake .. && make install
