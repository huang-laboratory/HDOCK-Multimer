# Bash script for installing the tools required by HDM

# Install FFT library
if ! ldconfig -p | grep -q "libfftw3.so.3"; then
    apt-get install -y -qq libfftw3-double3
fi

# Install HDOCKlite
wget -q https://github.com/huang-laboratory/HDOCKlite/archive/refs/heads/main.zip -O HDOCKlite-main.zip
unzip -q HDOCKlite-main.zip
cp HDOCKlite-main/createpl tools/
cp HDOCKlite-main/hdock tools/
rm -rf HDOCKlite-main.zip HDOCKlite-main/

# Install HSYMDOCK
wget -q https://github.com/huang-laboratory/HSYMDOCK/archive/refs/heads/main.zip -O HSYMDOCK-main.zip
unzip -q HSYMDOCK-main.zip
cp HSYMDOCK-main/chdock tools/
cp HSYMDOCK-main/compcn tools/
cp HSYMDOCK-main/dhdock tools/
cp HSYMDOCK-main/dhdock.sh tools/
cp HSYMDOCK-main/compdn tools/
rm -rf HSYMDOCK-main.zip HSYMDOCK-main/

# Install jq
wget -q https://github.com/jqlang/jq/releases/download/jq-1.8.1/jq-linux64
chmod +x jq-linux64
mv jq-linux64 tools/jq

# Install MMalign
wget -q https://github.com/pylelab/USalign/archive/refs/heads/master.zip -O USalign-master.zip
unzip -q USalign-master.zip
cd USalign-master/
g++ -O3 -ffast-math -lm -o MMalign MMalign.cpp
mv MMalign ../tools/
cd ../ && rm -rf USalign-master.zip USalign-master/

# Install STRIDE
wget -q https://github.com/heiniglab/stride/archive/refs/heads/main.zip -O stride-main.zip
unzip -q stride-main.zip
mv stride-main/stride tools/
rm -rf stride-main.zip stride-main/

