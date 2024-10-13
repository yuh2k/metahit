conda env create -f env.yaml

# BWA
wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2/download -O bwa-0.7.17.tar.bz2
tar -xjf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
cd ..
