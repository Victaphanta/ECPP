#!/bin/bash
cd /usr/local/src
apt-get update

apt-get --assume-yes install build-essential git wget python python-pip pigz libz-dev default-jre unzip libtbb-dev libncurses-dev curl libssl-dev gcc autoconf automake pkg-config cmake jellyfish salmon

pip install dendropy
pip install numpy
pip install biopython

echo export PATH=$PATH:/usr/local/bin/spade/bin:/usr/local/bin/trinity:/usr/local/bin/bbmap:/usr/local/bin/bowtie2:/usr/local/bin/blast+/bin:/usr/local/share/util >> ~/.bashrc
echo "alias varscan='java -jar /usr/local/bin/varscan.jar'" >> ~/.bashrc
echo "alias trimmomatic='java -jar /usr/local/bin/trimmomatic.jar'" >> ~/.bashrc

git clone https://github.com/Victaphanta/ECPP.git
chmod -R 777 /usr/local/src/ECPP
mv ECPP/ECPP_1.2.1.sh /usr/local/bin
mv ECPP/EssentialScripts/* /usr/local/bin
rm -rf ECPP

git clone https://github.com/lh3/seqtk.git
make -C seqtk/
mv seqtk/seqtk /usr/local/bin/seqtk
#rm -rf seqtk

wget 'http://mafft.cbrc.jp/alignment/software/mafft-7.310-without-extensions-src.tgz'
pigz -dc mafft-7.310-without-extensions-src.tgz | tar xf - 
make -C mafft-7.310-without-extensions/core/
make install -C mafft-7.310-without-extensions/core/
rm -rf mafft-7.310-without-extensions

wget 'https://nchc.dl.sourceforge.net/project/fastuniq/FastUniq-1.1.tar.gz'
pigz -dc FastUniq-1.1.tar.gz | tar xf -
make -C FastUniq/source/
mv FastUniq/source/fastuniq /usr/local/bin/fastuniq
rm -rf FastUniq

wget 'https://nchc.dl.sourceforge.net/project/bbmap/BBMap_36.92.tar.gz'
pigz -dc BBMap_36.92.tar.gz | tar xf -
mv  bbmap /usr/local/bin/bbmap/
rm -rf bbmap 

wget 'http://cab.spbu.ru/files/release3.10.1/SPAdes-3.10.1-Linux.tar.gz'
pigz -dc SPAdes-3.10.1-Linux.tar.gz  | tar xf -
mv SPAdes-3.10.1-Linux /usr/local/bin/spade/
rm -rf SPAdes-3.10.1-Linux

wget 'https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.8.4.tar.gz'
pigz -dc Trinity-v2.8.4.tar.gz| tar xf -
make -C trinityrnaseq-Trinity-v2.8.4/ 
make plugins -C trinityrnaseq-Trinity-v2.8.4/ 
mv trinityrnaseq-Trinity-v2.8.4 /usr/local/bin/trinity/
rm -rf trinityrnaseq-Trinity-v2.8.4

wget 'https://github.com/COMBINE-lab/salmon/releases/download/v0.12.0-alpha/salmon-latest_linux_x86_64.tar.gz'
pigz -dc salmon-latest_linux_x86_64.tar.gz| tar xf -
sudo mv salmon-latest_linux_x86_64/bin/* /usr/local/bin/
sudo mv salmon-latest_linux_x86_64/lib/* /usr/local/lib/

curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
unzip Trimmomatic-0.36.zip
mv Trimmomatic-0.36/trimmomatic-0.36.jar /usr/local/bin/trimmomatic.jar
mv /usr/local/src/Trimmomatic-0.36/adapters /usr/local/share/trimmomatic/
rm -rf Trimmomatic-0.36 

wget 'https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.3/bowtie2-2.3.4.3-source.zip'
unzip  bowtie2-2.3.4.3-source.zip
make -C bowtie2-2.3.4.3/
mv bowtie2-2.3.4.3 /usr/local/bin/bowtie2/
rm -rf bowtie2-2.3.4.3

git clone https://github.com/dkoboldt/varscan.git
mv varscan/VarScan.v2.4.3.jar /usr/local/bin/varscan.jar

wget 'https://nchc.dl.sourceforge.net/project/samtools/samtools/1.3.1/samtools-1.3.1.tar.bz2'
tar jxf samtools-1.3.1.tar.bz2
make -C samtools-1.3.1/
make install -C samtools-1.3.1/
rm -rf samtools-1.3.1

git clone https://github.com/vcftools/vcftools.git
cd vcftools
./autogen.sh
./configure
make
make install
cd ..
rm -rf vcftools

wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.7.1+-x64-linux.tar.gz'
pigz -d ncbi-blast-2.7.1+-x64-linux.tar.gz 
tar xvpf ncbi-blast-2.7.1+-x64-linux.tar
mv ncbi-blast-2.7.1+ /usr/local/bin/blast+/
#rm -rf ncbi-blast-2.7.1+-x64-linux.tar 

wget http://last.cbrc.jp/last-869.zip
unzip last-869.zip
make -C last-869/
make install -C last-869/
rm -rf last-869

##############################################
#rm -rf /var/lib/apt/lists/*
#/bin/bash -c 'source .bashrc'
##############################################
