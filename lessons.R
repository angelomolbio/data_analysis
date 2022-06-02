library(devtools)
install_github("kendomaniac/docker4seq", ref="master")
library(docker4seq)
downloadContainers(group="docker")
install_github("kendomaniac/rCASC", ref="master")
library(rCASC)
downloadContainers(group="docker")

#Fastqcsetting 
#http://users.path.ox.ac.uk/~pcook/w1/NGS_workshop2.htm
#sudo ln -s /home/angeloscarciglia/Desktop/calogero/FastQC/fastqc /usr/local/bin/fastqc
#chmod + fastqc in FastQC folder downloaded
#./fastqc
# in cd dataset0 cartella del file -> fastqc ra100k.R1.fastq

#MultiQC directly from the terminal with docker container 
#docker run -t -v `pwd`:`pwd` -w `pwd` ewels/multiqc

