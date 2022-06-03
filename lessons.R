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
#cd folder containing the file 
#docker run -t -v `pwd`:`pwd` -w `pwd` ewels/multiqc

#create image of ubuntu with specific requisiti
#docker pull ubuntu:21.04
#docker run -i -t ubuntu:21.04
#apt-get update/ upgrade/ apt-get install sudo / apt-get install nano
#once done exit and docker commit "istance in which we have done what we need" ubuntu:21.04.angelo (this final point provided the tag included)

#create an R file with the command install.packages('umap', repos='http://cran.us.r-project.org')#saved as umap_command.R in the same folder of Dockerfile
#create the dockerfile and the container from it 
#mkdir DockerFile (generate the folder)
#touch Dockerfile (generate the file with no extension)
#vim Dockerfile (go in the file)
#FROM rocker/r-ubuntu:20.04	
#RUN apt-get -y update
#RUN apt-get -y upgrade
#RUN apt-get -y install sudo 
#RUN apt-get -y install nano
#COPY ./umap_command.R/ home/
#RUN chmod +x /home/umap_command.R
#RUN sudo apt-get update && sudo apt-get -y install libssl-dev

#i to write and esc :wq! to exit saving
###After creating the Dockerfile goes on the terminal and type the command below.
#docker build Dockerfile -t r4:v.0.01 #build the docker 
#docker run -i -t r4:v.0.01 #run the docker
#sudo ls #check that sudo is present and work
#Then go in the home dir
#Rscript ./umap_command.R #run the R file to install umap exit
#docker ps -a #find the instance name ID of the instance previously closed 
#docker commit instaceID (find to be in my case d7a30dac5e11 through docker ps -a) r4:v.0.01 #commit the changes to the docker
#docker run -i -t r4:v.0.01 #run again the docker
#sudo ls #check that the changes were saved
#exit