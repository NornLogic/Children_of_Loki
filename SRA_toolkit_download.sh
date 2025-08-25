#install SRA toolkit 
#download from: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

#update these
toolkit_version="sratoolkit.3.2.1-ubuntu64"
toolkit_untar_file="sratoolkit.current-ubuntu64.tar.gz"
toolkit_path="/home/rstudio"
wget_path="https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz"

#download, i don't like using sudo though you can do this with sudo apt-get install
wget $wget_path

#extract it
tar -xvzf $toolkit_untar_file

#add to path
echo 'export PATH=$PATH:$toolkit_path/$toolkit_untar_file/bin' >> ~/.bashrc
source ~/.bashrc

#test it
fastq-dump --version

# Create config directory, make sure Remote Access is enabled 
vdb-config --interactive

#make download dir
mkdir ~/sra_downloads

#set up where downloads are going to go
export SRA_CACHE=~/sra_downloads
export VDB_CONFIG=~/sra_downloads

#add to bashrc
echo 'export SRA_CACHE=~/sra_downloads' >> ~/.bashrc
echo 'export VDB_CONFIG=~/sra_downloads' >> ~/.bashrc
source ~/.bashrc




