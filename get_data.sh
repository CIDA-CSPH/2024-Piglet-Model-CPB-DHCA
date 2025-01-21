home_dir=$(dirname $0)
cd $home_dir


wget https://zenodo.org/records/14713666/files/data_results.zip
unzip data_results.zip
rm -rf __MACOSX
mv data_results/data .
mv data_results/results .
rm -rf data_results
rm data_results.zip