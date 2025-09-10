cd src/
make clean
make
cd ..
python3 spikeinbam.py --variants scripts/sv_config.tsv --output /raw-data/projectes/SPIKE_IN_BAM/SIMULATED/ --suffix .simulated --threads 6 --reference /raw-data/ANNOTATIONS/REF_DIR/hg19/ucsc.hg19.fasta



