cd src/
make clean
make
cd ..
python3 bamulator.py --variants scripts/sv_config.tsv --output /raw-data/projectes/SPIKE_IN_BAM/SIMULATED/ --input_vcf_list /raw-data/projectes/SPIKE_IN_BAM/BAM_PANELS --map_dir /raw-data/ANNOTATIONS/ANN_DIR/1000Genomes/hg19/beagle/genetic_maps/plink.GRCh37.map --ref_dir /raw-data/ANNOTATIONS/ANN_DIR/1000Genomes/hg19/beagle/references --suffix .simulated --threads 6 --reference /raw-data/ANNOTATIONS/REF_DIR/hg19/ucsc.hg19.fasta