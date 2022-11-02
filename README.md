# NIPH_TB_pipeline
The pipeline developed for analysis of MTB analysis at NIPH

Supports:
- Quality assessment (FastQC, mash (prev kaiju), read depth)
- Typing using Coll (2014) scheme
- AMR prediction using mykrobe predictor
- Phylogenetic tree of nearest neighbors

Developed by Ola Brynildsrud (olbb@fhi.no)

## Setting appropriate storage location for containers
```
docker system prune -a
sudo service docker stop
# Create file /etc/docker/daemon.json:
{
    "graph": "/Path/to/docker/containers"
}
sudo service docker start
```
## Some DBs exist in separate volume. Need to build
Currently, the relevant file lives in "F:/Omr2/Felles/Avdeling for Ola, Jon og Vegard/Seksjon for Ola/TB_pipeline/Kaijudb". From that directory:
```
docker run -it --rm -v kaijudb:/volume -v "$(pwd)":/backup alpine sh -c "rm -rf /volume/* /volume/..?* /volume.[!.]* ; tar -C /volume/ -xvf /backup/backup.tar; cp Trimal* /volume/Reference/"
```

## Building

```
git clone https://github.com/folkehelseinstituttet/niph_tb_pipeline
docker build .
```

## Usage
From directory with all FASTQ folders. (Note - fill in <Path to FELLES GLOBAL> - That is the directory with all previously run files. Currently "N:/NGS/TB_pipeline/TB_pipeline_database/DB"
```
docker run -v kaijudb:/mnt -v "$(pwd)":"/data" -v '</Path to FELLES GLOBAL>':'/mnt/global_collection' folkehelseinstituttet/niph_tb_pipeline:latest /bin/bash -c "cd /data && niph_tb_pipeline"
```
