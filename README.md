# Foldseek 
Foldseek enables fast and sensitive comparisons of large structure sets.

<p align="center"><img src="https://github.com/steineggerlab/foldseek/blob/master/.github/foldseek.png" height="250"/></p>

## Publications

[van Kempen M, Kim S, Tumescheit C, Mirdita M, Söding J, and Steinegger M. Foldseek:  fast and accurate protein structure search. bioRxiv, doi:10.1101/2022.02.07.479398  (2022)](https://www.biorxiv.org/content/10.1101/2022.02.07.479398)

## Webserver 
Search your protein structures against the [AlphaFoldDB](https://alphafold.ebi.ac.uk/) and [PDB](https://www.rcsb.org/) in seconds using our Foldseek webserver: [search.foldseek.com](https://search.foldseek.com) 🚀

## Installation

    # static Linux AVX2 build (check using: cat /proc/cpuinfo | grep avx2)
    wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz; tar xvzf foldseek-linux-avx2.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH
    # static Linux SSE4.1 build (check using: cat /proc/cpuinfo | grep sse4_1)
    wget https://mmseqs.com/foldseek/foldseek-linux-sse41.tar.gz; tar xvzf foldseek-linux-sse41.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH
    # static macOS build (universal binary with SSE4.1/AVX2/M1 NEON)
    wget https://mmseqs.com/foldseek/foldseek-osx-universal.tar.gz; tar xvzf foldseek-osx-universal.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH
    # conda installer 
    conda install -c conda-forge -c bioconda foldseek

Other precompiled binaries for ARM64, PPC64LE amd SSE2 are available at [https://mmseqs.com/foldseek](https://mmseqs.com/foldseek).

### Quick start
`easy-search` can search single or multiple query structures formatted in PDB/mmCIF format (flat or `.gz`) against a target database (`example/`) of protein structures. It outputs a tab-separated file of the alignments (`.m8`) the fields are `query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits`.

    foldseek easy-search example/d1asha_ example/ aln.m8 tmpFolder
    
The output can be customized with the `--format-output` option e.g. `--format-output "query,target,qaln,taln"` returns the query and target accession and the pairwise alignments in tab separated format. You can choose many different [output columns](https://github.com/soedinglab/mmseqs2/wiki#custom-alignment-format-with-convertalis).    

The target database can be pre-processed by `createdb`. This make sense if searched multiple times. 
 
    foldseek createdb example/ targetDB
    foldseek easy-search example/d1asha_ targetDB aln.m8 tmpFolder

### Important search parameters
    # sensitivity 
    -s                       adjust the sensitivity to speed trade-off.
                             lower is faster, higher more sensitive (fast: 7.5, highest sensitivity (default): 9.5)
    --max-seqs               adjust the amount of prefilter that are handed to the alignment. 
                             Increasing it can lead to more hits (default: 300)
    # other                         
    --alignment-type         0: 3Di Gotoh-Smith-Waterman (local, not recommended), 
                             1: TMalign (global), 
                             2: 3Di+AA Gotoh-Smith-Waterman (local, default)
    -c                       list matches above this fraction of aligned (covered) residues (see --cov-mode) (default: 0.0) 
    --cov-mode               0: coverage of query and target, 1: coverage of target, 2: coverage of query

### Databases 
The `databases` command downloads pre-generated databases like PDB or AlphaFoldDB.
    
    # pdb  
    foldseek databases PDB100 pdb tmp 
    # alphafold db
    foldseek databases Alphafold/Proteome afdb tmp 

We currently support the following databases: 
```
  Name                  Type            Taxonomy        Url
- Alphafold/Proteome    Aminoacid            yes        https://alphafold.ebi.ac.uk/
- Alphafold/Swiss-Prot  Aminoacid            yes        https://alphafold.ebi.ac.uk/
- PDB100                Aminoacid            yes        https://www.rcsb.org
```

### Main Modules
- `easy-search`       fast protein structure search  
- `createdb`          create a database from protein structures (PDB,mmCIF, mmJSON)
- `databases`         download pre-assembled databases

### Using TMalign for the alignment
Foldseek supports to realign hits using TMalign as well as rescoring alignments using TMscore. 
```
foldseek easy-search example/d1asha_ example/ aln tmp --alignment-type 1
```
In case of the alignment type (`--alignment-type 1`) tmalign we sort the results by the TMscore normalized by query length. We write the TMscore into the e-value(=TMscore) as well as into the score(=TMscore*100) field.


### Rescore aligments using TMscore
Easiest way to get the alignment TMscore normalized by min(alnLen,qLen,targetLen) as well as a rotation matrix is through the following command:
```
foldseek easy-search example/ example/ aln tmp --format-output query,target,alntmscore,u,t
```

Alternative, it is possible to compute TMscores for the kind of alignment output (e.g. 3Di/AA) using the following commands: 
```
foldseek createdb example/ targetDB
foldseek createdb example/ queryDB
foldseek search queryDB targetDB aln tmpFolder -a
foldseek aln2tmscore queryDB targetDB aln aln_tmscore
foldseek createtsv queryDB targetDB aln_tmscore aln_tmscore.tsv
```

Output format `aln_tmscore.tsv`: query and target identifier, TMscore, translation(3) and rotation vector=(3x3)


### Search result visualisations
Foldseek can locally generate a search result HTML similiar to the [webserver](https://search.foldseek.com) by specifying the format mode `--format-mode 3`

```
foldseek easy-search example/d1asha_ example/ result.html tmp --format-mode 3
```

<p align="center"><img src="./.github/results.png" height="400"/></p>

### Cluster structures 
The following command aligns the input structures all-against-all and keeps only alignments with 80% of the sequence covered by the alignment (-c 0.8) (read more about alignment coverage [here](https://github.com/soedinglab/MMseqs2/wiki#how-to-set-the-right-alignment-coverage-to-cluster)). It then clusters the results using greedy set cover algorithm. The clustering mode can be adjusted using --cluster-mode, read more [here](https://github.com/soedinglab/MMseqs2/wiki#clustering-modes). The clustering output format is described [here](https://github.com/soedinglab/MMseqs2/wiki#cluster-tsv-format).

```
foldseek createdb example/ db
foldseek search db db aln tmpFolder -c 0.8 
foldseek clust db aln clu
foldseek createtsv db db clu clu.tsv
```

### Query centered multiple sequence alignment 
Foldseek can generate a3m based multiple sequence alignments using the following commands. 
a3m can be converted to fasta format using [reformat.pl](https://raw.githubusercontent.com/soedinglab/hh-suite/master/scripts/reformat.pl) (`reformat.pl in.a3m out.fas`).
```
foldseek createdb example/ targetDB
foldseek createdb example/ queryDB
foldseek search queryDB targetDB aln tmpFolder -a
foldseek result2msa queryDB targetDB aln msa --msa-format-mode 6
foldseek unpackdb msa msa_output --unpack-suffix a3m --unpack-name-mode 0
```

### Compile from source
Compiling `foldseek` from source has the advantage of system-specific optimizations, which should improve its performance. To compile it `git`, `g++` (4.9 or higher) and `cmake` (3.0 or higher) are required. Afterwards, the foldseek binary will be located in the `build/bin` directory.

    git clone https://github.com/steineggerlab/foldseek.git
    cd foldseek
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
    make -j
    make install
    export PATH=$(pwd)/foldseek/bin/:$PATH

:exclamation: If you want to compile `foldseek` on macOS, please install and use `gcc` from Homebrew. The default macOS `clang` compiler does not support OpenMP (by default) and `foldseek` will not be able to run multi-threaded. Adjust the `cmake` call above to:

    CC="$(brew --prefix)/bin/gcc-11" CXX="$(brew --prefix)/bin/g++-11" cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
