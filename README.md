Installation
------------

    git clone http://github.com/txje/fastat.git
    cd fastat
    mkdir -p incl
    cd incl
    git clone http://github.com/attractivechaos/klib.git
    cd ..
    make

Usage
-----

    fastat [options] FASTA/Q[.gz] ...
    
    Options:
      -v, --verbose: verbose
      -j, --joint: do aggregate stats over all files
      -c, --compact: print parseable compact stats
      -a, --all print size/GC stats for all sequences
      -h, --help: show this
