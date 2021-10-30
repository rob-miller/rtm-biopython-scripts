# rtm-biopython-scripts
Python scripts to access or test Biopython code

* **bpbp-gist.py** - various test and profile code for exercising Biopython [internal_coords](https://github.com/rob-miller/biopython/blob/master/Bio/PDB/internal_coords.py)
code (latest with [PR #3774](https://github.com/biopython/biopython/pull/3774).  Run with `-h` for a help message, nice to configure the pdbDirs variable if you have a local repository of PDB files.  The name comes from 'Biopython buildprot(ein)',
and previously I was making it available as a gist.  Yes, the range of options and parameters is getting unwieldy.

* **ca-plot-7rsa.py** - displays a C-alpha distance plot for PDB file 7rsa chain A (ribonuclease A).  Requires matplotlib, code shows how to change to all backbone atom or all atom plots.

* **phe-pairs-2.py** - prints to stdout PDB format ATOM lines for aligned PHE residues from PDB files 3pbl and 3rze.  Demonstrates use of coordinate space transforms in Biopython internal_coords module to align specified atoms in PDB files and then inspect similarities/differences in environments of those atoms.  Suggest 'color by chain' for viewing.

* **test_all.py** - runs Biopython PDB or mmCIF parsers on file names or PDB idcodes specified on stdin, prints a line of 

    `file_no format idcode chain_id chain_count residues disordered_residues atoms disordered_atoms resolution missing_residues non-sequential_residues`

    for each chain.  As with **bpbp-gist.py**, run with `-h` for help info and useful to configure pdbDirs and cifDirs for idcode source data if you have a local repository.  
    
    General usage is something like `find <my-dir> | ./test_all.py` or `cat <my-list> | test_all.py>, or just start it and type idcodes at it.  End cleanly with ctrl-d or `q`,
    check the `-h` for more options.
