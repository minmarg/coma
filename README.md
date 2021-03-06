COMA software package for detection of distant evolutionary
relationships

(C)2010 Mindaugas Margelevicius
Institute of Biotechnology,Vilnius

# Available Platforms

   The COMA binaries are provided for the following platforms:

  *  Linux x86 (32 bit)
  *  Linux x64 (64 bit)
  *  Win32 (32 bit MS Windows OSs)

# Getting the COMA Software

   The package is available from:

   [https://github.com/minmarg/coma](https://github.com/minmarg/coma)

# Structure of the Package

   The main directories are shortly described below:

  *  binaries  --  contains the built and locally installed COMA package
     for Linux x86 (subdirectory, linux32), Linux x64 (linux64), and 
     Windows 32-bit OSs (win32). You can distribute the subdirectories
     to your prefered path and consider it as the installation path.

  *  src       --  is the main directory of the source files to be
     compiled and linked into the executables (see Basic Installation
     below).

  *  win32     --  contains the built sources for Windows platform.

  *  optimized --  the compiled and linked source files for Linux x64.
     It should be reconfigured and remade by make (see Basic 
     Installation below) if you want to built the sources in this
     directory.

# Basic Installation

   These are generic installation instructions.

   The `configure` shell script attempts to guess correct values for
various system-dependent variables used during compilation.  It uses
those values to create a `Makefile' in each directory of the package.
It may also create one or more `.h' files containing system-dependent
definitions.  Finally, it creates a shell script `config.status' that
you can run in the future to recreate the current configuration, a file
`config.cache' that saves the results of its tests to speed up
reconfiguring, and a file `config.log' containing compiler output
(useful mainly for debugging `configure').

   The file `configure.in` is used to create `configure` by a program
called `autoconf`.  You only need `configure.in` if you want to change
it or regenerate `configure` using a newer version of `autoconf`.

The simplest way to compile this package is:

  0. Sometimes it is useful to regenerate `configure` and other related
     files by running the command:

     make -f Makefile.cvs

  1. Make a new directory to keep the compiled and linked sources and
     `cd` to that directory:

     mkdir built; cd built

     Configure the package for your system giving your installation
     path where the package will be installed after compilation:

     ../configure --prefix=/your/installation/path/here

     The COMA package contains the MPI-dependent program `mpiscaler` for
     parallel calculations. The program is built optionally but
     requires the MPICH2 package installed in your system. If you have
     MPICH2 and want to built `mpiscaler`, configure the package as
     follows:

     ../configure --prefix=/your/installation/path/here  mpi=yes

     If you're using `csh` or an old version of System V, you might
     need to type `sh ../configure` instead to prevent `csh` from
     trying to execute `configure` itself.

     Running `configure` takes a while. While running, it prints some
     messages telling which features it is checking for.

  2. Type `make` to compile the package.

  3. Type `make install` to install the programs and any data files and
     documentation.

   For more information, see the INSTALL file.

# Getting Started

   There are three main programs sufficient for making of profiles,
profile databases, and running of the COMA method itself. These
programs: `makepro`, `makedb`, and `coma`, respectively, are in the
bin directory in your installation path.

   If you have a profile database (see Custom Databases below) and an
input multiple alignment file in the FASTA format (see below), the
simplest way to run coma is to type (names depend on your choice of 
naming):

     bin/coma -i myinput.fa -d mydb -o output

   It is convenient, however, to have an input multiple alignment
converted to the profile, especially in cases when coma is to be run
many times with that multiple alignment. To make a profile from the
multiple alignment, type:

     bin/makepro -i myinput.fa -o myinput.pro

   After the profile is made, `coma` may be run by using the profile:

     bin/coma -i myinput.pro -d mydb -o output

   If you want to mutually compare two profiles or multiple alignments,
you can indicate another profile or multiple alignment instead of
database name:

     bin/coma -i myinput.pro -d another.pro -o output

   To control the COMA search, you may want to set some options in file
to be passed to COMA:

     bin/coma -i myinput.pro -d another.pro -o output -p options

   By default the COMA options are read from file `var/options.txt`.

# Input Multiple Alignment

   The programs in the package processing input multiple alignment 
files (`coma`, `makepro`, etc.) recognize the FASTA format. Thus, input
multiple alignments should be prepared in this format.
   An example of a multiple sequence alignment in FASTA is shown below:

```
>d1qhka_ d.100.1.2 (A:) N-terminal domain of RNase HI...
GNFYAVRKGRE--T---G--------IYNTW---NECKNQVDGYG---GAIYKKFNSYEQAKSFLG
>gi|28379120|ref|NP_786012.1|:(2-47) ribonuclease H (putative)...
-KYYAVRKGRQ--P---G--------IYRTW---PETQKQVSGYP---QAQYKSFTSEKDAQDFMA
>gi|84386727|ref|ZP_00989753.1|:(2-47) hypothetical ribonuclease HI...
-KYYVVWKGRT--P---G--------IFTTW---NECKSQVDGFA---GARYKSFPTLGEAESAFG
>gi|116492108|ref|YP_803843.1|:(2-47) RNase H with double-stranded...
-KFYAVKKGRK--P---G--------LYLTW---DAAKQQVDGFA---GAVYKSFLTKAEAEEWMA
>gi|6323890|ref|NP_013961.1|:(1-47) Ribonuclease H1...
GNFYAVRKGRE--T---G--------IYNTW---NECKNQVDGYG---GAIYKKFNSYEQAKSFLG
```

   An Input multiple alignment from (PSI-)BLAST output can be obtained by 
running the following program:

     bin/blast2fa.pl -i myblast.aln -o myinput.fa

   where myblast.aln is a BLAST output file of pairwise
alignments, and resulting myinput.fa is the multiple alignment in
FASTA ready to be converted to profile:

     bin/makepro -i myinput.fa -o myinput.pro

# Custom Databases

   If you have a set of profiles constructed by `makepro` and want to
assemble them into a profile database to search against with `coma`,
use the `makedb` program:

     bin/makedb -o mydb  profil1.pro profile2.pro profile3.pro

   Alternatively, you can use wildcards to indicate a set of profiles:

     bin/makedb -o mydb  *.pro

   Or, if your profies are all in one directory, you may want to
point out the directory the profiles should be read from:

     bin/makedb -o mydb -d mydirectory

   It is IMPORTANT to note that the database will consist of several
files. The files will have the names exactly the same as given by
option `-o`. Thus, the parameter names used in `makedb` and `coma`
should be the same. It is convenient not to use extensions in these
names. The `coma` search against the database may be run by typing:

     bin/coma -i myinput.pro -d mydb -o output

# Comparison modes (Scoring schemes)

   There are two optional (option -S) comparison modes: profile and
global. The default profile comparison mode is to compare a pair of
profiles given reference statistical parameters. The other, global,
mode is to compare a query profile against the all profile
vectors from a database. The global mode is to be given reference
statistical parameters as well and thus, shoud NOT be confused with the
Global Score System which is constructed for composition-based
statistics to solve reference statistical parameters. The global
comparison mode and Global Score System are two different concepts.
Currently, the default profile comparison mode should be prefered to 
the global mode.

# References

M.Margelevicius, C.Venclovas (2010)
Detection of distant evolutionary relationships between protein families
using theory of sequence profile-profile comparisons.
BMC Bioinformatics 11:89.
http://www.biomedcentral.com/1471-2105/11/89

