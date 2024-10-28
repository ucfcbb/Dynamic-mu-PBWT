# Dynamic mu-PBWT 

Dynamic mu-PBWT is a dynamic compressed PBWT. This tool supports dynamic updates (insertion/deletion)
of the index in its compressed format and supports long match query.

Related paper: "Dynamic mu-PBWT: Dynamic Run-length Compressed PBWT for Biobank Scale Data"

Contact Author: Pramesh Shakya (pramesh.shakya@ucf.edu)

## System Requirements
This tool has been tested on Linux environment. 

GCC Compiler: `>= g++11` 

## Dependencies
DYNAMIC library by Nicola Prezza and Alan Kuhnle (https://github.com/xxsds/DYNAMIC)

## Installation
Open a terminal and follow these steps:
1. Clone the github repository:
`git clone https://github.com/ucfcbb/Dynamic-mu-PBWT`

2. `cd` into the cloned folder:
`cd ./Dynamic-mu-PBWT`

4. Configure the build:
`./configure`
5. Build the project:
`./build`

There will be a `build` directory which contains 
the following executables `dmupbwt`, `insert` and `delete`.

Test the installation by running the executables as follows:
`./build/dmupbwt`,
`./build/insert`, and
`./build/delete`.
It will print out the usage details for each of the executables.


## Usage
### Construction of Dynamic mu-PBWT and Long match query
To construct Dynamic mu-PBWT on the input vcf or run long match query algorithm.
Following are the usage details when you run `./build/dmupbwt`:

|     Flag     |        Description         |                                       Details                                        |
|:------------:|:--------------------------:|:------------------------------------------------------------------------------------:|
| `-i <file>`  | Path to reference VCF file |                                uncompressed VCF file                                 |
| `-q <file>`  |   Path to query VCF file   |                                uncompressed VCF file                                 |
| `-o <file>`  | Path to output match file  | matches between haplotypes of query panel and the reference panel |
| `-l <value>` |  length threshold (sites)  |                    minimum length threshold for long match query                     |
| `-v ` |          verbose           |            prints out memory usage and other information about the panel             |

## Construction
To construct Dynamic mu-PBWT on the input vcf:
`./build/dmupbwt -i ./test_data/ref.1.vcf -v`

## Long match query
To find long match queries on Dynamic mu-PBWT :
`./build/dmupbwt -i ./test_data/ref.1.vcf -q ./test_data/query.1.vcf. -o ./test_data/out_matches -v`

## Insertion
To insert haplotypes on empty Dynamic mu-PBWT:

`./build/insert -i ./test_data/ref.1.vcf -v`

To insert haplotypes on non-empty Dynamic mu-PBWT:

`./build/insert -i ./test_data/ref.1.vcf -q ./test_data/query.1.vcf -v`


## Deletion
To randomly delete all the haplotypes from a given Dynamic mu-PBWT:

`./build/delete -i ./test_data/ref.1.vcf -v`

This builds the Dynamic mu-PBWT on the input VCF file and randomly deletes all the haplotypes.
