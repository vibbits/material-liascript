<!--
author:   Alexander Botzki
email:    Alexander.Botzki@vib.be
version:  0.1.1
language: en
narrator: US English Female

comment:  slides Protein Structure Analysis

logo: img/Logo.png

link:     https://cdnjs.cloudflare.com/ajax/libs/animate.css/3.7.2/animate.min.css
link:     https://raw.githubusercontent.com/vibbits/material-liascript/master/img/org.css

debug: true

-->

# Protein Structure Analysis

- Sequences, structures and databases
- Experimental methods (X-rays, electrons and NMR)
- Finding and visualising structures from the  Protein Data Bank
- Comparing structures
- Modelling mutations
- Creating homology models

## Sequences and Structures

!?[What is a protein](https://youtu.be/wvTv8TqWC48)<!-- width="890px" height="488px" -->

[^1](https://pdb101.rcsb.org/learn/videos/what-is-a-protein-video)

## Amino acids and peptide structure

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/amino-acids.png)<!-- width="100%" -->

!?[Proteins are polypeptides of aminoacids](https://www.youtube.com/watch?v=iBV67yamqlA&list=PLuIpgNT2hMwRQKFy4okoNQKiJwM8li3Sz&index=20)

## The Structure-Function Connection

- Folded proteins provide a well-defined 3D arrangement of functional groups, creating microenvironments and active sites.
- Structural changes are often involved in functional  mechanisms  (motor proteins, ...)

!?[](https://material.bits.vib.be/topics/protein-structure-analysis/images/hemoglobin.png)<!-- height="600px" -->

!?[](https://www.youtube.com/watch?v=E6LHoKe-zRU&list=PLuIpgNT2hMwRQKFy4okoNQKiJwM8li3Sz&index=13)
!?[](https://www.youtube.com/watch?v=CViLhxiPq1k&list=PLuIpgNT2hMwRQKFy4okoNQKiJwM8li3Sz&index=18)

## Databases

                           {{0-1}}
*********************

[**Uniprot**](https://www.uniprot.org/) approx. 570,000 curated sequences (SwissProt)

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/uniprot-logo.png)<!-- width="50%" -->

- mainly determined by large-scale DNA sequencing of individual genes or whole genomes

- increasingly automated annotation

*******************

                            {{1}}
*******************

[**Protein Data Bank**](http://www.wwpdb.org) approx. 208,000 experimentally determined structures

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/pdb-logo.png)<!-- width="50%" -->

- Mainly determined by X-ray crystallography, electron microscopy and recently lesser high-resolution NMR spectroscopy
- Automation is increasing, but there are still significant limitations in the rate of solving new structures

*******************

## X-Ray Crystallography - Technical setup


                     {{0-1}}
*********************
![](https://material.bits.vib.be/topics/protein-structure-analysis/images/xray-tech-setup.png)<!-- width="100%" -->

![](https://www.youtube.com/watch?v=4HZoDjJ4A2k Lecture 01, concept 18: Protein structure determination - X-ray crystallography & Cryo-EM)
*******************

                       --{{0}}--
In a crystal, a large number  of macromolecules are packed together in a regular grid, with consistent  orientations and relative distances.
When exposed  to an X-ray beam, this  arrangement gives rise to diffracted rays in specific directions, resulting in discrete spots on the planar detector. By rotating the crystal, a series  of images is obtained. From these, the intensities of all the diffracted rays of the crystal can be derived.

                      {{1}}
*******************


![](https://material.bits.vib.be/topics/protein-structure-analysis/images/diffraction-pattern.png)<!-- width="50%" -->

*******************

## X-Ray Crystallography - from diffraction pattern to atomic model

![](img/X-ray-crystallography-1.png)<!-- width="100%" -->

|Diffraction Spot Intensities and Phases| | Electron density|
|:---------------------------------------:|:--------:|:----------------:|
| $F_{obs}(h,k,l)$ and $\phi_{obs}(h,k,l)$ | $$R_{cryst} = \frac{\sum_{h,k,l}F_{obs}-F_{calc}}{\sum_{h,k,l}F_{obs}}$$ |  $\rho(x,y,z)$ |

## The Protein Databank

                                  {{0-1}}
*********************
[http://www.wwpdb.org](http://www.wwpdb.org)

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/wwpdb-welcome-page.png)

*******************

 {{1}}
*******************
- contains structures of  proteins, nucleic acids  and complexes,  determined by X-ray  crystallography, NMR  spectroscopy
- No purely theoretical  or ab initio models  (since 2006)
- Also stores supporting  experimental data
- Full deposition now  required by all peer-reviewed journals

*******************

## Exercise 1: Search the PDB

- Use the UniProt site to search for “dnak”.  
- Use the PDB site to search for “dnak”.
- Compare the UniProt and PDB result lists.
- Use the sequence search function to see if there are structures with sequences similar to that of the  DnaK C-terminal domain.
- Look at the summary pages of a number of  structures and note some interesting properties.
- Select a number of structures and create a report page.

### Search for a structure

#### Via [UniProt](http://www.uniprot.org/)
The way of searching for a specific protein structure depends on the data you already have. You might already have the PDB ID (a unique identifier), that's an easy one. But mostly you have the protein name or you just have a sequence. In the last cases I recommend to start from the UniProt website at <http://www.uniprot.org>, which is the best annotated protein database in the world. Our first model protein will be the molecular chaperone DnaK from *E. coli*. Below is an image of the UniProt search box where you can start your search for proteins.

![uniprotsearchbox.png](../../images/uniprotsearchbox.png "Search box")

<i class="fa fa-question-circle" aria-hidden="true"></i>

> #### {% icon hands_on %} Explore a PDB structure on the Uniprot web site
>
> 1. Go to the UniProt website and search for the DnaK protein
> - The UniProt search engine returns a list of DnaK protein sequences from a variety of organisms. An entry with accession code **P0A6Y8** and entry name **DNAK_ECOLI** should be near the top of this list.
> 2. Click on the *accession code* (column Entry) to view the protein page of this DnaK from the model organism *Escherichia coli*.
> 3. Click on *Structure* in the left-side menu and then look at the *3D structure databases* table.

####  Guidelines which PDB structures to select

Which structures (give the 4-character PDB ID) of the C-terminal domain of DnaK should preferentially be use for analysis and why?

> <details markdown="1">
>  
> <summary>Solution</summary>
>
>>
>>  Usually, the recommended selection criteria are using an X-ray structure with low resolution and low $R_{free}$ factor. Furthermore, the PDB database has pre-calculated a validation report for all of the structures.
>>
>>  As an example, have a look at https://www.ebi.ac.uk/pdbe/entry/pdb/4EZX under the section 'Experiments and Validation'. 
>>  For many PDB structures, there is also a re-done structure available with a vast amount of information on the quality of the X-ray structure and suggested 'better' models e.g. (https://pdb-redo.eu/db/4ezx). 
>>  In our case, we could opt for the structures 1DKX and 4EZX.
>>
>>  This is a difficult example since there are so many high resolution structures available. So, it is recommended to study the articles and compare the available structures to find your favorite structure for further analysis.
>> 
>
>
> </details>

### Via the Protein Data Bank by PDB ID

You can find structural information directly at the PDB database. The web site of the PDB consortium is located at <http://www.wwpdb.org>. This web site provides links to all members of the PDB (left side). It is a question of taste which resource you start off with. For X-ray structures, it is currently PDBe, RCSB PDB, PDBj. For NMR structres, you find the BMRB. In today's course, we focus on the PDB resources only.

Below is an image of the RCSB search box <http://www.rcsb.org/pdb/home/home.do> where you can start your search for structures.

![Pdbsearchbox_RCSB.png](../../images/pdbsearchbox_RCSB.png "PDB Search Box")

The PDB file with ID **1DKX** contains the atomic coordinates of the molecular chaperone (DnaK) from *E. coli*.

#### Search a structure on the RCSB web site

1. Go to the PDB website and type 1DKX in the search box


This will lead you to the same page we got earlier through UniProt.

#### Via the Protein Data Bank by sequence

In lots of cases we only have a sequence of which we would like to find out if there is structural information. The PDB can be searched using a sequence as input. Here is the sequence of the C-terminal substrate binding domain of DnaK:
```
    DVKDVLLLDVTPLSLGIETMGGVMTTLIAKNTTIPTKHSQVFSTAEDNQSAVTIHVLQGE
    RKRAADNKSLGQFNLDGINPAPRGMPQIEVTFDIDADGILHVSAKDKNSGKEQKITIKAS
    SGLNEDEIQKMVRDAEANAEADRKFEELVQTRNQGDHLLHSTRKQVEEAGDKLPADDKTA
    IESALTALETALKGEDKAAIEAKMQELAQVSQKLMEIAQQQHAQQQTAGADASANNAKDD
    DVVDAEFEEVKDKK
```
The PDB allows sequence searches through the same search box we used before.

![Pdbsearchbox_RCSB.png](../../images/pdbsearchbox_RCSB.png "PDB Search Box")

There is also an Advanced Search section, with a Blast/Fasta option in the Sequence Features section.

![Blastpdb.png](../../images/blastpdb.png "BLAST")

> #### Hands-on: BLAST search for PDB structure
>
> 1. Go to the Advanced Search section
> 2. Please select 'Sequence BLAST/PSI-BLAST' in the Query type drop down.
>    This method allows you to change some parameters for the search.
> 3. Copy and paste the sequence in the ''Sequence'' field
> 4. Press ''Submit query''.
> 5. You should see the same structures popping up as you saw in the UniProt page of DnaK.

# The PDB file

## Introduction

A PDB (Protein Data Bank) file is a plain text file that contains the
atom coordinates of a solved 3D structure of a protein or even DNA. Such
coordinate files can be obtained at the Protein Data Bank at
<http://www.rcsb.org/pdb>. Each PDB file has a unique identifier (ID)
consisting of 4 characters, the first one is always a number. Note: It
has been announced that the 4 character code will change in the future
<https://www.wwpdb.org/news/news?year=2017\#5910c8d8d3b1d333029d4ea8>.

The PDB file with ID **1DKX** contains the atomic coordinates of the
molecular chaperone (DnaK) from *E coli*.

> ### {% icon hands_on %} Hands-on: BLAST search for PDB structure
>
> 1. Go to the PDB website at <http://www.rcsb.org/pdb>
> 2. Type 1DKX in the search and try to answer the following questions.
{: .hands_on}

### {% icon question %} Questions

1. How many molecules were solved in this PDB file? What kind of molecules are these (proteins, peptides, DNA, ...)?
2. Does the structure represent the full protein? If not, how many residues are missing? Hint: Click on the UniProt KB link in the Sequence tab to see the full sequence.
3. Was this structure solved by X-Ray or NMR?
4. What is the atomic resolution and R-factor?

> <details markdown="1">
> <summary>{% icon solution %} Solution
> </summary>
> 1. Two, called polymers or chains: they are polypeptides ![Type](../../images/Mol_desc_1DKX.png)
> 2. To answer this question you can go to the sequence tab at the top:
>    - ![Uniprot view](../../images/Pdb_firstresiduesmissing_1dkx.png)
>    - Summary: a large chunk of the N-terminus is missing from the structure, the C-terminus is virtually complete.
> 3. X-RAY diffraction, as shown by Experimental Details
> 4. Atomic resolution: 2.00 Ångstrom and R-factor of 0.206
> </details>
{: .question}


## Downloading the structure

The file that holds the 3D coordinates can be downloaded by clicking on
*Download files* in the top right corner and then choosing *PDB file (text)*.
For convenience, save this file on your desktop. The filename is the
4-character unique PDB ID.

![Pdbdownloadfile1.png](../../images/pdbdownloadfile1.png)

> ### {% icon hands_on %} Hands-on: Open downloaded PDB file in an editor
> 1.   Open this file with a text editor, e.g. WordPad is an excellent tool for that.
> 2. Do you see the different sections in the PDB file? Analyse some ATOM lines and try to explain what kind of data is in each column.


## PDB File Format

                          {{0-1}}
******************

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/pdb-file-format.png)<!-- height="600px" -->

******************

                 --{{0}}--
The Protein Data Bank (PDB) format provides a standard representation for macromolecular structure data derived from X-ray diffraction and NMR studies. This representation was created in the 1970's and a large amount of software using it has been written.


                 {{1-2}}
************************************************

Atomic descriptions in classical PDB format

```
ATOM     21  N  ALEU A 392      64.898 118.959  40.827  0.60 24.19           N  
ATOM     22  N  BLEU A 392      64.981 118.958  40.769  0.45 24.73           N  
ATOM     23  CA ALEU A 392      64.818 117.251  41.235  0.60 22.53           C  
ATOM     24  CA BLEU A 392      64.908 117.527  41.044  0.45 23.55           C  
ATOM     25  C  ALEU A 392      63.999 117.260  40.013  0.60 21.73           C  
ATOM     26  C  BLEU A 392      64.972 116.794  39.706  0.45 22.95           C  
ATOM     27  O  ALEU A 392      64.789 117.304  38.508  0.60 19.91           O  
ATOM     28  O  BLEU A 392      64.270 117.164  38.760  0.45 21.80           O  
ATOM     29  CB ALEU A 392      62.994 117.229  41.239  0.60 22.01           C  
ATOM     30  CB BLEU A 392      63.588 117.170  41.739  0.45 23.26           C  
ATOM     31  CG ALEU A 392      62.020 116.866  42.907  0.60 22.43           C  
ATOM     32  CG BLEU A 392      63.411 117.414  43.240  0.45 24.38           C  
ATOM     33  CD1ALEU A 392      63.906 117.128  42.987  0.60 22.36           C  
ATOM     34  CD1BLEU A 392      62.004 117.019  43.655  0.45 24.55           C  
ATOM     35  CD2ALEU A 392      62.660 117.298  44.405  0.60 22.73           C  
ATOM     36  CD2BLEU A 392      64.432 116.611  44.023  0.45 24.44           C  
ATOM     37  N  AASP A 393      65.422 116.659  39.809  0.60 21.96           N  
ATOM     38  N  BASP A 393      65.827 115.776  39.625  0.45 22.93           N  
ATOM     39  CA AASP A 393      65.972 115.093  38.433  0.60 22.44           C  
ATOM     40  CA BASP A 393      65.967 114.977  38.409  0.45 23.84           C  
ATOM     41  C  AASP A 393      64.830 113.969  38.617  0.60 21.62           C  
ATOM     42  C  BASP A 393      64.900 113.889  38.543  0.45 22.45           C  
ATOM     43  O  AASP A 393      64.810 113.091  39.785  0.60 21.11           O  
ATOM     44  O  BASP A 393      64.962 113.054  39.456  0.45 21.89           O  
ATOM     45  CB AASP A 393      66.589 114.675  39.315  0.60 24.09           C  
ATOM     46  CB BASP A 393      67.374 114.360  38.335  0.45 26.75           C  
ATOM     47  CG AASP A 393      66.965 114.032  38.106  0.60 25.18           C  
ATOM     48  CG BASP A 393      67.888 114.196  36.899  0.45 29.40           C  
ATOM     49  OD1AASP A 393      66.303 115.375  36.780  0.60 25.41           O  
ATOM     50  OD1BASP A 393      67.190 114.588  35.938  0.45 31.28           O  
ATOM     51  OD2AASP A 393      68.508 115.272  37.824  0.60 25.36           O  
ATOM     52  OD2BASP A 393      69.014 113.671  36.724  0.45 32.22           O  
```



************************************************

               --{{1}}--
The original short term goal of the working group was to fulfill the mandate set by the IUCr: to define mmCIF data names that needed to be included in the CIF dictionary in order to adequately describe the macromolecular crystallographic experiment and its results. In January 1997, the mmCIF dictionary was completed and submitted to COMCIFS for review, and version 1.0 was released in June 1997.

                 {{2-3}}
************************************************

Atomic descriptions in mmCIF format [^1](http://mmcif.wwpdb.org/docs/pdb_to_pdbx_correspondences.html)

```
ATOM   21   N N   A LEU A 1 4   ? 64.898 118.959 40.827 0.60 24.19 ? 392  LEU A N   1
ATOM   22   N N   B LEU A 1 4   ? 64.981 118.958 40.769 0.45 24.73 ? 392  LEU A N   1
ATOM   23   C CA  A LEU A 1 4   ? 64.818 117.251 41.235 0.60 22.53 ? 392  LEU A CA  1
ATOM   24   C CA  B LEU A 1 4   ? 64.908 117.527 41.044 0.45 23.55 ? 392  LEU A CA  1
ATOM   25   C C   A LEU A 1 4   ? 63.999 117.260 40.013 0.60 21.73 ? 392  LEU A C   1
ATOM   26   C C   B LEU A 1 4   ? 64.972 116.794 39.706 0.45 22.95 ? 392  LEU A C   1
ATOM   27   O O   A LEU A 1 4   ? 64.789 117.304 38.508 0.60 19.91 ? 392  LEU A O   1
ATOM   28   O O   B LEU A 1 4   ? 64.270 117.164 38.760 0.45 21.80 ? 392  LEU A O   1
ATOM   29   C CB  A LEU A 1 4   ? 62.994 117.229 41.239 0.60 22.01 ? 392  LEU A CB  1
ATOM   30   C CB  B LEU A 1 4   ? 63.588 117.170 41.739 0.45 23.26 ? 392  LEU A CB  1
ATOM   31   C CG  A LEU A 1 4   ? 62.020 116.866 42.907 0.60 22.43 ? 392  LEU A CG  1
ATOM   32   C CG  B LEU A 1 4   ? 63.411 117.414 43.240 0.45 24.38 ? 392  LEU A CG  1
ATOM   33   C CD1 A LEU A 1 4   ? 63.906 117.128 42.987 0.60 22.36 ? 392  LEU A CD1 1
ATOM   34   C CD1 B LEU A 1 4   ? 62.004 117.019 43.655 0.45 24.55 ? 392  LEU A CD1 1
ATOM   35   C CD2 A LEU A 1 4   ? 62.660 117.298 44.405 0.60 22.73 ? 392  LEU A CD2 1
ATOM   36   C CD2 B LEU A 1 4   ? 64.432 116.611 44.023 0.45 24.44 ? 392  LEU A CD2 1
ATOM   37   N N   A ASP A 1 5   ? 65.422 116.659 39.809 0.60 21.96 ? 393  ASP A N   1
ATOM   38   N N   B ASP A 1 5   ? 65.827 115.776 39.625 0.45 22.93 ? 393  ASP A N   1
ATOM   39   C CA  A ASP A 1 5   ? 65.972 115.093 38.433 0.60 22.44 ? 393  ASP A CA  1
ATOM   40   C CA  B ASP A 1 5   ? 65.967 114.977 38.409 0.45 23.84 ? 393  ASP A CA  1
ATOM   41   C C   A ASP A 1 5   ? 64.830 113.969 38.617 0.60 21.62 ? 393  ASP A C   1
ATOM   42   C C   B ASP A 1 5   ? 64.900 113.889 38.543 0.45 22.45 ? 393  ASP A C   1
ATOM   43   O O   A ASP A 1 5   ? 64.810 113.091 39.785 0.60 21.11 ? 393  ASP A O   1
ATOM   44   O O   B ASP A 1 5   ? 64.962 113.054 39.456 0.45 21.89 ? 393  ASP A O   1
ATOM   45   C CB  A ASP A 1 5   ? 66.589 114.675 39.315 0.60 24.09 ? 393  ASP A CB  1
ATOM   46   C CB  B ASP A 1 5   ? 67.374 114.360 38.335 0.45 26.75 ? 393  ASP A CB  1
ATOM   47   C CG  A ASP A 1 5   ? 66.965 114.032 38.106 0.60 25.18 ? 393  ASP A CG  1
ATOM   48   C CG  B ASP A 1 5   ? 67.888 114.196 36.899 0.45 29.40 ? 393  ASP A CG  1
ATOM   49   O OD1 A ASP A 1 5   ? 66.303 115.375 36.780 0.60 25.41 ? 393  ASP A OD1 1
ATOM   50   O OD1 B ASP A 1 5   ? 67.190 114.588 35.938 0.45 31.28 ? 393  ASP A OD1 1
ATOM   51   O OD2 A ASP A 1 5   ? 68.508 115.272 37.824 0.60 25.36 ? 393  ASP A OD2 1
```

************************************************

                 {{3}}
************************************************

mmCIF file looks like a paired collection of data item names and values.

```
data_1DKX
#
_entry.id   1DKX
#
_audit_conform.dict_name       mmcif_pdbx.dic
_audit_conform.dict_version    5.279
_audit_conform.dict_location   http://mmcif.pdb.org/dictionaries/ascii/mmcif_pdbx.dic
#
loop_
_database_2.database_id
_database_2.database_code
PDB   1DKX         
WWPDB D_1000172826
#
_pdbx_database_status.status_code                     REL
_pdbx_database_status.entry_id                        1DKX
_pdbx_database_status.recvd_initial_deposition_date   1996-06-03
_pdbx_database_status.deposit_site                    ?
_pdbx_database_status.process_site                    ?
_pdbx_database_status.SG_entry                        .
_pdbx_database_status.pdb_format_compatible           Y
_pdbx_database_status.status_code_mr                  ?
_pdbx_database_status.status_code_sf                  ?
_pdbx_database_status.status_code_cs                  ?
#
loop_
_audit_author.name
_audit_author.pdbx_ordinal
'Zhu, X.'           1
'Zhao, X.'          2
'Burkholder, W.F.'  3
'Gragerov, A.'      4
'Ogata, C.M.'       5
'Gottesman, M.E.'   6
'Hendrickson, W.A.' 7
#
_citation.id                        primary
_citation.title                     'Structural analysis of substrate binding by the molecular chaperone DnaK.'
_citation.journal_abbrev            Science
_citation.journal_volume            272
_citation.page_first                1606
_citation.page_last                 1614
_citation.year                      1996
_citation.journal_id_ASTM           SCIEAS
_citation.country                   US
_citation.journal_id_ISSN           0036-8075
_citation.journal_id_CSD            0038
_citation.book_publisher            ?
_citation.pdbx_database_id_PubMed   8658133
_citation.pdbx_database_id_DOI      ?
#
```

************************************************

## Occupancy

                 {{0}}
************************************************

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/occupancy.png)<!-- width="80%" -->
[^1](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/dealing-with-coordinates)

************************************************

                 --{{0}}--
Macromolecular crystals are composed of many individual molecules packed into a symmetrical arrangement. In some crystals, there are slight differences between each of these molecules. For instance, a sidechain on the surface may wag back and forth between several conformations, or a substrate may bind in two orientations in an active site, or a metal ion may be bound to only a few of the molecules. When researchers build the atomic model of these portions, they can use the occupancy to estimate the amount of each conformation that is observed in the crystal. For most atoms, the occupancy is given a value of 1, indicating that the atom is found in all of the molecules in the same place in the crystal. However, if a metal ion binds to only half of the molecules in the crystal, the researcher will see a weak image of the ion in the electron density map, and can assign an occupancy of 0.5 in the PDB structure file for this atom. Occupancies are also commonly used to identify side chains or ligands that are observed in multiple conformations. The occupancy value is used to indicate the fraction of molecules that have each of the conformations. Two (or more) atom records are included for each atom, with occupancies like 0.5 and 0.5, or 0.4 and 0.6, or other fractional occupancies that sum to a total of 1.

## Related Web sites

- Nucleic Acid Database: DNA and RNA structures

  [http://ndbserver.rutgers.edu/](http://ndbserver.rutgers.edu/)

- PDB-REDO: automatically re-refined deposited structures, using the latest methods

  [http://www.cmbi.ru.nl/pdb_redo/](http://www.cmbi.ru.nl/pdb_redo)

- EBI: many online tools for structure analysis

  [http://www.ebi.ac.uk/Tools/structure/](http://www.ebi.ac.uk/Tools/structure/)

- Replaced Electron Density Server: convenient overview of quality parameters for crystal structures

  [http://www.ebi.ac.uk/pdbe/litemol](http://www.ebi.ac.uk/pdbe/litemol)

## High-Resolution NMR Spectrometry - Principles


                 {{0}}
************************************************

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/nmr-peaks-to-structure.png)<!-- height="307px"-->

************************************************

                              --{{0}}--
Many atomic nuclei, including the ubiquitous hydrogen nuclei,  resonate at specific radio frequencies when placed in a  strong, uniform magnetic field. The chemical environment of each individual atom slightly modulates its exact resonance  frequency. In macromolecules with thousands of atoms, many different  effects combine to generate an extremely complicated pattern of chemical shifts, which therefore more or less uniquely  identify each atom. Multidimensional spectra allow these  frequencies to be assigned to specific atoms.


## High-Resolution NMR Spectroscopy - Nuclear Overhauser Effect

                 {{0}}
************************************************

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/nmr-noe.jpg)<!-- height="600px" -->

************************************************

                          --{{0}}--
When two atoms are near each other in 3D space, they can exchange magnetization, giving rise to crosspeaks at the  intersection of their respective frequencies. This nuclear Overhauser effect (NOE) is used to effectively measure the distances between pairs of atoms, at least qualitatively.

## High-Resolution NMR Spectroscopy

                 {{0}}
************************************************

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/nmr-model-example.png)<!-- height="600px" -->

************************************************

                --{{0}}--
After identification of the atoms by means of their unique chemical shifts, distance restraints are derived from the Overhauser crosspeaks. An extended model of the protein is generated, and then condensed into a shape that is consistent with as many of distance restraints as possible.

## Other methods

- Electron microscopy (and especially Cryo-electron microscopy): Electron crystallography and single particle reconstruction
- Small-angle X-ray and neutron scattering (SAXS and SANS)

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/saxs.png)<!-- height="600px" -->

## Related Web sites

- BioMagResBank: experimental data for NMR-  derived structures (lists of distance restraints  and other experimentally derived properties)

  [http://www.bmrb.wisc.edu/](http://www.bmrb.wisc.edu/)

- BioIsis: database of SAXS-derived structures

  [http://www.bioisis.net](http://www.bioisis.net)

- EMBL database of SAXS-derived structures

  [http://www.sasbdb.org](http://www.sasbdb.org)

- EM Databank for cryo-EM structures

  [http://www.emdatabank.org](http://www.emdatabank.org.left)

## Assessing Structure Quality

General geometric properties (bond lengths and angles, Ramachandran distribution, …):
[MolProbity](http://molprobity.biochem.duke.edu/)

**Crystal Structures**

- Diffraction data resolution and completeness (PDB)
- Final $ R_{cryst} $ and $ R_{free} $ factors (PDB)
- Local correspondence to electron density [EDS](https://www.ebi.ac.uk/pdbe/eds)

**NMR Structures**

- Number of distance restraints and additional experimental data sources [BMRB](http://www.bmrb.wisc.edu/)
- Local restraint density [on-line NMR constraint analyser](http://molsim.sci.univr.it/bioinfo/tools/constraint/index.html)

**Other techniques**

- Difficult to generalise: carefully read and consider the  published protocol

## Molecular Graphics Software

- [PyMOL](http://www.pymol.org/): high-quality output,  good examples on [Wiki](http://www.pymolwiki.org/)
- [ChimeraX](https://www.cgl.ucsf.edu/chimerax/): good  documentation on website
- [VMD](https://www.ks.uiuc.edu/Research/vmd/): excellent for the analysis of MD trajectories
- [Yasara](http://www.yasara.org)  
- [SwissPDBViewer](http://spdbv.vital-it.ch/)

## ChimeraX

- ChimeraX is the next-generation molecular visualization program from the [UCSF RBVI](https://www.rbvi.ucsf.edu/). 

- Provides high-performance rendering of large structures and data

- Calulate fast and robust solvent-excluded surfaces

- Integration with AlphaFold and ESMFold AI-predicted protein structures (fetch existing or calculate new) 

## Exercise 2: Show a structure

[Link to exercise 2](https://material.bits.vib.be/topics/protein-structure-analysis/tutorials/visualise-structures/tutorial.html)

- Load PDB entry 1TRZ using the File menu.

- Compare the default representations and use the various controls to change the view  of the protein.

- Explore the Edit and View menus to change  various aspects of the graphical representation  of the structure.

- Examine the hydrogen bonding patterns. Display a molecular surface.

- Create an interesting view of the molecule and save it as a high-resolution image.

### Working with ChimeraX

There are great videos on working with ChimeraX that you can access via [a playlist from Brown Lab](https://www.youtube.com/playlist?list=PL4eF1KHNgDfIYSKCS3_S0PTRYtYTV9Myi).

We suggest to start with the first and the second video of this playlist.

### Scene styles

Open the PDB with code 1TRZ in ChimeraX via pasting `open 1TRZ` into the command field at the bottom. This will fetch the structure from the RCSB Protein Data Bank and load it into ChimeraX. Please have a look at the log window where more information about the structure is shown next to the full command `open 1TRZ format mmcif fromDatabase pdb`.

The default view is a quite nice one highlighting the interactions of the side chains of three residues with zinc and sodium atoms as well as selected water molecules.

Changing the molecule display can be easily achieved by toggling the icons *Show* and *Hide* from the sections **Atoms**, **Cartoons** and **Styles** in the `Home` menu.

### Showing and hiding residues

If you click the *Hide* icon followed by the *Show* icon of the **Atoms** section, all atoms and residues are shown by default. Toggling the Cartoons to state *Show* doen not explicitly show atoms and residues but are merely a impressionistic representation of the structure. The menu command **Actions** --> **Atoms/Bonds** --> **Show Sidechain/Base**, to a certain extent, shows atoms, but only of side chains, not main chain atoms. Mostly to do structure analysis, we want to show only the most interesting residues, the ones we want to analyze, and hide all the others.

The structure of insulin was crystallized together with some water molecules. In many cases, it is no problem to permanently delete those waters via `Select > Residues > HOH` followed by `Actions > Delete`. To visualize the water molecules, got to `Select > Residues > HOH` followed by `Atoms/Bonds > Show`. Do you see the red water (oxygen) atoms floating around the surface?

There are several ways to show the residues of interest:

  From the menu: View > Show atoms in > Residue. Select Cys7 from Molecule A and press OK.

  From the sequence selector. Hover the mouse on the bottom of the screen, you will see the sequence selector opening. Open it permanently by pressing the blue nailpin on the left side of it. Search for Cys7 from Molecule B, right-click and select Show > Residue.

Now show the atoms of His5 in Molecule B using a method of choice.

And now that we’re on it, what is special about the two cysteines we just visualized?

Hiding individual atoms or residues works in the same way as showing them, only now you should go to Hide atoms in the menu.

### Showing and Hiding Secondary Structure

Most published molecular images show a detailed active site and all the rest is hidden for clarity. From the previous exercise we show the atoms of 3 residues (let’s assume this is our active site). Now secondary structure of the rest of the molecule is also still visible. To hide all that, we do not have to hide atoms, but hide the secondary structure (the F5 tube view) from the rest of the structure. Atoms and residues in YASARA are not the same as the term ‘secondary structure’. Atoms and residues are balls and sticks, ‘secondary structure’ is an artistic impression of the structure (beta sheet arrows, helix ribbons, …). If you get this concept, you are a YASARA master.

So let’s hide many of the secondary structure, but keep just a few stretches around our active site. Our active site is Cys7 (A), Cys7 (B) and His 5 (B). This can be done in several ways. Since we would have to hide almost everything, I propose to hide first everything and then show again those stretches that we want. But if you have a better idea, I would like to hear it.

Hide all secondary structure via View > Hide secondary structure of > All.

Then show stretches of residues 2-10 in Mol B and residues 4-10 in Mol A in tube view via View > Show secondary structure > Tube through > Residue.

Then select the correct stretches of residues by keeping the CTRL key pressed to select multiple residues.

There are still some metal-bound histidines flying around that weren’t hidden because they are metal bound (a YASARA specific thing). Hide those histidines by clicking on one of the sidechain atoms, then
right-click and select Hide atoms > Residue.

The nasty dative bonds and metals can be removed simply by deleting all of them via Edit > Delete > Residue > Name.

In the name column select all the metals and ions you can find.

Et voilà, a publication ready image!

### Labels

You can put labels on the residues you want to highlight by selecting an atom from a residue (right-click). Subsequently, you select `Label > Residue` and choose the formatting of the label.

Via the same menu, you can change the height to 0.7 A or e.g. 20 pixels.

### Colors

You can color on all levels: atoms, residues, molecules and objects. So be careful, if you color a residue, all of its atoms will get that color. If you color a molecule, all atoms in that molecule will get that color.

Let’s color the secondary structure (the backbone in our case) of our active site in orange. But the sidechains should keep their Element colors. So we shouldn’t color entire residues, but only a selected atom set. Therefore our selection will be at the atom level, not the residue level. Go to `Actions > Color > All Options` and unselect the objects to which coloring applies.

Then select the orange color and select `Apply`. Hopefully, it is a satisfying result.

### Saving the work 

It would be a pitty that you spent hours creating fancy molecular graphics for that next paper while you can’t continue on the work the next day. That’s why ChimeraX can save the entire session including orientations, colors, views, everything. To save the current session, go to File > Save and select `ChimeraX session` as File type.

Choose a filename such as `MyInsulin.cxs`. To load the work again in ChimeraX go to `File > Open`.

### Save as high-quality image

To save the current view to a high quality publication ready image file, click on the Snapshot icon which will save a PNG file at the default location. 

### Distances

Distances between atoms are calculated as follows:

    Select the first atom
    Keep CTRL pressed and select the second atom.
    Left of the screen indicates the Marked Distance in Angstrom.

### Hydrogen bonds

To show hydrogen bonds, you can go to Analysis section in the Molecule Display tab. Before you click on the `H-bonds` icon, it is recommended to define a search space by selecting residues or molecules between which you would like to calculate hydrogen bonds. 
    
To hide or remove the hydrogen bonds, you have multiple choices: use `Molecule Display > Hide H-bonds` or select the hydrogen bonds object in the Models overview on the right side of the ChimeraX window and click on *Close*. 

### Surfaces

It can be very useful and informative to show the molecular surface of a protein. you can visualize cavities, ligand binding sites, etc … To show the molecular surface of one monomer of dimeric insulin, go to `Molecule Display > Surfaces > Show`.

In our case, 4 surfaces will be calculated which display status and color can be individually changed using the toggles in the Models window.

### Another exercise

represenation of the insulin structure – pdb code: 1TRZ

Hints:

    Choose the proper secondary structure scene style (F6 was used here);
    Find the correct orientation first;
    Color all backbone atoms in gray;
    Find the residue numbers of the 2 colored helices;
    Color those residues magenta;
    Show the sidechain atoms and the CA of the two histidines and the
    glutamate;
    Color the sidechain atoms of all residues in the Element color;
    Label the histidines and the glutamate;
    If you need some help how to change the parameters for the label,
    please have a look at Help > User guide and search for the label command.

### More coloring

Download GroEL via PDB code 1WE3 in YASARA.

Try to reproduce (approximately) the following image (hints below):

Hints:

> Load the PDB via `open 1we3` 
> Zoom out and find the correct orientation
> Delete the ADP, DMS and Mg molecules (are treated as residues in YASARA). So Edit > Delete > Residue > Adp …
> Color by molecule (every molecule will get another color) and color by gradient (now you need to specify 2 colors, the begin and end color).
> Choose a first color (eg. color with code 0)
> Choose a second color (eg. color with code 300, so you go over the entire color wheel spectrum)

### Conclusion

Now, you have explored the ChimeraX interface and acquainted with basic visualisations. You have identified how you can visualise secondary structure elements, surfaces, and hydrogen bonds. And most importantly, you can create publication-ready figures using ChimeraX.

## Protein folds are the structures of domains

- Similarities in assembly of secondary structure elements

- So not based on sequence like motifs but on 3D structure

- Folds represent the shapes of protein domains

!?[https://youtu.be/nJVkwdNm_HY](https://youtu.be/nJVkwdNm_HY) !?[https://youtu.be/1ljeUscKE9Y](https://youtu.be/1ljeUscKE9Y) !?[https://youtu.be/wWHtcbAnOmQ](https://youtu.be/wWHtcbAnOmQ)


* SCOP (http://scop.mrc-lmb.cam.ac.uk/scop/)

  Class, Fold, Superfamily, Family
  The human hemoglobin protein (1HHO-A), is a domain classfied by SCOP as the following.
  Class: all helical proteins; fold: globin-like; superfamily: globin-like; and family: globins.
  
!?[https://youtu.be/Wz6Ugn7BZjg](https://youtu.be/Wz6Ugn7BZjg)

* CATH (http://www.cathdb.info/)
  class (equivalent to SCOP class), architecture, topology (equivalent to SCOP fold), and homologous superfamily (equivalent to SCOP superfamily).
  1HHO-A is a domain classied by CATH as the following.
  Class: mainly helical; architecture: orthogonal bundle; topology: globin-like; and superfamily: globins.

!?[https://youtu.be/mb898MT3eLc](https://youtu.be/mb898MT3eLc)

## Similarity searches based on 3D structure

Similarity on structural level: aligning 3D structures

Structure of query protein is known and aligned to PDB structures

- [Foldseek](https://search.foldseek.com/search) - preferred method
- [DALI](http://ekhidna.biocenter.helsinki.fi/dali_server/)
- [VAST+](https://www.ncbi.nlm.nih.gov/Structure/vastplus/vastplus.cgi)

!?[Presentation by Martin Steinegger about FoldSeek](https://www.youtube.com/watch?v=k5Rbi22TtOA)

Compare proteins with low sequence similarity: similar structure implies homology -> same function

## Exercise 3: Use FoldSeek to find similar structures from distinct organisms

[Link to exercise 3](https://material.bits.vib.be/topics/protein-structure-analysis/tutorials/compare-structures/tutorial.html)

- Compare directed search via BlastP of protein A

- Download the five provided PDB files and open them in Yasara.

- Use the `Analyze|Align|Objects` with MUSTANG function to align the four last objects with the first one.

- Use the space bar to open the text console and  see the reported root mean square deviations as well as the number of residues matched.

$$ rmsd = \sqrt{\frac{1}{N}\sum_{i=1}^{N}R_{i}^{2}} $$

- Color all structures by B factor and compare the  distribution to the local variability of the structures.

[Additional information](https://web.stanford.edu/class/cs279/lectures/lecture5.pdf)

## Structure Alignment

Structures can be structurally aligned to achieve the best possible match between corresponding atoms. It is  possible to consider all atoms, or only a subset  (such as the Cα atoms).

When the structures are very similar, determining  which atoms should be matched to each other is trivial. When there are larger differences, it takes more preparation.

Different algorithms use different combinations of sequence- and secondary of tertiary structure-based information.

- [DALI](http://ekhidna.biocenter.helsinki.fi/dali_server/) [^1](Holm and Sander, 1993 and Holm 2019 https://doi.org/10.1093/bioinformatics/btz536)
- [TM-Align](https://zhanglab.ccmb.med.umich.edu/TM-align/)
- [DeepAlign](https://www.nature.com/articles/srep01448)
- [LGA](https://dx.doi.org/10.1093%2Fnar%2Fgkg571)
- [CE](https://doi.org/10.1093/protein/11.9.739)
- [FatCat](https://doi.org/10.1093/bioinformatics/btg1086)
- [MUSTANG](http://lcb.infotech.monash.edu.au/mustang/mustang_psfb-final.pdf)
- [MMLigner](https://doi.org/10.1093/bioinformatics/btw757)

## Exercise 4: Compare Structures

- Download the five provided PDB files and open them in ChimeraX.

- Use the `Analyze|Align|Objects` with MUSTANG function to align the four last objects with the first one.

- Use the space bar to open the text console and  see the reported root mean square deviations as well as the number of residues matched.

$$ rmsd = \sqrt{\frac{1}{N}\sum_{i=1}^{N}R_{i}^{2}} $$

- Color all structures by B factor and compare the  distribution to the local variability of the structures.

### Structural comparison and RMSD 
We compare structures by structurally aligning them on top of each other. That is, we
align structurally equivalent atoms. For now, we will only use CA atoms as a representation of the backbones. 
But Yasara also can align on any type of atom you want. You always need to specify:

-  source object(s): the structure(s) that needs to be rotated and translated to superpose on anoth
er structure
-  target object: the structure to superpose on

An optimal alignment is found when the root-mean-square deviation (RMSD) is at a minimum. 
The RMSD is given as:
![RMSD](../../images/RMSD.gif "calculation of RMSD")
where R is the distance between two structurally equivalent atom pairs (CA in our case) and n is the total number of atom pairs.

> #### {% icon hands_on %} Hands-on: Data download
>
> 1. Download the following adapted PDB files from [Zenodo](https://zenodo.org/record/3550492#.XdeNL1dKiUk) 
>
>    ```
>     1DKX_1.pdb 1DKY_1.pdb 1DKZ_1.pdb 3DPO_1.pdb 3DPP_1.pdb 
>    ```
>

### Aligning multiple structures using YASARA

Now load all of them in YASARA:

```
File > Load > PDB File
```

and select the CA (C-alpha) view (F4) and superpose with the MUSTANG algorithm:


```
Analyze > Align > Objects with MUSTANG
```


In the first window you have to select the source objects that will be repositioned. Select Objects 2 till 5. In the second window you select the target Object to superpose on. That would then be the first object.

Notice that YASARA prints the RMSD of every structural alignment in the lower Console. Open the Console by pressing the spacebar once or twice to extend it.

Color the atoms by their B-factor:

```
View > Color > Atom > Belongs to or has > All
Then choose BFactor in the next window and press 'Apply unique color'.
```

High BFactors are yellow, low BFactors are blue.

#### Questions

Do you see a correlation between the BFactors and the variability in the structure?

> <details markdown="1">
> <summary>{% icon solution %} Solution
> </summary>
> 
> ![Structural alignemnt](../../images/aligned-structures.png) 
>
> 
> </details>


## PDB File Surprises

- Missing atoms, residues and loops  
- Multiple molecules in the asymmetric unit
- Incomplete oligomers due to coincident crystal and  oligomer symmetry
- Cα-only models
- alternate conformations
- multiple models, especially for NMR ensembles
- use of B factors to represent other properties  
- other non-standard extensions (PDBQT, ...)

## Exercise 5: Model a Mutation

- Load the 2AC0.sce Yasara scene file.
- Set an appropriate structure representation.
- Locate residue Ala159 using the sequence view, and right-click to access the `FoldX|Mutate` residue function. Change it to a Trp residue.
- Look at the effect of the substitution on the structure, and use the space bar to open the text  console and read the output of the FoldX calculation.
- Mutate Arg273 to an alanine side chain. Discuss  the effects of this substitution.

## Exercise 6: Study Protein-Ligand Interactions

### Introduction
{:.no_toc}

<!-- This is a comment. -->

The goal of this exercise is appreciate how protein interactions can be studied through visual inspection and other software tools. Protein interactions can be classified into different groups regarding the molecular properties and functions of the interacting partners. (These groups are intertwined in several cases.) Some examples include:

- the interactions of proteins with other proteins, small molecules, carbohydrates, lipids or nucleic acids;
- Receptor-ligand interactions;
- Antigen-antibody interactions;
- Enzymatic interactions, enzyme-inhibitor interactions.

#### Exploring the structure of a nanobody-stabilized active state of the β2 adrenoceptor - the ligand 

We will start with exploring one crystal structure of the β2 adrenoceptor. Together with the Steyaert lab from VIB, Kobilka published several crystal structures of the β2 adrenoceptor in its various activation states (Rasmussen et al. Nature 2011, 477)


> #### {% icon hands_on %} Get the structure
>
> 1. Download the crystal structure 3P0G from the PDB into YASARA. 
>
>    ```
>    File - Load - PDB file from internet    
>    ```
>    As you can immediately appreciate, it is a bigger crystal structure with more than one molecule. 
>

> #### {% icon question %} Questions
>
> 1. How many molecules are present in the crystallized structures? 
> 2. And how many chain identifiers are used? 
>
> <details markdown="1">
> <summary>{% icon solution %} Solution
> </summary>
> 
> 1. There are three molecules, chain A Beta-2 adrenergic receptor; Endolysin, chain B Camelid Antibody Fragment, and a small molecule ligand. 
     Also have a look at PDBe [3P0G](https://www.ebi.ac.uk/pdbe/entry/pdb/3p0g) which gives a very nice overview of the structure and its composition.
> 2. Only two chain identifiers A and B. Sometimes, this leads to issues depending on the software you might want to use for downstream processing.
> 
> </details>
>

Some software routines need seperate chain identifiers for molecular entities to work correctly, so we suggest to rename the small molecule to chain L.

> #### {% icon hands_on %}  
>
> 1. Activate the Head-up display
> 2. Select Rename
> 3. Enter 'L' to proceed with the renaming. 
>

We first have a look whether we can find out if there are specific interactions of the small molecule ligand with the adrenoreceptor.

In order to do so, we first have to add Hydrogens to all present molecules.

> #### {% icon hands_on %}  
>
> 1. Edit - Add - hydrogens to : All 
> 2. Change the display of the ligand to Sticks
> 3. Select the amino acids of the binding pocket i.e. a sphere of 10 Angstrom around the ligand
>    ```
>    Select – in sphere around – Residue and drag with the mouse until the display says 10 Å
>    ``` 
> 4. ```
>    View – show interactions – hydrogen bonds of - Residues
>    ```
>    select 'Selected' in the panel Belongs to or has
>    and press OK in the subsequent window
>

Given that hydrogen bonding is dependent on the definition of a hydrogen bond in the program, it is not a bad idea to use other tools to compare the analysis. There are many options to do this online if you look at published crystal structures. Next to the tools which are directly linked out from the web site of the crystal structure at the PDB database you can use the [ProteinPlus server](http://proteinsplus.zbh.uni-hamburg.de/)

Go to the web site of ProteinPlus and enter the PDB code 3P0G into the search box. After clicking on Go, you should be presented with on overview of tools the ProteinPlus server provides.

We do not go into great detail on all the tools but only mention PoseView. With this tool, you can prepare an automatic sketch of the small molecule-protein interactions.

![Protein Plus Server](../../images/ProteinPlusPoseView.png "Overview of 3P0G")
![Zoom on ligand of 3P0G](../../images/3P0G_A_PoseView_Input.png "Zoom on ligand co-crystallized with 3P0G")


> #### {% icon question %} Questions
>
> 1. Between which amino acids and the ligand do you see hydrogen bonds using YASARA? 
> 2. According to PoseView, between which amino acids and the ligand do you see hydrogen bonds?
> 3. What other interactions are presented in the sketch?
> 4. Inspect the visualisation in Yasara: Do you see the interactions in Yasara as well?
>
> <details markdown="1">
> <summary>{% icon solution %} Solution
> </summary>
> 
> 1. In YASARA, you observe hydrogen bonds between Asp113A as well as the carbonyl function of Asn312A and the charged amine function.
>    
> 2. PoseView indicates hydrogen bonds between Asp113A as well as the carbonyl function of Asn312A and the charged amine function. Furthermore, hydrogen bonds are indicated between the phenolic OH and Ser207A and Ser203A as well as the amine function and Ser203A.
> 
> 3. Furthermore, hydrophobic interactions are indicated for the methylbenzyl moiety and pi-pi interactions of Phe290A and the phenolic moiety.
>
> 4. With YASARA Structure license, those hydrophobic interactions can also be visualised. 
> </details>
{: .question}


### Exploring the structure of a nanobody-stabilized active state of the β2 adrenoceptor - the nanobody 

In order to estimate the binding energy between the nanobody and the β2 adrenoceptor, we can use the FoldX tool AnalyseComplex. It is recommended to calculate these binding energies on energy-minimized structures. To illustrate the effect of the energy minimization, we compare the interaction energy of the current crystal structure and its minimized structure.

### Use the tool FoldX tool AnalyseComplex 

> #### {% icon hands_on %} 
>
> 1. Given that energy-minimization takes a while for this rather large complex,
>     please download the Yasara scene [here](http://data.bits.vib.be/pub/trainingen/PSA/3P0G_1.sce)  
>    
>    Calculate the interaction energies between the chain A and B of the object 3P0G 
>    and the RepairObj1, respectively. 
>
>    ```
>    Analyze - FoldX - Interaction energy of molecules
>    ```
{: .hands_on}

> #### {% icon question %} Questions
>
> 1. What is the dG in the two cases? 
> 2. Any idea why the difference is rather hugh?
>
> <details markdown="1">
> <summary>{% icon solution %} Solution
> </summary>
> 
> 1. first case (X-ray structure): Interaction energy between molecule(s) A and B in object 1 = -9.86 (kcal/mol)
>    second case: 
>    Interaction energy between molecule(s) A and B in object 2 = -20.19 (kcal/mol)
> 2. Through the energy minimisation of the Repair Object function, the interactions of the amino acids are optimised.  
> </details>


This command also creates a list of residues forming the interface of the two proteins. Hit space to see the list of residues in the interface.

Tip: This list can also be useful if you want to make visualisations of the interaction site.

```
Plugin>interface residues between A and B
Plugin>TA66 TA68 IA72 IA127 RA131 AA134 IA135 TA136 SA137 PA138 FA139 KA140 QA142 YA219 VA222 EA225 AA226 LA266 KA267 EA268 AA271 LA272 TA274 LA275 IA278 IA325 YA326 RA328 SA329 PA330 SB27 IB28 FB29 SB30 IB31 TB33 AB50 IB51 eB52 SB56 TB57 NB58 YB100 AB102 VB103 LB104 YB105 EB106 YB107
```

### Comparing the active and the inactive conformation of the β2 adrenoceptor 

In case, there is still time, I would recommend to try to use some of your capabilities you learned today and create a superposition of the inactive and active conformation of the β2 adrenoceptor. We take one of the crystal structures which are available: 3SN6

```
File - Load - PDB file from Internet
```

You will be kind of overwhelmed once the structure is loaded into YASARA. In order to get a first quick overview, click on the 'Center' buttom in the menu of YASARA (5th buttom from the left). Then, it is time to look at the PDB entry of 3SN6 in the PDB database to have a first idea on what molecules are in the PDB file.

As you see on the website [3SN6](http://www.rcsb.org/pdb/explore/explore.do?structureId=3SN6i), the chain R consists of 2 molecules, the β2 adrenoceptor and lysozyme. 
In the corresponding article, it is stated that 'the unstructured amino terminus of the β2AR is replaced with T4 lysozyme (T4L)'.

Since this is an extra molecule in the crystal structure which disturbes our view, we will delete it.

After the manipulation, the overall picture should look roughly like this.

![Superposition](../../images/3SN6_withoutLysozyme.png "Overview of 3SN6 without lysozyme")

In the following step, we structurally align only the receptors. The rest of the structures will move along.
It is suggested to use the first chain A from 3P0G as target. In order to do a structural alignment, it is suggested to use the first chain A from 3P0G as target.

```
Analyze - Align - Pairwise, based on structure - Molecules with MUSTANG
```

Investigate the differences in TM helices and the binding of the nanobody compared to the subunit of the G protein.

Tip: Color the secondary structures to better identify the individual chains/units of G protein.

## On-Line References

- [Crystallography 101 (Bernhard Rupp)](http://www.ruppweb.org/Xray/101index.html)

- [Protein NMR, A Practical Guide (Vicky Higman)](http://www.protein-nmr.org.uk/)

- [Model validation course](http://xray.bmc.uu.se/gerard/embo2001/modval/index.html)

- [Assessing model quality](http://spdbv.vital-it.ch/TheMolecularLevel/ModQual/)

- [Lectures by Burkhard Rost on protein structure prediction](https://www.youtube.com/channel/UCU6j8BG4RbEtTgyIZJ6Vpow)

## Annex 

### Older methods before we had AlphaFold to predict protein structures by fold recognition

1. Search SCOP/CATH for protein with same fold and known 3D structure
2. Align each amino acid of query sequence to a position in the template structure
3. Evaluate how well the sequence fits the fold and select best-fit fold
4. Build structural model of query based on alignment with selected fold

- HHpred (http://toolkit.lmb.uni-muenchen.de/hhpred)
- Phyre (http://www.sbg.bio.ic.ac.uk/phyre2/html/page.cgi?id=index)
- DescFold (http://202.112.170.199/DescFold/)

Works because:

- Number of different folds in nature is fairly small (approximately 2000-3000)
- 90% of new submissions in PDB have similar folds to those already in PDB
- Not always accurate [^1](https://onlinelibrary.wiley.com/doi/full/10.1002/prot.25823)

[^2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3003448/#!po=87.5000)

### Guidelines to improve fold recognition results

- Run as many methods as you can
- Run each method on many sequences from your homologous protein family
- After all of these runs, build up a consensus picture of the likely fold
- Compare function of your protein to function of the proteins with the likely fold
- Compare secondary structure of your protein to that of the likely fold

### Homology Modelling

- When a structure is available for a protein with a  similar sequences, it is possible to predict the  structure of a new sequence with varying degrees of  confidence.
- Use PSI-BLAST/DELTA-BLAST to detect sequences with similar structures.
- All homology modelling procedures start from an  alignment of the template and target sequences.  The quality of this alignment will have a major  impact on the resulting model.
- Available applications include stand-alone programs  (Modeller, FoldX, …) and web-based services (such as SwissModel).

### Exercise: Make a Homology Model using Swiss Model

[Link to exercise 5](https://material.bits.vib.be/topics/protein-structure-analysis/tutorials/homology-modeling/tutorial.html)

### Force Fields

- Energy terms representing physical interactions

  - Covalent bond lengths and angles
  - Dihedral angles and van der Waals forces (steric effects)
  - Electrostatic interactions and hydrogen bonds
  - …

- Energy can be minimized, and forces can literally be derived from the potential function.

- Consistency and careful consideration of the properties to be simulated are essential.

### Force Field Terms

Each energy term has a functional form, which includes one or more parameters:

- Covalent bond energy term

  $ E_{bond} = k(r-r_{0})^{2} $ with $r_0$ being the equilibrium length

- Van der Waals contact energy term

  $ E_{VdW} = C_6r^{-6} + C_{12}r^{-12} $ with $C_6$ attractive, $C_{12}$ repulsive term

The parameters are collectively optimized to  reproduce a chosen set of experimentally observed parameters.

A given force field should be used as a  consistent system, and can only be used to predict properties that are covered by the  training set.

### FoldX

Is designed for quantitative modelling of the contribution of structural interactions to the stability of proteins and protein complexes. It also supports  protein/DNA complexes.

The force field describes the different interactions in a protein structure or complex in analytical terms. It has  been calibrated using a set of experimentally determined  stabilities.

Applications include the optimisation of structures, the calculation of the stability of complexes, and predicting  the effect of amino-acid or base-pair mutations on these  properties.

![FoldX Force Field](img/foldx-formula.PNG)

[FoldX](https://dx.doi.org/10.1093%2Fnar%2Fgki387)

### The FoldX Plugin for YASARA

In order to make FoldX more accessible and integrate its functions into Yasara, dr. Joost van  Durme (SWITCH laboratory) made a Yasara plugin  module that can apply FoldX functions to structures  that are loaded as Yasara objects.

This greatly simplifies the use of FoldX, and allows  for a quick visual analysis of the resulting changes  in the structures.

More information can be found at [wiki](http://foldxyasara.switchlab.org/index.php/) [FoldX](http://foldxsuite.crg.eu/)

### Exercise 4a: Repair a PDB File

[Link to exercise 4](https://material.bits.vib.be/topics/protein-structure-analysis/tutorials/mutate-structure/tutorial.html)

- Load the 1CRN PDB file.
- Use the “Repair object” option in the `Analysis|FoldX` menu to activate the corresponding FoldX function.
- Select the object to repair.

This exports the object as a temporary PDB file,  starts FoldX with the appropriate options, and loads the repaired PDB file as a new object in Yasara.

- Compare the original and repaired objects.  
- Describe the changes that were introduced.

### Example tutorials

**Docking runs**

https://ringo.ams.stonybrook.edu/index.php/2023_DOCK_tutorial_2_with_PDBID_3WZE#ChimeraX_.28optional.29

**Assemblies**

https://static-bcrf.biochem.wisc.edu/tutorials/chimerax/Chimera_X_overview_document.html#12_biological_assembly

https://pubmed.ncbi.nlm.nih.gov/35001912/

**Foldseek**

https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac072/6749558
