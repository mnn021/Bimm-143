Untitled
================

## The PDB Database for biomolecular structure data

### 1.1

> Q1: Download a CSV file from the PDB site (accessible from “Analyze”
> -\> “PDB Statistics” \> “by Experimental Method and Molecular Type”.
> Move this CSV file into your RStudio project and determine the
> percentage of structures solved by X-Ray and Electron Microscopy.

> Also can you determine what proportion of structures are protein? Aim
> to have a rendered GitHub document with working code that yields your
> answers.

download CSV file from PDB website(“Analyze” -\> “PDB Statistics” \> “by
Experimental Method and Molecular Type”)

``` r
# Read CSV
data <- read.csv("Data Export Summary.csv")
head(data)
```

    ##   Experimental.Method Proteins Nucleic.Acids Protein.NA.Complex Other
    ## 1               X-Ray   131278          2059               6759     8
    ## 2                 NMR    11235          1303                261     8
    ## 3 Electron Microscopy     2899            32                999     0
    ## 4               Other      280             4                  6    13
    ## 5        Multi Method      144             5                  2     1
    ##    Total
    ## 1 140104
    ## 2  12807
    ## 3   3930
    ## 4    303
    ## 5    152

``` r
# Total number of Entries 
sum(data$Total)
```

    ## [1] 157296

``` r
# Proportions of entries from each method
(data$Total/sum(data$Total)) * 100
```

    ## [1] 89.0702879  8.1419744  2.4984742  0.1926305  0.0966331

Proportion that are protein

``` r
sum(data$Proteins)
```

    ## [1] 145836

``` r
round(sum(data$Proteins)/sum(data$Total) *100, 2)
```

    ## [1] 92.71

One method is about 90% (X-Ray microscopy)

Q2: Type HIV in the PDB website search box on the home page and
determine how many HIV-1 protease structures are in the current PDB?

### 1.2: HIV-Pr Structure

Now download the “PDB File” for the HIV-1 protease structure with the
PDB identifier 1HSG. On the website you can “Display” the contents of
this “PDB format” file. Alternatively, you can examine the contents of
your downloaded file in a suitable text editor.

> NOTE: You can type 1HSG in the PDB search box to jump to its entry and
> then click “Download Files” to the right of the top display. Selecting
> “Display Files” will allow you to view the PDB file directly in your
> browser window.

Will use VMD… DOn’t HAVE find a way to get it …. File \>\> new molecule,
file name\>\> load

## HIV-Pr Structure analysis

Here we will read the 1HSG PDB structure and select the protein
component and write out a new **protein-only** PDB format file. We then
do **ligand-only** PDB file

``` r
library(bio3d)
```

``` r
pdb <- read.pdb("1hsg.pdb")

# Select chain A
a.inds <- atom.select(pdb, chain="A")
# Select C-alphas of chain A
ca.inds <- atom.select(pdb, "calpha", chain="A")
# We can combine multiple selection criteria to return their intersection
cab.inds <- atom.select(pdb, elety=c("CA","CB"), chain="A", resno=10:20)
```

### Functions

`library(bio3d)` `read.pdb` `write.pdb` `atom.select()` `trim.pdb()`

``` r
?atom.select

ligand <- atom.select(pdb, "ligand", value=TRUE)
ligand
```

    ## 
    ##  Call:  trim.pdb(pdb = pdb, sele)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 45,  XYZs#: 135  Chains#: 1  (values: B)
    ## 
    ##      Protein Atoms#: 0  (residues/Calpha atoms#: 0)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 45  (residues: 1)
    ##      Non-protein/nucleic resid values: [ MK1 (1) ]
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

``` r
write.pdb(ligand, file = "1hsg_ligand.pdb")



ptn <-atom.select(pdb, "protein", value=TRUE)
ptn
```

    ## 
    ##  Call:  trim.pdb(pdb = pdb, sele)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1514,  XYZs#: 4542  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 0  (residues: 0)
    ##      Non-protein/nucleic resid values: [ none ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

``` r
write.pdb(ptn, file = "1hsg_ptn.pdb")
```

with “ligand”, we see only 45 atoms, with a non-protein nuclic residue
of 1 (MK1)
