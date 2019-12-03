lec12
================
Michael Nguyen
11/7/2019

## R Markdown

# Section 1: In silico docking of drugs to HIV-1 protease

## 1.1 Obtaining and inspecting our input structure.

### We want to produce a protein-only PDB file and a drug only PDB file

``` r
library(bio3d)

# download the PDB file 
get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

    ## [1] "./1hsg.pdb"

## 1.2 Prepare initial protein and ligand input files

Use the `trim.pdb()` function to make protein-only and ligand only
objects called prot and lig that you can then write out to new PDB
format files in your RStudio project directory.

``` r
pdb <- read.pdb("1hsg.pdb")
protein <- atom.select(pdb, "protein", value = TRUE)
write.pdb(protein, file ="1hsg_protein.pdb")
protein
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
ligand <- atom.select(pdb, "ligand", value = TRUE)
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
```

Q1: What is the name of the two non protein resid values in this
structure? What does resid correspond to and how would you get a listing
of all reside values in this structure? –\> HOH (water) and (wait no
HoH…) MK. Corresponds to the ID of the residue// residuals is a
generic function which extracts model residuals from objects returned by
modeling functions.

# 1.3 Using AutoDockTools to setup protein docking input

Docking algorithms require each atom to have a charge and an atom type
that describes its properties. However, typical PDB structures don’t
contain this information. We therefore have to ‘prep’ the protein and
ligand files to include these values along with their atomic
coordinates. All this will be done in a tool called AutoDock Tools (ADT
for short). Launch AutoDock Tools (ADT) from the Desktop icon on Windows
or using the Finder on Mac (from Applications \> MGLTools \>
AutoDocTools ).

Bonds and atoms are shown in white. For better visualization, color the
structure by atom type - Color \> By Atom Type. Click All Geometries and
then OK.

Q2: Can you locate the binding site visually? Note that crystal
structures normally lack hydrogen atoms, why? –\> No not really. Because
they are so small that in x-ray crystallopgraphy, H won’t be detected.

As we have already noted (remind me if we haven’t) crystal structures
normally lack hydrogen atoms. However, hydrogen atoms are required for
appropriate treatment of electrostatics during docking so we need to add
hydrogen atoms to the structure using Edit \> Hydrogen \> Add. Click OK.
You should see a lot of white dashes where the hydrogens were added.

Now we need to get ADT to assign charges and atom type to each atom in
the protein. We do this with Grid \> Macromolecule \> Choose…. Choose
1hsg\_protein in the popup window and click Select Molecule. ADT will
assign charges and prompt you to save the macromolecule. Save this
molecule to your class12 folder with the default file extension of
“pdbqt” (i.e. a new file called 1hsg\_protein.pdbqt). Open this in a
text editor and look at the last two columns - these should be the
charge and atom type for each atom. Note how they are different from the
1hsg\_protein.pdb file we input to ADT.

Q3: Look at the charges. Does it make sense (e.g. based on your
knowledge of the physiochemical properties of amino acids)? –\> sure

Section 2: Docking ligands into HIV-1 protease For this section, we will
use the program called Autodock Vina \[4\]. Autodock Vina is a fast
docking program that requires minimal user intervention and is often
employed for high- throughput virtual screening. We will run it from the
command line. Note: If you are working on your own laptop (i.e. NOT the
classroom machines) you will need to first install the AutoDoc vina
software for protein ligand docking from <http://vina.scripps.edu/>
download.html chose the appropriate version for your OS. Once downloaded
the Windows versions can be installed by double clicking the resulting
autodock\_vina\_1\_1\_2\_win32.msi file and following the install
prompts. For Mac simply double click the downloaded
autodock\_vina\_1\_1\_2\_mac.tar file in your Downloads directory to
yield a new directory called autodock\_vina\_1\_1\_2\_mac . Within this
new directory there is a bin directory containing the vina program.

## 2.2 Docking indinavir into HIV-1 protease

Now we are ready to run autodoc vina to perform the docking. We will
keep a log of all program output in a file `log.txt`

> ~/Downloads/autodock\_vina\_1\_1\_2\_mac/bin/vina –config config.txt
> –log log.txt

Once the run is complete, you should have two new files all.pdbqt, which
contains all the docked modes, and log.txt, which contains a table of
calculated affinities based on AutoDock Vina’s scoring function \[4\].
The best docked mode, according to AutoDock Vina, is the first entry in
all.pdbqt.

## 2.3 Inspecting your docking results

In order to visualize the docks and compare to the crystal conformation
of the ligand we will process the all.pdbqt to a PDB format file that
can be loaded into VMD. To do this we will use R and the Bio3D package.
Back in RStudio running in your class12 directory and with the bio3d
package loaded issue the following commands:

### Process our docking results

``` r
library(bio3d)
res <- read.pdb("all.pdbqt", multi=TRUE)
write.pdb(res, "results.pdb")
```
