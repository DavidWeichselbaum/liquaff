# Liquid-like Affinity Landscape (LAL) 

LAL predicts RNA/protein interaction propensities and interacting regions from physicochemical first principles.

It so by finding regions of minimal interaction energies in a RNA/protein "affinity landscape".
A landscape is a table of propensities between the monomers of the two polymers (see *'./example_data/escherichia102013.CDS.db'*) arranged by sequence. 
The table is populated using a base/residue affinity scale (eg. *'./example_data/set2p'*).
A dynamic programming algorithm is employed to generate energy-sums of all possible sub-tables in quadratic time, instead of O(m³n³).

## Example
Compile *'lal'*:
> make

Run standard analysis:
> ./lal example_data/escherichia102013.CDS.db example_data/set2p

## Options
Usage: *affinity [options] \<db\> \<scale\>*
  -e print energy matrices
  -i specify rna/prot pair by index
  -n specify rna/prot pair by name (uniprot ID) (overrides -i)
  -h display this and exit
  -s schedule file

## Formats
Input format (tab separated):     protein/RNA ID; protein/RNA sequence; sequence length; identifier [*'RNA'* or **'PROT'**] 

Output format (tab separated):		protein ID; protein length; RNA ID; RNA length; protein start index; protein region length; RNA start index; RNA region length; energy of region; boltzmann probability

Format schedule:	Integer slices of protein/RNA sequence indexes with exclusive upper bound divided by a dash. Schedules are separated by newline.
Examples for valid schedule formats: 'p:p-r:r' ':-r:r' 'p:p-'  '-r:' ':-:'  '-'

## Thanks

The sequence file *'./example_data/escherichia102013.CDS.db'* was provided by Mario Hlevnjak¹, the affinity scale file *'./example_data/set2p'* by Anton A. Polyansky².

## References
1) Hlevnjak,M. et al. (2012) Sequence signatures of direct complementarity be-
tween mRNAs and cognate proteins on multiple levels. *Nucleic Acids Research*, 40, 8874–8882
2) Polyansky A. and Zagrovic B. (2013) Evidence of direct complementary interactions between messenger RNAs and their cognate proteins. *Nucleic Acid Research*, 41, 8434–8443
