# Telomere Search

Python script to identify telomeric repeats at the start and end of contigs from a genome assembly

Jamie McGowan, 2021

**Requirements**

- Python >= 3
- BioPython

**Usage**

```
python TelomereSearch.py -i input.fasta -l length -t threshold
```

- `-i --input` fasta file containing genome contigs/scaffolds (**required**)
- `-l --length` scan the first and last N bp of each sequence [default = 200 bp]
- `-t --threshold` a telomere is reported if >= X% of bp scanned is composed of telomeric repeats [default = 0.4]. Cannot be used with `-m`
- `-m --min_copies` a telomere is reported if >= N copies of the telomeric repeat sequence is identified. Cannot be used with `-t`
- `-f --forward` regular expression to search start of contig
- `-r --reverse` regular expression to search end of contig 


For the `-l --length` parameter, something in the range between 100 bp and 600 bp seem like a sensible choice but will depend on the assembly quality and species. Similarily, `-t --threshold` or `-m --min_copies` should be adjusted depending on the assembly/species.

**Example Usage**

```
python TelomereSearch.py -i genome_assembly.fasta -l 100 -t 0.4
```


In this example, the script scans the first and last 100 bp of each contig for telomeric repeats. If the total length of telomeric repeats identified is >= 40 bp (i.e., 40% of the length scanned) a telomere is reported. The output is a `TSV` file showing the telomeric repeats identified for each contig. The parameters should be modified depending on the assembly or species. By examining the output, you can tell if the threshold or length parameters are too strict.


```
python TelomereSearch.py -i genome_assembly.fasta -l 100 -m 10
```

In this example, a telomere is reported if at least 10 copies of the repeat sequence is identified within the first/last 100 bp of the contig.


The script should be able to detect the most common telomeric sequences (`TTAGGG`/`CCCTAA`). The following regular expressions are used for searching `T{1,3}A{1,2}G{2,4}` and `C{2,4}T{1,2}A{1,3}`, which allows for some variation.

The regular expression used for searching for telomers can be modified depending on your organism of interest using the parameters `-f` and `-r`. For example, for ciliates, you might use `C{3,4}A{2,4}`
and `T{2,4}G{3,4}`


See [A Broad Phylogenetic Survey Unveils the Diversity and Evolution of Telomeres in Eukaryotes](https://academic.oup.com/gbe/article/5/3/468/582374) for a list of known telomeres.

See also [http://telomerase.asu.edu/sequences_telomere.html](http://telomerase.asu.edu/sequences_telomere.html) for a list of known telomeres.