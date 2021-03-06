This file explains how to use the provided code and gives some calculation details


=== Gene length data ===
The gene length data can be taken from
https://www.sdbonline.org/sites/fly/segment/hunchbk1.htm
or
https://www.ncbi.nlm.nih.gov/gene/41032
I used the first option to get Snail and Knirps lengths



=== MS2 loops length ===
Figure 4 of Desponds 2016 suggests that two MS2 loops contains
- 2 loops of 4 nt,
- 2 stems of 8 nt,
- 2 stems of 7 nt.
- 1 short linker of 20 nt,
- 1 long linker of 50 nt.
Total for 2 MS2 loops: 8 + 16 + 14 + 20 + 50 = 108 nt.
Total for 24 MS2 loops: 12 * 108 = 1296 nt.
It will take a polymerase around 49.8 s to pass this region at an elongation rate of 26 bp/s (Garcia et al., 2013).



== output_data_table_for_latex.m ==

This function outputs the mean normalized values of the slope and steady-state polymerase numbers to be used in the LaTeX article.
Launch and wait to finish




