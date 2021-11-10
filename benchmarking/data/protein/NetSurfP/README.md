Original NetSurfP 2.0 data was obtained from PISCES. Get data from PISCES at the highest identity threshold and partition from there.  

http://dunbrack.fccc.edu/pisces/download/  

Drop empty lines and fix headers with a | :  
```
awk 'NF' cullpdb_pc95.0_res0.0-2.5_noBrks_len40-10000_R0.3_Xray_d2021_10_29_chains30046.fasta | sed 's/ /|/' >netsurfp.fasta
```