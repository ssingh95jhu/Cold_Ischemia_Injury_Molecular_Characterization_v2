This folder contains the code demonstrating how to integrate all the 17 spatial transcriptomics (ST) datasets mentioned in the current study.
The ST datasets belong to the following experimental groups namely:
1. Cold Ischemia Injury (CIS)
2. Warm Ischemia-Reperfusion Injury (AKI) (female mice)
3. Native kidneys (CTRL)
4. 24-hours Warm Ischemia-Reperfusion Injury (AKI24) (male mice)


Note: 
1.AKI24 is same as IRL group (as mentioned in the code)

2.For AKI Dataset, refer to the following article:
Dixon+, J Am Soc Nephrol., 2022 Feb;33(2):279-289.  doi: 10.1681/ASN.2021081150 (PMID: 34853151)

3. CTRL and IRL datasets were obtained upon request from authors of the following article:
Gharaie+, Sci Rep., 2023 Nov 28;13(1):20888.  doi: 10.1038/s41598-023-48213-2 (PMID: 38017015)

===============================================================================================

The code makes use of Harmony package to remove batch effects among the ST datasets and perform subsequent clustering which facilitates idenitification
of shared compartments compartments across all the 17 ST datasets namely: inner medulla, outer medulla and cortex.

Note: 
1. inner medulla is referred to as medulla (as mentioned in the code)
2. outer medulla is referred to as interface (as mentioned in the code)
