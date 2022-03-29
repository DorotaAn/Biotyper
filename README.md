# Raw MS data pipeline for distinguishing bacterial subspecies - Biotyper

## A tool for processing raw MS data of bacterial proteome 

## Motivation

To distinguish between bacterial sub types, Mass Spectrum (MS) data can be used. It allows identification of differences in protein expression between bacterial samples treated and untreated with antibiotic. This pipeline allows processing of raw MS data required for analysis and visualization of the results, which in turn can help distinguish between the bacterial sub types based on their proteome behaviour in response to antibiotic.

## Project description

The pipeline allows processing of the raw Mass Spectrum data directly from .dat files produced by MALDI TOF MS. Obtained, preccessed data can be then used for further analysis of the results. Each step of the pipeline processes spectra and metadata simultaniously. 
Initially the focus was on filtering out spectra of low quality. That was followed by transformation of itensities using square root, wave smoothing and removal of the spectrum baseline. Since the spectra ID was lost at the baseline removal, retenetion step was included prior calibration of the itensity. Based on the obtained spectra the range can be trimmed.
The second step involves detection and alignment of peaks based on pre set parameters, which should be adjusted according to the data. At this point parameter for rare peak selection should be considered for their correct removal if that is considered to be useful. Additionally, outliers can be identified and removed.
Finally, the metadata and produced peak itensities can be merged and processed according to the needs. In this case, data was saved for later analysis as a whole, as well as merged by Strain, BioRep, Antibiotic and Time of Antibiotic treatment for analysis by BioRep with retention of other required information. 

The pipeline was based on ilustrative pipeline (Palarea-Albaladejo et al., 2018) and was builed using very popular and useful package MALDIquant (Gibb and Strimmer, 2012) as well as MALDIrppa (Palarea-Albaladejo et al., 2018). MALDIquant package contains useful packages for processing of MS data. MALDIrppa was created for the purpose of this project and was fine tuned to maximase efficiency of the tasks performed in this pipeline.

## Credit
This pipeline was created in collaboration between Moredun Research institute and Biomathematics and Statistics Scotland. 

## References

Gibb S, Strimmer K (2012). “MALDIquant: a versatile R package for the analysis of mass spectrometry data.” Bioinformatics, 28(17), 2270–2271. doi: 10.1093/bioinformatics/bts447, https://academic.oup.com/bioinformatics/article-pdf/28/17/2270/682111/bts447.pdf, https://academic.oup.com/bioinformatics/article-abstract/28/17/2270/246552.

Palarea-Albaladejo J., McLean K., Wright F. and Smith (2018). MALDIrppa: quality control and robust analysis for mass spectrometry data. Bioinformatics 34(3):522–523. <doi: http://dx.doi.org/10.1093/bioinformatics/btx628>
