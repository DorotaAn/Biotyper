# Raw MS data pipeline for distinguishing bacterial subspecies - Whole Cell MALDI

## A tool for processing raw MS data of bacterial proteome 

## Motivation

To distinguish between bacterial sub types, Mass Spectrum (MS) data can be used. It allows identification of differences in protein expression between bacterial samples treated and untreated with antibiotic. This pipeline allows processing of raw MS data required for analysis and visualization of the results, which in turn can help distinguish between the bacterial sub types based on their proteome behaviour in response to antibiotic.

## Project description

The pipeline allows processing of the raw Mass Spectrum data directly from .dat files produced by MALDI TOF MS. Processed data can then be used for further analysis of the results. Each step of the pipeline processes spectra and metadata simultaneously. 
Initially the focus is on filtering out spectra of low quality. This is followed by transformation of intensities using square root, wave smoothing and removal of the spectrum baseline. As the spectra ID was lost at the baseline removal, a retention step was included prior to normalisation of intensity. Based on the obtained spectra the range can be trimmed.
The second step involves detection and alignment of peaks based on pre-set parameters, which should be adjusted according to the dataset being analysed. At this point parameters for selecting and removing rarely detected peaks can be optimised. Additionally, outliers can be identified and removed.
Finally, the metadata and processed peak intensities can be merged and analysed according to the user needs. 
The pipeline was based on illustrative pipeline (Palarea-Albaladejo et al., 2018) and was built using popular package MALDIquant (Gibb and Strimmer, 2012) and MALDIrppa (Palarea-Albaladejo et al., 2018). MALDIquant  contains useful packages for processing of MS data. MALDIrppa was created for Whole Cell MALDI analysis and was fine tuned to maximise efficiency of the tasks performed in this pipeline.

## Credit
This pipeline was created in collaboration between Moredun Research institute and Biomathematics and Statistics Scotland. 

## References

Gibb S, Strimmer K (2012). “MALDIquant: a versatile R package for the analysis of mass spectrometry data.” Bioinformatics, 28(17), 2270–2271. doi: 10.1093/bioinformatics/bts447, https://academic.oup.com/bioinformatics/article-pdf/28/17/2270/682111/bts447.pdf, https://academic.oup.com/bioinformatics/article-abstract/28/17/2270/246552.

Palarea-Albaladejo J., McLean K., Wright F. and Smith (2018). MALDIrppa: quality control and robust analysis for mass spectrometry data. Bioinformatics 34(3):522–523. <doi: http://dx.doi.org/10.1093/bioinformatics/btx628>

