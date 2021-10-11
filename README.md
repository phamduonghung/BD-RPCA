# BD-RPCA

This MATLAB package is a collection of scripts allowing to generate figures (Fig. 1 and Fig. 2a-2e) in the paper [1]. This paper addresses the problem of high-resolution Doppler blood flow estimation from an ultrafast sequence of ultrasound images. Formulating the separation of clutter and blood components as an inverse problem has been shown in the literature to be a good alternative to spatio-temporal singular value decomposition (SVD)-based clutter filtering. In particular, a deconvolution step has recently been embedded in such a problem to mitigate the influence of the experimentally measured point spread function (PSF) of the imaging system. Deconvolution was shown in this context to improve the accuracy of the blood flow reconstruction. However, measuring the PSF requires non-trivial experimental setups. To overcome this limitation, we propose herein a blind deconvolution method able to estimate both the blood component and the PSF from Doppler data. Numerical experiments conducted on simulated and in vivo data demonstrate qualitatively and quantitatively the effectiveness of the proposed approach in comparison with the previous method based on experimentally measured PSF and two other state-of-the-art approaches.


## Instructions
1. Download the package in .zip file (click green Code above) and then unzip it. Note that the **name** of the unzipped folder should be **BD-RPCA**.  
2. Set **Current Folder** of MATLAB to this unzipped folder, i.e. **BD-RPCA**.  
3. Download all **simulation** data from the following link: 
https://cloud.irit.fr/index.php/s/846gUKURnYbehVl and then put them into the folder **Data**
4. Run each file **Fig?.m** corresponding to each figure (from Fig. 1 to Fig. 2a-2e) in [1]. 
5. To print nice pdf figures, the **export_fig** package was used, which required a supporting software downloaded from the following link (there is also a portable version of this software): https://www.ghostscript.com/download/gpcldnld.html. In the codes, change **FigFeatures.print= 1** if you want to print the .pdf figure using this package. 


[1] D.-H. Pham, A. Basarab, I. Zemmoura, JP. Remenieras, and D. Kouame, "Joint Blind Deconvolution and Robust Principal Component Analysis for Blood Flow Estimation in Medical Ultrasound Imaging," in *IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control*,vol. 68, no. 4, pp. 969-978, April 2021. Available: https://arxiv.org/pdf/2007.05428.pdf.


## Copyright
RobustPCA_Doppler.m comes from https://github.com/dlaptev/RobustPCA.

BD-RPCA is copyright reserved. If any issue related this package appears, please contact Duong Hung PHAM at duong-hung.pham@irit.fr.
