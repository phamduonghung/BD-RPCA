# BD-RPCA

This MATLAB package is a collection of scripts allowing to generate figures in the paper [1]. This paper addresses the problem of high-resolution Doppler blood flow estimation from an ultrafast sequence of ultrasound images. Formulating the separation of clutter and blood components as an inverse problem has been shown in the literature to be a good alternative to spatio-temporal singular value decomposition (SVD)-based clutter filtering. In particular, a deconvolution step has recently been embedded in such a problem to mitigate the influence of the experimentally measured point spread function (PSF) of the imaging system. Deconvolution was shown in this context to improve the accuracy of the blood flow reconstruction. However, measuring the PSF requires non-trivial experimental setups. To overcome this limitation, we propose herein a blind deconvolution method able to estimate both the blood component and the PSF from Doppler data. Numerical experiments conducted on simulated and in vivo data demonstrate qualitatively and quantitatively the effectiveness of the proposed approach in comparison with the previous method based on experimentally measured PSF and two other state-of-the-art approaches.


Requirements: 
1. To print nice pdf figures, the "export_fig" package was used, whic required software support downloaded from the following link: https://www.ghostscript.com/download/gpcldnld.html




[1] D.-H. Pham, A. Basarab, I. Zemmoura, JP. Remenieras, and D. Kouame, "Joint Blind Deconvolution and Robust Principal Component Analysis for Blood Flow Estimation in Medical Ultrasound Imaging," Accepted for publication in IEEE Ultrasonics, Ferrolectrics, and Frequency Control, September 2020. Available: https://arxiv.org/abs/2007.05428.

COPYRIGHT

BD-RPCA is copyright reserved. For further information, please contact Duong Hung PHAM at duong-hung.pham@irit.fr.
