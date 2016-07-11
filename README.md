# SynQuant
SynQuant is a Fiji plugin that automatically quantify synapses from multi-channel fluorescence microscopy images. Both synapse and corresponding dendrite are detected. Synapses are detected on synapse channel, where they act as puncta surrounded by highly inhomogeneous interference signals. Dendrite is extracted from the reference dendrite channel.

SynQuant detect synapses through a totally unsupervised probability principled framework. In this framework, analysis is conducted on salient regions rather than pixels. All synapse candidates are scored by their own local contrast and compared fairly with each other. What’s more, false discover rate (FDR) control is utilized to determine synapse selection, which not only controls the false positive rate but also provides a statistical evidence of the detected synapse. The parameter used in this framework is only the value of FDR which is easy to tune. The framework of synapse detection algorithm now is based on the idea of component tree [1][2]. SynQuant extract dendrite by steerable filter [3]. Extracted dendrite then are segmented into roughly homogeneous pieces by branch points and end points. Based on the dendrite pieces and synapses, linear regression is used to find the effects of dendrite’s properties to the number of synapses on it.

# Update
Add the component tree structure for synapse detection suggested by Dr.Petter Ranefall. The version 1.1 now is 3 times faster than before.

The paper for this algorithm is in preparation.
# How to cite
The journal paper for this algorithm is in preparation. You can cite the conference paper first:

Yizhi Wang, Guilai Shi, Yinxue Wang, Lin Tian, Guoqiang Yu*, “PPSD: Probability Principled Synapse Detection,” BioImage Informatics Conference 2015, October 2015.

# Tutorial
Please see the user manual.

[1].	Ranefall, Petter, Sajith Kecheril Sadanandan, and Carolina Wählby. "Fast Adaptive Local Thresholding Based on Ellipse fit." International Symposium on Biomedical Imaging (ISBI'16), Prague, Czech Republic, April 13-16, 2016. 2016.

[2].  Najman, Laurent, and Michel Couprie. "Building the component tree in quasi-linear time." IEEE Transactions on Image Processing 15.11 (2006): 3531-3539.

[3].  Meijering E, Jacob M, Sarria J, Steiner P, Hirling H, Unser M. "Design and valida-tion of a tool for neurite tracing and analysis in fluorescence microscopy images." Cytometry A 2004;58:167–76.
