# SynQuant
SynQuant is a Fiji plugin that automatically quantify synapses from multi-channel fluorescence microscopy images. Both synapse and corresponding dendrite are detected. Synapses is detected on synapse channel, where they act as puncta surrounded by highly inhomogeneous interference signals. Dendrite is extracted from the reference dendrite channel.

SynQuant detect synapses through a totally unsupervised probability principled framework. In this framework, analysis is conducted on salient regions rather than pixels. All synapse candidates are scored by their own local contrast and compared fairly with each other. What’s more, false discover rate (FDR) control is utilized to determine synapse selection, which not only controls the false positive rate but also provides a statistical evidence of the detected synapse. The parameter used in this framework is only the value of FDR which is easy to tune. SynQuant extract dendrite by steerable filter [1]. Extracted dendrite then are segmented into roughly homogeneous pieces by branch points and end points. Based on the dendrite pieces and synapses, linear regression is used to find the effects of dendrite’s properties to the number of synapses on it.

The paper for this algorithm is in preparation.



[1].	Meijering E, Jacob M, Sarria J, Steiner P, Hirling H, Unser M. Design and valida-tion of a tool for neurite tracing and analysis in fluorescence microscopy im-ages. Cytometry A 2004;58:167–76.
