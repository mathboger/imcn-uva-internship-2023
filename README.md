# Disclaimers

This repository contains the result of work done from February through July of 2023 by Matheus Boger in the IMCN (Integrative Model-based Cognitive Neuroscience) research unit at the University of Amsterdam led by Birte Forstmann under the supervision of Niek Stevenson and Steven Miletic as a research intern as part of the coursework of the Brain and Cognitive Sciences Research Master at the University of Amsterdam.

Unless stated otherwise, all code regarding the volatile Kalman filter (VKF) was based on the mathematical definitions proposed by Piray and Daw, 2020 (full citation in text). 

The original data from Miletic et al., 2021 used in the analysis is available at: https://osf.io/txgb3 .

# Abstract

Copied directly from the final report:

Cognitive modeling is an important tool to investigate behavior and the brain. To this end, two large classes of models have been developed in the past: reinforcement learning (RL) models to portray the changes in behavior as a response to its consequences, and evidence accumulator models (EAMs) to portray the underlying processes of decision-making, including how long it takes to reach a decision. However, the combination of the two methodologies into one is recent, and few is known on how more complex RL learning rules that incorporate environmental volatility in their definitions could be used along EAMs. In this study, the volatile Kalman filter (VKF), as proposed by Piray and Daw, 2020, is tested as such a possible learning rule on a reversal learning task conducted by Miletić et al., 2021. A parameter recovery analysis was also conducted in order to ascertain the robustness of the VKF itself. Results indicated that the VKF has poor recovery of its parameters, especially the ones regarding the modeling of volatility, and that the simplest integration possible of it to an LBA is not enough to make it a good descriptor of participant accuracy and reaction time (RT) data. To solve this issue, it is recommended that future work focus on how the VKF can be improved to be more discrete in nature, or that a more complex linkage to an EAM with an urgency signal be attempted. Nevertheless, the increase of model complexity brought by the VKF was demonstrated to not worsen the BIC of the overall model, which gives hope to the overarching class of multi-layered learning rules.

# Folder structure

In the root folder you can find this README and the final report itself. The rest is as follows:

- Data
    - Fitted models (RData files with the estimated parameters for each model and participant combination)
    - Simulated data (RData files with data generated for the profile plot, parameter recovery and simulation for the descriptive adequacy test for each model and participant combination)
    - Original data (empty; place here the original from Miletic et al., 2021 for the R scripts to work)
- Figures (all figures that went into the final report as separate files)
- R Scripts (R scripts developed as the final delivery of the internship)