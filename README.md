# GP_RTAP
Path-based gradient projection algorithm for solving the range-constrained traffic assignment problem (RTAP), which considers the heterogeneity of range anxiety (i.e., Discrete and Continuous) in the electrified transportation network. 
 
Several components are included in this package:

(1) Path: .../Source/DllProj/evNet
--ev.h: Definition of classes and data structures.
--EvTap.cpp: Definition of functions. For example: 1) "LabelSettingBSPTree" is equivalent to Algorithm 5 in Appendix A; 2) "SolveContinuousDCTAP_GP" is used to solve the continuous RTAP.
(2) Path: .../Source/ExeDriver/transitDriver
--evDriver.cpp: Main function   

(3) Path: .../Source/Networks/evNet
--Two-OD pair small network (Section 5.1 in the paper)
2odNet_net.dat: Network data
2odNet_trp.dat: Demand data

--Winnipeg network (Section 5.2 in the paper)
winnipeg_net.dat: Network data
winnipeg_trp.dat: Demand data

