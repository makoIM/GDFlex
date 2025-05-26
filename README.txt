GDFlex README file

==================================================
The datasets (./data) used was downloaded from the sources below and processed to make it more accessible for the program.

Dau, A. H., et al. 2018. The UCR Time Series Archive. Retrieved
from https://www.cs.ucr.edu/~eamonn/time_series_data_2018/
==================================================

- To execute GDFlex, specify the target dataset ID from ML/UCR and the mode used in the Ablation study (baseline, znormBias, intra, locMis, inter, interNoise) described in the paper, then run PUB_GDFlex_exec from the command line. Among these modes, 'interNoise' represents the final version of GDFlex.

- To specify the target data, use either a data ID or "all". 

- When a specific data ID is specified, the following outputs are generated: the spreadsheet (result_all.csv), an intermediate result file (InfoDemo_(dataId).mat), and plotted five figures explaining the intermediate results.

- Note that selecting "all" will execute all 250 datasets, which may take several hours. With "interNoise", execution takes approximately 7 hours on a laptop with an Intel(R) Core(TM) i5-1335U CPU (1.30GHz) and 16GB of RAM. The baseline mode does not restrict the target data, so it takes more than a day. The other modes take less time than "interNoise".

% Change the directory to 'src' in the directory obtained by unzipping 'Code_forSupportingPage.zip'.
>> cd ./src 

[Example 1: Execution with a Specified Data ID]
>> Pub_GDFlex_execute(66, "interNoise")

[Example 2: Execution with All Data]
>> Pub_GDFlex_execute("all","interNoise")

==================================================
For details, please refer to HowToRunGDFlex.pdf.
It explains the plotted figures and the details of the output spreadsheet.
==================================================
