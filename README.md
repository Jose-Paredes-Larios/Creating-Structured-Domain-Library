# Creating a Structured Domain Library

## Introduction
Biological condensates are membrane-less compartments that concentrate biomolecules to regulate processes, such as gene expression, RNA processing, signal transduction, and stress responses. Because of their central role in cellular organization, dysregulation of condensate formation and dynamics has been implicated in diseases, including cancer and neurodegenerative disorders. Thus, understanding how biomolecules intrinsically encode their propensity to condense is important at both fundamental and translational levels. 

Recently, CondenSeq has been developed and used to screen the effects of thousands of protein sequences on nuclear condensates [(Kappel et al., 2025)](https://www.nature.com/articles/s41592-025-02726-y). The study primarily found that aromatic residues and higher net charge per residue promote condensate formation across sequence contexts. However, only intrinsically disordered regions were screen. Certain structured domains are known to affect a protein’s propensity to condense. Here, we are generating a library of structured domains to perform CondenSeq, enabling us to assess how these domains contribute to condensate formation.

## Step 0: Downloading the Human Proteome
To maintain physiological relevance, our library will only contain structured domains found in the human proteome. Thus, we begin by downloading the human proteome from [UniProt]( https://www.uniprot.org/proteomes/UP000005640). When downloading:
-	Only download reviewed proteins
-	Download in TSV format
-	Download the compressed file
-	Under UniProt Data, only select: Entry name, Length, Sequence, Gene Name, and Domain [FT]

Note that this pipeline was built on the assumption that all the data I have mentioned has been included in the download. To decompress your file, go to your terminal and run `gunzip -c [file path of the tsv downloaded from Uniprot] > [file path of where you want the decompressed tsv to be]`

## Step 1: Extracting Domains from the Human Proteome
The domain info stored in the human proteome download is quite messy. Thus, run `1_domainInfoExtractor` on your human proteome download by changing the input and output pathnames in the file. After running the code, a new TSV file will be generated listing every domain sequence in the human proteome, including its name, length, and start and end positions within the parent protein. The parent protein’s UniProt accession number, gene name, and total length is also a part of each entry. 

Because of the technical limitations of CondenSeq, the code will also filter out any domain sequence that is longer than 66 amino acids long. However, you can comment this part out if desired.

## Step 2: Filtering out Disordered Domains
Uniprot provides us with domains found throughout the proteome; however, their annotated domains are not guaranteed to be structured. We can use a tool to predict how disordered a domain sequence is; if it is not too disordered, we can consider it structured. I decided to use [AIUPred](https://aiupred.elte.hu) to make disordered predictions for our domain sequences because of its high accuracy and how easy it is to import as a python library.

Before running `2_disorderedPredictions` on the output of step 1, go to the [AIUPred GitHub](https://github.com/doszilab/AIUPred) and download their repository into a folder. Place the folder into the same environment where you are running the pipeline. Specify pathnames (input, output, and AIUPred directory), and run `2_disorderedPredictions`. Like the [AIUPred paper](https://academic.oup.com/nar/article/52/W1/W176/7673484), any residue with a predicted value above 0.5 is considered disordered. Lastly, the file can either keep structured domains or disordered domains, depending on the line you filter out towards the end. 

## Step 3: Quantifying Domain Interactions
For the next steps of the pipeline, we need structures of each domain sequence and its parent protein. However,most of these sequences have not been experimentally solved. Thus, we will utlilize [AlphaFoldDB](https://alphafold.ebi.ac.uk) and work with predicted structures.

Because we want structured domains that can fold independent of its parent protein, we next quantify domain interactions with `3_alphaFoldDomainInteractions`. The idea is the less a domain interacts with its parent, the more likely it is to fold on its own. To help us quantify domain interactions, we calculate:
- **Anchoring Index:** A measure of how constrained a domain is by its parent. Helps us determine if domain B has to be around residues X-Y of the parent on average. Note that this is not a measure of interaction strength.
- **Fraction Buried:** Calculates the fraction of a domain that is buried within its parent protein.
- **Contanct Density:** Quantifies the number of interactions a domain has with its parent protein by calculating the fraction of total possible contacts found between the domain and protein. 

We take a weighted average of these three metrics to obtain a single score that reflects how much a domain is interacting with its parent. Note that before any of these calculations we take the pLDDT value (AlphaFold's confidence in a model) of each domain sequence and filter out low confidence sequences (<80). Lastly, to run `3_alphaFoldDomainInteractions` specify the pathnames of the input (the output of the last step), output, and where you want to store the AlphaFoldDB files that will be downloaded. 

## Step 4: Determing Physical Properties of Domains
For structured domains, we want physical properties relating to surface residues, as they are the ones that will interact with other biomolecules to influence condensation. Thus, `4_physicalPropertyDomainStruct` will calculate the following metrics for each domain sequence:
- Rg (Å)(Compactness)
- Fraction of Surface Residues
- Fraction of Aromatic Surface Residues
- Fraction of Positive Surface Residues
- Fraction of Negative Surface Residues

To run `4_physicalPropertyDomainStruct` specify the pathnames of the input (the output of the last step), output, and where you downloaded the AlphaFoldDB files to in the last step.

## Step 5: Obtaining and Visualizing Final Candidate Sequences
For each domain sequence, we now have a handful of metrics to describe it. We use them in `5_obtainFinalCandidateSequences` to identify sequences to experimentally test based on the following: 
- We want our structured domains to fold on their own and not have a drastic difference in behavior when isolated from its parent protein. Thus, we are interested in representative sequences with low interaction indexes
- Protein condensation can be driven through aromatic interactions and electrostatic interactions. Thus, we want sequences with a high fraction of aromatic surface residues and a high net charge of surface residues
- We don’t want sequences that are too compact (little surface residues) or too loose (IDR like behavior)
- We do not have any gold standards to set thresholds yet. Thus, thresholds are set relative to the dataset:
  - Bottom 30% for interaction threshold
  - Top 25% for aromatic threshold
  - Top 25% for charge threshold
  - Middle 50% for Rg threshold

To run `5_obtainFinalCandidateSequences` specify the pathnames of the input (the output of the last step) and output. Note that this is an R file, not a python file. 

## Step 6: Imaging Candidate Sequences with Pymol
To confirm that our final candidate sequences are not nonsensical, we can image each one in their parent protein with annotations and combine these images to create a tiled image with `6_pymolImages`. 

To run `6_pymolImages` specify the pathnames of the input (the output of the last step), the directory for individual images, and the final tiled image output.




