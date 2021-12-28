# kaplan-meier-sarcoma

This script is a r shiny GUI that generates a Kaplan Meier plot based on CbioPortal sarcoma data set. user chooses two genes and the gene expression level, the software divides the data into four groups High/Low expression/Gene x/Gene Y and generates a Kaplan Meier plot. 

The "table_generator_sarcoma.R" receives as input a list of genes and one driver gene, thus outputs a list with the Kaplan Meier p-value per gene for the same four groups categorization explained above.
