- For **Countmatrix.txt :** A data frame of non-negative integers(the read counts).The first column must be `genename` (generally could be `ENSEMBL` or `SYMBOL` id. And please note that: the column `genename` must be unique), and the other columns are your `sample names`.
- For **ExperimentalDesign.txt :** Data frame, Experimental design at least with 'label', 'ID', 'condition' and 'replicate' columns. Of course, other information columns that you are interested in could be included.
  + **label :** Label names
  + **ID :** The same as label
  + **condition :** Experimental conditions
  + **replicate :** Replicate number
  + **Other information columns that you are interested in**