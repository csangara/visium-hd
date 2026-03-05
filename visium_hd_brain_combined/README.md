## Compare results of FFPE and fresh frozen VisiumHD brain sections

This folder creates Figure 4.7 and other supplementary figures. It also downloads the MERFISH sections that were used as ground truth (see 
from the [Allen Brain Cell Atlas - Data Access](https://alleninstitute.github.io/abc_atlas_access/intro.html).

* 0: contains the color palettes of the cell types and the regions from the Allen Brain Institute
* 1_1 and 1_2 contain the same code to download the MERFISH sections, but are just written as a Python script and notebook, respectively
* 1_3: Creates rds files of the MERFISH dataset, creates Supplementary Figure 4.9
* 1_4: Creates Supplementary Figure 4.12 and Figure 4.2
* 2_1: Creates Figure 4.7b
* 2_2: Creates Figure 4.7a (other components are added via Inkscape), Figure 4.7c, Supplementary Figure 4.11. Note that Supp Figure 4.11 was combined in Inkscape. The FPR values are also manually added in Inkscape (these values can be found in the `deconv_fpr_per_region.csv` file).

