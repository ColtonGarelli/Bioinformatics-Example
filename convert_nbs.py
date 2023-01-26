import sos_rmarkdown
from sos_rmarkdown.converter import RmarkdownToNotebookConverter

rc = RmarkdownToNotebookConverter()

rc.convert("./GSE119207_analysis.Rmd","./GSE119207_analysis.ipynb")

