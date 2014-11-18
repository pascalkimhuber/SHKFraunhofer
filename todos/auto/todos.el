(TeX-add-style-hook "todos"
 (lambda ()
    (LaTeX-add-labels
     "sec-1"
     "sec-2")
    (TeX-add-symbols
     '("alert" 1))
    (TeX-run-style-hooks
     "hyperref"
     "amssymb"
     "latexsym"
     "wasysym"
     "marvosym"
     "textcomp"
     "soul"
     "wrapfig"
     "float"
     "longtable"
     "graphicx"
     "fixltx2e"
     ""
     "fontenc"
     "T1"
     "inputenc"
     "utf8"
     "latex2e"
     "art10"
     "article"
     "a4paper")))

