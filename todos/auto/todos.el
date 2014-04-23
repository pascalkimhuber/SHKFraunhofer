(TeX-add-style-hook "todos"
 (lambda ()
    (LaTeX-add-labels
     "sec-1"
     "sec-2"
     "sec-2-1"
     "sec-2-2"
     "sec-2-3"
     "sec-2-3-1"
     "sec-2-3-2")
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
     "art11"
     "article"
     "11pt")))

