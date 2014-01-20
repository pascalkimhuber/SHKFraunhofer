(TeX-add-style-hook "todos"
 (lambda ()
    (LaTeX-add-labels
     "sec-1"
     "sec-1.1"
     "sec-1.2"
     "sec-2"
     "sec-3"
     "sec-4")
    (TeX-run-style-hooks
     "hyperref"
     "amssymb"
     "soul"
     "wrapfig"
     "float"
     "longtable"
     "graphicx"
     "fontenc"
     "T1"
     "inputenc"
     "utf8"
     "latex2e"
     "art11"
     "article"
     "11pt")))

