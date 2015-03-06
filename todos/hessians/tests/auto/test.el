(TeX-add-style-hook
 "test"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("fontenc" "T1") ("ulem" "normalem")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art11"
    "inputenc"
    "fontenc"
    "fixltx2e"
    "graphicx"
    "longtable"
    "float"
    "wrapfig"
    "rotating"
    "ulem"
    "amsmath"
    "textcomp"
    "marvosym"
    "wasysym"
    "amssymb"
    "hyperref")
   (LaTeX-add-labels
    "sec-1"
    "sec-1-1"
    "eq:1"
    "eq:2"
    "sec-1-1-1"
    "sec-1-1-2"
    "sec-1-1-3"
    "sec-1-2"
    "sec-1-3"
    "sec-1-4"
    "sec-1-4-1"
    "sec-1-4-1-1"
    "sec-1-4-1-2"
    "sec-1-4-1-3"
    "sec-1-4-1-4"
    "sec-1-4-1-5"
    "sec-1-4-1-6"
    "sec-1-4-2"
    "sec-1-4-2-1"
    "sec-1-4-2-2"
    "sec-1-4-2-3"
    "sec-1-4-2-4"
    "sec-1-4-2-5"
    "sec-1-4-2-6"
    "sec-1-4-3"
    "sec-1-4-3-1"
    "sec-1-4-3-2"
    "sec-1-4-3-3"
    "sec-1-4-3-4"
    "sec-1-4-3-5"
    "sec-1-4-3-6"
    "sec-1-4-4"
    "sec-1-4-4-1"
    "sec-1-4-4-2"
    "sec-1-4-4-3"
    "sec-1-4-4-4"
    "sec-1-4-4-5"
    "sec-1-4-4-6"
    "sec-1-4-5"
    "sec-1-4-5-1"
    "sec-1-4-5-2"
    "sec-1-4-6"
    "sec-1-4-6-1"
    "sec-1-4-6-2"
    "sec-1-5"
    "sec-1-5-1"
    "sec-1-5-1-1")))

