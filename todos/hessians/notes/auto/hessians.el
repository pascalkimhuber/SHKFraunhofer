(TeX-add-style-hook
 "hessians"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8")))
   (TeX-run-style-hooks
    "latex2e"
    "scrartcl"
    "scrartcl10"
    "inputenc"
    "amsmath")
   (LaTeX-add-labels
    "sec:kraftberechnung"
    "sec:hesse-matrix"
    "sec:1.-fall"
    "sec:2.-fall"
    "sec:zusammenfassung"
    "sec:implementierung")))

