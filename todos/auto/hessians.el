(TeX-add-style-hook "hessians"
 (lambda ()
    (LaTeX-add-labels
     "sec:kraftberechnung"
     "sec:hesse-matrix"
     "sec:1.-fall"
     "sec:2.-fall"
     "sec:zusammenfassung"
     "sec:implementierung")
    (TeX-run-style-hooks
     "amsmath"
     "inputenc"
     "utf8"
     "latex2e"
     "scrartcl10"
     "scrartcl"
     "")))

