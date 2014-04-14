(TeX-add-style-hook "overview"
 (lambda ()
    (LaTeX-add-labels
     "sec:basics"
     "sec:ensembles"
     "sec:micr-ensemble-nve"
     "sec:canon-ensemble-nvt"
     "sec:isoth-isob-ensemble"
     "sec:grand-canon-ensemble")
    (TeX-run-style-hooks
     "amssymb"
     "mathtools"
     "latex2e"
     "scrartcl10"
     "scrartcl"
     "")))

